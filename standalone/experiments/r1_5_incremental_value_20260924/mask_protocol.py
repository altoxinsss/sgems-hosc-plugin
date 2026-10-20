"""Immutable bit-mask bundles used by the formal spatial-analysis protocol."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


MASK_ORDER = "C"
MASK_BITORDER = "little"


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_mask_bundle(
    base_path: str | Path,
    mask: np.ndarray,
    *,
    coordinate_basis: str,
    extra: dict[str, Any] | None = None,
) -> tuple[Path, Path, dict[str, Any]]:
    """Write a C-order/little-endian packed Boolean mask and its metadata."""
    base = Path(base_path)
    values = np.asarray(mask, dtype=bool)
    if values.ndim != 3:
        raise ValueError("Formal spatial masks must have (z, y, x) shape")
    flat = values.reshape(-1, order=MASK_ORDER)
    packed = np.packbits(flat, bitorder=MASK_BITORDER)
    npz_path = base.with_suffix(".npz")
    json_path = base.with_suffix(".json")
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, packed=packed)
    metadata: dict[str, Any] = {
        "shape_zyx": [int(value) for value in values.shape],
        "flatten_order": MASK_ORDER,
        "bitorder": MASK_BITORDER,
        "valid_bits": int(flat.size),
        "valid_count": int(flat.sum()),
        "coordinate_basis": coordinate_basis,
        "dtype": "bool",
        "npz_file": npz_path.name,
        "npz_sha256": sha256(npz_path),
    }
    if extra:
        overlap = set(metadata).intersection(extra)
        if overlap:
            raise ValueError(f"Mask metadata cannot overwrite protocol fields: {sorted(overlap)}")
        metadata.update(extra)
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return npz_path, json_path, metadata


def read_mask_bundle(base_path: str | Path) -> tuple[np.ndarray, dict[str, Any]]:
    """Read, verify and exactly reconstruct a formal bit-mask bundle."""
    base = Path(base_path)
    npz_path = base.with_suffix(".npz")
    json_path = base.with_suffix(".json")
    metadata = json.loads(json_path.read_text(encoding="utf-8"))
    if metadata.get("flatten_order") != MASK_ORDER or metadata.get("bitorder") != MASK_BITORDER:
        raise ValueError("Mask bundle does not use the frozen C/little protocol")
    if metadata.get("npz_sha256") != sha256(npz_path):
        raise ValueError("Mask bundle SHA-256 does not match metadata")
    shape = tuple(int(value) for value in metadata.get("shape_zyx", []))
    if len(shape) != 3 or any(value <= 0 for value in shape):
        raise ValueError("Mask bundle has an invalid (z, y, x) shape")
    bit_count = int(metadata.get("valid_bits", -1))
    if bit_count != int(np.prod(shape)):
        raise ValueError("Mask bundle bit count does not match shape")
    with np.load(npz_path) as archive:
        if "packed" not in archive:
            raise ValueError("Mask bundle has no packed payload")
        unpacked = np.unpackbits(np.asarray(archive["packed"], dtype=np.uint8), bitorder=MASK_BITORDER)
    if len(unpacked) < bit_count:
        raise ValueError("Mask bundle payload is shorter than the declared bit count")
    result = unpacked[:bit_count].astype(bool).reshape(shape, order=MASK_ORDER)
    if int(result.sum()) != int(metadata.get("valid_count", -1)):
        raise ValueError("Mask bundle valid-count metadata does not match payload")
    return result, metadata


def write_packed_boolean_bundle(
    base_path: str | Path,
    values: np.ndarray,
    *,
    coordinate_basis: str,
    extra: dict[str, Any] | None = None,
) -> tuple[Path, Path, dict[str, Any]]:
    """Archive an arbitrary Boolean array with the frozen bit-packing protocol.

    Spatial common-support masks use :func:`write_mask_bundle`, whose 3-D
    shape is deliberate.  Matched-ESS and same-weight permutation controls
    additionally need to archive one Boolean selection per window/sample
    pair.  This helper uses exactly the same C-order/little-endian contract
    without pretending that a two-dimensional window-by-sample array is a
    physical volume.
    """
    base = Path(base_path)
    mask = np.asarray(values, dtype=bool)
    if mask.ndim < 1:
        raise ValueError("Packed Boolean bundles require at least one dimension")
    flat = mask.reshape(-1, order=MASK_ORDER)
    packed = np.packbits(flat, bitorder=MASK_BITORDER)
    npz_path = base.with_suffix(".npz")
    json_path = base.with_suffix(".json")
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, packed=packed)
    metadata: dict[str, Any] = {
        "shape": [int(value) for value in mask.shape],
        "flatten_order": MASK_ORDER,
        "bitorder": MASK_BITORDER,
        "valid_bits": int(flat.size),
        "true_count": int(flat.sum()),
        "coordinate_basis": coordinate_basis,
        "dtype": "bool",
        "npz_file": npz_path.name,
        "npz_sha256": sha256(npz_path),
    }
    if extra:
        overlap = set(metadata).intersection(extra)
        if overlap:
            raise ValueError(f"Packed Boolean metadata cannot overwrite protocol fields: {sorted(overlap)}")
        metadata.update(extra)
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return npz_path, json_path, metadata


def read_packed_boolean_bundle(base_path: str | Path) -> tuple[np.ndarray, dict[str, Any]]:
    """Read and hash-verify a bundle written by :func:`write_packed_boolean_bundle`."""
    base = Path(base_path)
    npz_path = base.with_suffix(".npz")
    json_path = base.with_suffix(".json")
    metadata = json.loads(json_path.read_text(encoding="utf-8"))
    if metadata.get("flatten_order") != MASK_ORDER or metadata.get("bitorder") != MASK_BITORDER:
        raise ValueError("Packed Boolean bundle does not use the frozen C/little protocol")
    if metadata.get("npz_sha256") != sha256(npz_path):
        raise ValueError("Packed Boolean bundle SHA-256 does not match metadata")
    shape = tuple(int(value) for value in metadata.get("shape", []))
    if not shape or any(value <= 0 for value in shape):
        raise ValueError("Packed Boolean bundle has an invalid shape")
    bit_count = int(metadata.get("valid_bits", -1))
    if bit_count != int(np.prod(shape)):
        raise ValueError("Packed Boolean bundle bit count does not match shape")
    with np.load(npz_path) as archive:
        if "packed" not in archive:
            raise ValueError("Packed Boolean bundle has no packed payload")
        unpacked = np.unpackbits(np.asarray(archive["packed"], dtype=np.uint8), bitorder=MASK_BITORDER)
    if len(unpacked) < bit_count:
        raise ValueError("Packed Boolean bundle payload is shorter than the declared bit count")
    result = unpacked[:bit_count].astype(bool).reshape(shape, order=MASK_ORDER)
    if int(result.sum()) != int(metadata.get("true_count", -1)):
        raise ValueError("Packed Boolean bundle true-count metadata does not match payload")
    return result, metadata
