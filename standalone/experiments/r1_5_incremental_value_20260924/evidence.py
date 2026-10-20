from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: str | Path, value) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")


def collect_environment(cli_path: Path, source_files: list[Path], seed: int) -> dict:
    versions = {}
    for module_name in ("numpy", "pandas", "scipy", "sklearn", "matplotlib"):
        try:
            module = __import__(module_name)
            versions[module_name] = module.__version__
        except Exception as error:  # pragma: no cover - diagnostic only
            versions[module_name] = f"unavailable: {error}"
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "python_executable": sys.executable,
        "python_executable_sha256": sha256(sys.executable),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "packages": versions,
        "random_seed": seed,
        "cli_path": str(cli_path.resolve()),
        "cli_sha256": sha256(cli_path),
        "source_files": [{"path": str(path.resolve()), "sha256": sha256(path)} for path in source_files],
    }


def write_manifest(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_sha256_manifest(run_dir: Path) -> list[dict]:
    rows = []
    manifest_path = run_dir / "sha256_manifest.csv"
    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path != manifest_path:
            rows.append({"path": path.relative_to(run_dir).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_manifest(manifest_path, rows, ["path", "bytes", "sha256"])
    return rows


def verify_sha256_manifest(run_dir: Path) -> bool:
    manifest_path = run_dir / "sha256_manifest.csv"
    with manifest_path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return all((run_dir / row["path"]).is_file() and sha256(run_dir / row["path"]) == row["sha256"] for row in rows)


def initialize_run(output_root: Path) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_root / f"run_{timestamp}"
    suffix = 1
    while run_dir.exists():
        run_dir = output_root / f"run_{timestamp}_{suffix:02d}"
        suffix += 1
    for relative in ("configs", "logs", "maps/cpp", "maps/python_unweighted", "maps/python_weighted", "metrics", "models", "figures", "tables"):
        (run_dir / relative).mkdir(parents=True, exist_ok=True)
    return run_dir


def finalize_gate(
    run_dir: Path,
    checks: list[dict],
    status: str,
    error: str | None = None,
    metadata: dict | None = None,
) -> None:
    payload = {"status": status, "checks": checks, "error": error}
    if metadata:
        overlap = set(payload).intersection(metadata)
        if overlap:
            raise ValueError(f"Gate metadata may not overwrite reserved fields: {sorted(overlap)}")
        payload.update(metadata)
    write_json(run_dir / "gate.json", payload)


def sha256_manifest_recheck_gate() -> dict:
    """Return the final integrity gate shared by reports and bundle sealing."""
    return {
        "name": "sha256_manifest_recheck",
        "passed": True,
        "details": "All archived files match sha256_manifest.csv after final gate serialization.",
    }


def seal_evidence_bundle(
    run_dir: Path,
    checks: list[dict],
    status: str,
    error: str | None = None,
    metadata: dict | None = None,
) -> list[dict]:
    """Write the final gate and a verified manifest without leaving it stale.

    The manifest intentionally excludes itself.  The gate is written before
    each manifest build so its final contents are covered by the returned
    manifest.  The second verification is performed after all writes.
    """
    finalized = list(checks)
    finalize_gate(run_dir, finalized, status, error, metadata)
    build_sha256_manifest(run_dir)
    if not verify_sha256_manifest(run_dir):
        raise RuntimeError("Initial SHA-256 manifest verification failed")
    manifest_checks = [
        index for index, check in enumerate(finalized)
        if str(check.get("name", "")) == "sha256_manifest_recheck"
    ]
    if len(manifest_checks) > 1:
        raise ValueError("SHA-256 manifest recheck may be declared at most once")
    if manifest_checks:
        expected = sha256_manifest_recheck_gate()
        if finalized[manifest_checks[0]] != expected:
            raise ValueError("Predeclared SHA-256 manifest recheck does not match the sealing contract")
    else:
        finalized.append(sha256_manifest_recheck_gate())
    finalize_gate(run_dir, finalized, status, error, metadata)
    build_sha256_manifest(run_dir)
    if not verify_sha256_manifest(run_dir):
        raise RuntimeError("Final SHA-256 manifest verification failed")
    return finalized
