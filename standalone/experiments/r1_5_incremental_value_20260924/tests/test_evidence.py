from pathlib import Path

import csv
import json

from evidence import build_sha256_manifest, seal_evidence_bundle, sha256, sha256_manifest_recheck_gate, verify_sha256_manifest


def test_sha256_manifest_can_be_rechecked(tmp_path: Path):
    (tmp_path / "a.txt").write_text("alpha", encoding="utf-8")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "b.txt").write_text("beta", encoding="utf-8")
    build_sha256_manifest(tmp_path)
    assert verify_sha256_manifest(tmp_path)
    (tmp_path / "a.txt").write_text("changed", encoding="utf-8")
    assert not verify_sha256_manifest(tmp_path)


def test_sealed_bundle_covers_the_final_gate(tmp_path: Path):
    (tmp_path / "result.csv").write_text("metric,value\npr_auc,0.5\n", encoding="utf-8")
    checks = [{"name": "numerical_gate", "passed": True, "details": "aligned"}]
    final_checks = seal_evidence_bundle(tmp_path, checks, "PASS")
    gate = json.loads((tmp_path / "gate.json").read_text(encoding="utf-8"))
    assert gate["status"] == "PASS"
    assert final_checks[-1]["name"] == "sha256_manifest_recheck"
    assert gate["checks"][-1]["passed"] is True
    assert verify_sha256_manifest(tmp_path)


def test_sealed_bundle_reuses_the_predeclared_final_manifest_gate(tmp_path: Path):
    (tmp_path / "result.csv").write_text("metric,value\npr_auc,0.5\n", encoding="utf-8")
    checks = [{"name": "numerical_gate", "passed": True, "details": "aligned"}, sha256_manifest_recheck_gate()]

    final_checks = seal_evidence_bundle(tmp_path, checks, "PASS")

    assert [check["name"] for check in final_checks].count("sha256_manifest_recheck") == 1
    assert verify_sha256_manifest(tmp_path)


def test_sealed_bundle_hashes_a_final_report_written_before_sealing(tmp_path: Path):
    (tmp_path / "RUN_REPORT.md").write_text("# Final report\n", encoding="utf-8")

    seal_evidence_bundle(tmp_path, [sha256_manifest_recheck_gate()], "SMOKE_PASS")

    with (tmp_path / "sha256_manifest.csv").open("r", newline="", encoding="utf-8-sig") as stream:
        hashes = {row["path"]: row["sha256"] for row in csv.DictReader(stream)}
    assert hashes["RUN_REPORT.md"] == sha256(tmp_path / "RUN_REPORT.md")
