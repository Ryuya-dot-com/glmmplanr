"""Read-only or fresh-directory verifier for saved-response replay CSVs.

Usage:
  # 1. Read-only check of archived replay outputs (no files modified):
  python3 validation/verify_replay_against_archive.py

  # 2. Verify a fresh recomputation directory against the archived baseline:
  python3 validation/verify_replay_against_archive.py /path/to/fresh_replay

  # 3. Optional explicit reference directory and JSON write:
  python3 validation/verify_replay_against_archive.py /path/to/fresh_replay validation/review_contract_replay --write
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

Z_975 = 1.959963984540054


def fail(message: str) -> None:
    print(f"VERIFICATION FAILED: {message}", file=sys.stderr)
    raise SystemExit(1)


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        fail(f"Missing CSV file: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        headers = list(reader.fieldnames or [])
        rows = list(reader)
    return headers, rows


def finite_diff_exceeds(a_str: str, b_str: str, tol: float) -> bool:
    if a_str == "NA" and b_str == "NA":
        return False
    try:
        a = float(a_str)
        b = float(b_str)
    except ValueError:
        return True
    if not (math.isfinite(a) and math.isfinite(b)):
        return True
    return abs(a - b) > tol


def max_finite_diff(
    tgt_map: dict[tuple[str, ...], dict[str, str]],
    ref_map: dict[tuple[str, ...], dict[str, str]],
    col: str,
) -> float:
    diffs = []
    for k in ref_map:
        a_str = tgt_map[k][col]
        b_str = ref_map[k][col]
        if a_str == "NA" and b_str == "NA":
            diffs.append(0.0)
        else:
            diffs.append(abs(float(a_str) - float(b_str)))
    return max(diffs) if diffs else 0.0


def verify_internal_fit_decisions(name: str, rows: list[dict[str, str]]) -> None:
    for idx, r in enumerate(rows, start=1):
        valid = r["valid"] == "TRUE"
        if not valid:
            continue
        est = float(r["estimate"])
        se = float(r["se"])
        if not (math.isfinite(est) and math.isfinite(se) and se > 0):
            fail(f"{name} row {idx}: valid=TRUE with non-positive/non-finite SE ({r['se']})")
        det = abs(est) > Z_975 * se
        wid = (2.0 * Z_975 * se) <= 20.0
        jnt = det and wid
        if (r["detected"] == "TRUE") != det:
            fail(f"{name} row {idx}: detected={r['detected']} inconsistent with estimate/se")
        if (r["width_ok"] == "TRUE") != wid:
            fail(f"{name} row {idx}: width_ok={r['width_ok']} inconsistent with se")
        if (r["joint"] == "TRUE") != jnt:
            fail(f"{name} row {idx}: joint={r['joint']} inconsistent with detected/width_ok")
        if "truth" in r:
            cov = abs(est - float(r["truth"])) <= Z_975 * se
            if (r["covered"] == "TRUE") != cov:
                fail(f"{name} row {idx}: covered={r['covered']} inconsistent with estimate/se/truth")
        if "logit" in r and "logit_se" in r and "logit_truth" in r and "logit_covered" in r:
            l_est = float(r["logit"])
            l_se = float(r["logit_se"])
            l_tr = float(r["logit_truth"])
            if not (math.isfinite(l_est) and math.isfinite(l_se) and l_se > 0 and math.isfinite(l_tr)):
                fail(f"{name} row {idx}: valid=TRUE with invalid logit/logit_se/logit_truth")
            l_cov = abs(l_est - l_tr) <= Z_975 * l_se
            if (r["logit_covered"] == "TRUE") != l_cov:
                fail(f"{name} row {idx}: logit_covered={r['logit_covered']} inconsistent with logit/logit_se/logit_truth")


def main(argv: list[str]) -> int:
    write_json = "--write" in argv
    pos = [a for a in argv if a != "--write"]
    val_root = Path(__file__).resolve().parent

    if len(pos) == 0:
        target = val_root / "review_contract_replay"
        reference = val_root / "curvature_contract_replay"
    elif len(pos) == 1:
        target = Path(pos[0]).resolve()
        reference = val_root / "review_contract_replay"
    elif len(pos) == 2:
        target = Path(pos[0]).resolve()
        reference = Path(pos[1]).resolve()
    else:
        fail("Usage: verify_replay_against_archive.py [TARGET_DIR [REFERENCE_DIR]] [--write]")

    specs = [
        (
            "comparison.csv",
            ("case", "replicate"),
            {
                "valid": None,
                "early_valid": None,
                "before_valid": None,
                "estimate": 1e-4,
                "se": 1e-4,
                "estimate_change_pp": 1e-4,
                "se_change_pp": 1e-4,
                "covariance_max_change": 1e-6,
                "covered": None,
                "detected": None,
                "width_ok": None,
                "joint": None,
                "coverage_changed": None,
                "detection_changed": None,
                "width_changed": None,
                "joint_changed": None,
                "response_hash": None,
            },
        ),
        (
            "reference_fits.csv",
            ("case", "replicate", "method"),
            {
                "layout": None,
                "chance": 1e-12,
                "correlation": None,
                "effect_pp": 1e-8,
                "n_subj": None,
                "n_item": None,
                "per_person": None,
                "valid": None,
                "estimate": 1e-4,
                "se": 1e-4,
                "truth": 1e-8,
                "covered": None,
                "detected": None,
                "width_ok": None,
                "joint": None,
                "logit": 1e-5,
                "logit_se": 1e-5,
                "logit_truth": 1e-8,
                "logit_covered": None,
                "response_hash": None,
            },
        ),
        (
            "reference_summary.csv",
            ("case", "replicate"),
            {
                "native_valid": None,
                "profile_valid": None,
                "early_valid": None,
                "native_reported_score": 1e-5,
                "native_logsd_score": 1e-5,
                "native_rho_score": 1e-5,
                "profile_score": 1e-5,
                "objective_improvement": 1e-5,
                "min_eigen": 1e-5,
                "boundary": None,
                "optimizer_code": None,
                "fine_se_difference_pp": 1e-4,
                "alternate_objective_difference": 1e-5,
                "evaluations": None,
                "failed_evaluations": None,
                "error": None,
                "before_profile_valid": None,
                "before_score_evaluable": None,
            },
        ),
    ]

    results: dict[str, object] = {
        "target_dir": str(target),
        "reference_dir": str(reference),
    }

    for name, keys, col_tols in specs:
        ref_headers, ref_rows = read_table(reference / name)
        tgt_headers, tgt_rows = read_table(target / name)
        if ref_headers != tgt_headers:
            fail(f"Header mismatch in {name}: target={tgt_headers} vs reference={ref_headers}")
        missing_spec = [c for c in col_tols if c not in ref_headers]
        if missing_spec:
            fail(f"Expected column(s) {missing_spec} missing from {name}")

        ref_map = {tuple(r[k] for k in keys): r for r in ref_rows}
        tgt_map = {tuple(r[k] for k in keys): r for r in tgt_rows}
        if not (len(ref_map) == len(ref_rows) == len(tgt_rows) == len(tgt_map)):
            fail(
                f"Row/key count mismatch in {name}: "
                f"ref_rows={len(ref_rows)}, ref_keys={len(ref_map)}, "
                f"tgt_rows={len(tgt_rows)}, tgt_keys={len(tgt_map)}"
            )
        if ref_map.keys() != tgt_map.keys():
            fail(f"Key mismatch in {name}")

        if name in ("comparison.csv", "reference_fits.csv"):
            verify_internal_fit_decisions(name, tgt_rows)

        changes: dict[str, int] = {}
        for col, tol in col_tols.items():
            if tol is not None:
                changes[col] = sum(
                    finite_diff_exceeds(tgt_map[k][col], ref_map[k][col], tol=tol)
                    for k in ref_map
                )
            else:
                changes[col] = sum(tgt_map[k][col] != ref_map[k][col] for k in ref_map)

        entry: dict[str, object] = {"rows": len(ref_rows), "changed_beyond_tolerance": changes}
        if "estimate" in col_tols:
            entry["max_point_difference_pp"] = max_finite_diff(tgt_map, ref_map, "estimate")
            entry["max_se_difference_pp"] = max_finite_diff(tgt_map, ref_map, "se")
        if "logit" in col_tols:
            entry["max_logit_point_difference"] = max_finite_diff(tgt_map, ref_map, "logit")
            entry["max_logit_se_difference"] = max_finite_diff(tgt_map, ref_map, "logit_se")
        if name == "reference_summary.csv":
            entry["failed_evaluations"] = sum(int(r["failed_evaluations"]) for r in tgt_rows)
        results[name] = entry

        if any(changes.values()):
            nonzero = {k: v for k, v in changes.items() if v > 0}
            fail(f"Column mismatch beyond tolerance in {name}: {json.dumps(nonzero)}")

    cpp_copy = target / "source_glmm_laplace.cpp"
    if cpp_copy.exists():
        results["source_sha256"] = hashlib.sha256(cpp_copy.read_bytes()).hexdigest()

    if write_json:
        (target / "verification.json").write_text(json.dumps(results, indent=2) + "\n")

    print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
