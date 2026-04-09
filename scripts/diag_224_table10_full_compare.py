#!/usr/bin/env python3
"""
Compare baseline vs a diagnostic beta-variant across the full Table 10.

Default outputs:
- out/validation/diag_224_table10_full_compare.csv
- out/validation/DIAG_224_TABLE10_FULL_COMPARE.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import shlex
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_ROOT = ROOT / "out"
DEFAULT_BUILD_BIN = ROOT / "build" / "helmvec3_fig13_rect"
SECTION = "2.2.4"
CASE = "Figure13_Table10"
DEFAULT_BACKEND = "closed-form"
DEFAULT_VARIANT = "diag_eq141_eps_mass_qtt_plus_qzz_half"
RUNS_SUBDIR = ".diag_224_table10_full_compare"
DELTA_TOL = 1.0e-12
EXPECTED_BLOCK_COUNTS = {
    0.0: 6,
    0.167: 6,
    0.286: 6,
    0.375: 6,
    0.5: 7,
    0.6: 7,
    0.8: 7,
}


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _label_number(value: float) -> str:
    return f"{value:.12g}".replace(".", "_")


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _run_root(out_root: Path, variant: str) -> Path:
    return out_root / RUNS_SUBDIR / variant / DEFAULT_BACKEND.replace("-", "_")


def _table10_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "csv" / "helmvec3_fig13_rect_table10.csv"


def _raw_spectrum_path(out_root: Path, d_over_a: float, br_over_lambda0: float) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(d_over_a)
        + "_br"
        + _label_number(br_over_lambda0)
        + "_raw_spectrum.csv"
    )
    return out_root / "helmvec3" / "fig13_rect" / "csv" / filename


def _run_log_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "run.log"


def _run_timing_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "run_timing.csv"


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    return rows


def _to_float(raw: str, *, field: str, path: Path) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError) as exc:
        raise SystemExit(
            f"Invalid numeric value for {field!r} in {_rel_from_root(path)}: {raw!r}"
        ) from exc


def _load_table10_lookup(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    rows = _read_csv_rows(path)
    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    block_counts: dict[float, int] = {}
    for row in rows:
        d_over_a = _to_float(row["d_over_a"], field="d_over_a", path=path)
        br_over_lambda0 = _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path)
        key = _point_key(d_over_a, br_over_lambda0)
        grouped.setdefault(key, []).append(row)
        block_counts[round(d_over_a, 12)] = block_counts.get(round(d_over_a, 12), 0) + 1

    if len(grouped) != 45:
        raise SystemExit(
            f"Expected 45 Table 10 rows in {_rel_from_root(path)}, found {len(grouped)}."
        )

    for d_over_a, expected in EXPECTED_BLOCK_COUNTS.items():
        observed = block_counts.get(round(d_over_a, 12), 0)
        if observed != expected:
            raise SystemExit(
                f"Unexpected Table 10 block count for d/a={d_over_a} in {_rel_from_root(path)}: "
                f"expected {expected}, found {observed}."
            )

    out: dict[tuple[float, float], dict[str, str]] = {}
    for key, matches in grouped.items():
        if len(matches) != 1:
            raise SystemExit(
                f"Expected exactly one Table 10 row for {key} in {_rel_from_root(path)}, found {len(matches)}."
            )
        out[key] = matches[0]
    return out


def _run_variant(
    build_bin: Path,
    out_root: Path,
    variant: str,
    nx: int,
    ny: int,
) -> tuple[list[str], Path]:
    run_root = _run_root(out_root, variant)
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(run_root)
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM"] = "1"
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM_TABLE10_ALL"] = "1"
    env["TP3485_HELMVEC3_DIAG_BETA_VARIANT"] = variant
    cmd = [
        str(build_bin),
        "--d-over-a-preview",
        "0.20",
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--backend",
        DEFAULT_BACKEND,
    ]
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    return cmd, run_root


def _raw_summary(raw_path: Path) -> dict[str, float | int]:
    raw_rows = _read_csv_rows(raw_path)
    useful = [row for row in raw_rows if row["kept_after_dedup"] == "1"]
    if not useful:
        raise SystemExit(f"No useful spectrum rows in {_rel_from_root(raw_path)}")

    def beta(row: dict[str, str]) -> float:
        return _to_float(row["beta_ratio_if_real_positive"], field="beta_ratio_if_real_positive", path=raw_path)

    def ez_ratio(row: dict[str, str]) -> float:
        return _to_float(row["ez_ratio"], field="ez_ratio", path=raw_path)

    ez_like = [row for row in useful if ez_ratio(row) >= 0.5]
    et_like = [row for row in useful if ez_ratio(row) < 0.5]
    highest_ez = max((beta(row) for row in ez_like), default=float("nan"))
    lowest_et = min((beta(row) for row in et_like), default=float("nan"))
    gap = lowest_et - highest_ez if ez_like and et_like else float("nan")

    selected_rows = [row for row in raw_rows if row["selected_after_matching"] == "1"]
    if len(selected_rows) != 1:
        raise SystemExit(
            f"Expected exactly one selected raw-spectrum row in {_rel_from_root(raw_path)}, found {len(selected_rows)}."
        )
    selected = selected_rows[0]

    return {
        "physical_kept_count": sum(row["kept_after_physical_filter"] == "1" for row in raw_rows),
        "dedup_kept_count": len(useful),
        "highest_ez_like_beta": highest_ez,
        "lowest_et_like_beta": lowest_et,
        "ez_to_et_gap": gap,
        "selected_solver_index": selected["solver_index"],
        "selected_candidate_rank_after_dedup": selected["candidate_rank_after_dedup"],
    }


def _variant_data(
    build_bin: Path,
    out_root: Path,
    variant: str,
    nx: int,
    ny: int,
) -> tuple[dict[tuple[float, float], dict[str, str]], dict[tuple[float, float], dict[str, float | int]], list[str], Path]:
    cmd, run_root = _run_variant(build_bin, out_root, variant, nx, ny)
    table10_path = _table10_path(run_root)
    table_lookup = _load_table10_lookup(table10_path)
    raw_lookup: dict[tuple[float, float], dict[str, float | int]] = {}
    for key in table_lookup:
        d_over_a, br_over_lambda0 = key
        raw_lookup[key] = _raw_summary(_raw_spectrum_path(run_root, d_over_a, br_over_lambda0))
    return table_lookup, raw_lookup, cmd, run_root


def _fmt_intish(value: str) -> int:
    return int(value)


def _build_rows(
    baseline_table: dict[tuple[float, float], dict[str, str]],
    baseline_raw: dict[tuple[float, float], dict[str, float | int]],
    variant_table: dict[tuple[float, float], dict[str, str]],
    variant_raw: dict[tuple[float, float], dict[str, float | int]],
    baseline_root: Path,
    variant_root: Path,
    variant_name: str,
    nx: int,
    ny: int,
) -> list[dict[str, str]]:
    keys = sorted(baseline_table.keys())
    if keys != sorted(variant_table.keys()):
        raise SystemExit("Baseline and variant Table 10 points do not match.")

    rows: list[dict[str, str]] = []
    for key in keys:
        d_over_a, br_over_lambda0 = key
        base = baseline_table[key]
        cand = variant_table[key]
        base_raw = baseline_raw[key]
        cand_raw = variant_raw[key]
        if (
            base["beta_over_k0_analytic"] != cand["beta_over_k0_analytic"]
            or base["beta_over_k0_helmvec3"] != cand["beta_over_k0_helmvec3"]
        ):
            raise SystemExit(f"Reference mismatch between baseline and variant for point {key}.")

        base_err_primary = _to_float(base["error_percent_analytic"], field="error_percent_analytic", path=_table10_path(baseline_root))
        cand_err_primary = _to_float(cand["error_percent_analytic"], field="error_percent_analytic", path=_table10_path(variant_root))
        base_err_secondary = _to_float(base["error_percent_helmvec3"], field="error_percent_helmvec3", path=_table10_path(baseline_root))
        cand_err_secondary = _to_float(cand["error_percent_helmvec3"], field="error_percent_helmvec3", path=_table10_path(variant_root))
        base_gap = float(base_raw["ez_to_et_gap"])
        cand_gap = float(cand_raw["ez_to_et_gap"])

        rows.append(
            {
                "backend": DEFAULT_BACKEND,
                "section": SECTION,
                "case": CASE,
                "variant": variant_name,
                "mode": f"d/a={d_over_a},br/lambda0={br_over_lambda0}",
                "mesh_label": f"{nx}x{ny}",
                "nx": str(nx),
                "ny": str(ny),
                "d_over_a": f"{d_over_a:.12g}",
                "br_over_lambda0": f"{br_over_lambda0:.12g}",
                "ref_primary": base["beta_over_k0_analytic"],
                "ref_secondary": base["beta_over_k0_helmvec3"],
                "baseline_fem": base["beta_over_k0_fem_matched"],
                "variant_fem": cand["beta_over_k0_fem_matched"],
                "baseline_err_primary_pct": base["error_percent_analytic"],
                "variant_err_primary_pct": cand["error_percent_analytic"],
                "delta_err_primary_pct": f"{cand_err_primary - base_err_primary:.16g}",
                "baseline_err_secondary_pct": base["error_percent_helmvec3"],
                "variant_err_secondary_pct": cand["error_percent_helmvec3"],
                "delta_err_secondary_pct": f"{cand_err_secondary - base_err_secondary:.16g}",
                "baseline_ez_ratio": base["ez_ratio"],
                "variant_ez_ratio": cand["ez_ratio"],
                "baseline_selected_rank": base["selected_candidate_rank"],
                "variant_selected_rank": cand["selected_candidate_rank"],
                "baseline_selected_eig_index": base["selected_eig_index"],
                "variant_selected_eig_index": cand["selected_eig_index"],
                "selected_eig_changed": "1" if base["selected_eig_index"] != cand["selected_eig_index"] else "0",
                "baseline_match_status": base["match_status"],
                "variant_match_status": cand["match_status"],
                "baseline_field_status": base["field_status"],
                "variant_field_status": cand["field_status"],
                "baseline_physical_kept_count": str(base_raw["physical_kept_count"]),
                "variant_physical_kept_count": str(cand_raw["physical_kept_count"]),
                "baseline_dedup_kept_count": str(base_raw["dedup_kept_count"]),
                "variant_dedup_kept_count": str(cand_raw["dedup_kept_count"]),
                "baseline_highest_ez_like_beta": str(base_raw["highest_ez_like_beta"]),
                "variant_highest_ez_like_beta": str(cand_raw["highest_ez_like_beta"]),
                "baseline_lowest_et_like_beta": str(base_raw["lowest_et_like_beta"]),
                "variant_lowest_et_like_beta": str(cand_raw["lowest_et_like_beta"]),
                "baseline_ez_to_et_gap": str(base_raw["ez_to_et_gap"]),
                "variant_ez_to_et_gap": str(cand_raw["ez_to_et_gap"]),
                "delta_ez_to_et_gap": (
                    f"{cand_gap - base_gap:.16g}" if math.isfinite(base_gap) and math.isfinite(cand_gap) else "nan"
                ),
                "improved_primary": "1" if cand_err_primary < (base_err_primary - DELTA_TOL) else "0",
                "improved_secondary": "1" if cand_err_secondary < (base_err_secondary - DELTA_TOL) else "0",
                "baseline_raw_spectrum_csv": _rel_from_root(_raw_spectrum_path(baseline_root, d_over_a, br_over_lambda0)),
                "variant_raw_spectrum_csv": _rel_from_root(_raw_spectrum_path(variant_root, d_over_a, br_over_lambda0)),
                "baseline_table10_csv": _rel_from_root(_table10_path(baseline_root)),
                "variant_table10_csv": _rel_from_root(_table10_path(variant_root)),
                "baseline_run_log": _rel_from_root(_run_log_path(baseline_root)),
                "variant_run_log": _rel_from_root(_run_log_path(variant_root)),
                "baseline_run_timing_csv": _rel_from_root(_run_timing_path(baseline_root)),
                "variant_run_timing_csv": _rel_from_root(_run_timing_path(variant_root)),
            }
        )
    return rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "variant",
        "mode",
        "mesh_label",
        "nx",
        "ny",
        "d_over_a",
        "br_over_lambda0",
        "ref_primary",
        "ref_secondary",
        "baseline_fem",
        "variant_fem",
        "baseline_err_primary_pct",
        "variant_err_primary_pct",
        "delta_err_primary_pct",
        "baseline_err_secondary_pct",
        "variant_err_secondary_pct",
        "delta_err_secondary_pct",
        "baseline_ez_ratio",
        "variant_ez_ratio",
        "baseline_selected_rank",
        "variant_selected_rank",
        "baseline_selected_eig_index",
        "variant_selected_eig_index",
        "selected_eig_changed",
        "baseline_match_status",
        "variant_match_status",
        "baseline_field_status",
        "variant_field_status",
        "baseline_physical_kept_count",
        "variant_physical_kept_count",
        "baseline_dedup_kept_count",
        "variant_dedup_kept_count",
        "baseline_highest_ez_like_beta",
        "variant_highest_ez_like_beta",
        "baseline_lowest_et_like_beta",
        "variant_lowest_et_like_beta",
        "baseline_ez_to_et_gap",
        "variant_ez_to_et_gap",
        "delta_ez_to_et_gap",
        "improved_primary",
        "improved_secondary",
        "baseline_raw_spectrum_csv",
        "variant_raw_spectrum_csv",
        "baseline_table10_csv",
        "variant_table10_csv",
        "baseline_run_log",
        "variant_run_log",
        "baseline_run_timing_csv",
        "variant_run_timing_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_markdown(
    rows: list[dict[str, str]],
    baseline_cmd: list[str],
    variant_cmd: list[str],
    variant_name: str,
) -> str:
    improved_primary = sum(row["improved_primary"] == "1" for row in rows)
    improved_secondary = sum(row["improved_secondary"] == "1" for row in rows)
    changed_eig = sum(row["selected_eig_changed"] == "1" for row in rows)
    worsened_primary = sum(float(row["delta_err_primary_pct"]) > DELTA_TOL for row in rows)
    worsened_secondary = sum(float(row["delta_err_secondary_pct"]) > DELTA_TOL for row in rows)
    max_base_primary = max(float(row["baseline_err_primary_pct"]) for row in rows)
    max_var_primary = max(float(row["variant_err_primary_pct"]) for row in rows)
    max_base_secondary = max(float(row["baseline_err_secondary_pct"]) for row in rows)
    max_var_secondary = max(float(row["variant_err_secondary_pct"]) for row in rows)

    top_rows = sorted(rows, key=lambda row: float(row["variant_err_primary_pct"]), reverse=True)[:10]
    regressions = [row for row in rows if float(row["delta_err_primary_pct"]) > DELTA_TOL][:5]

    lines = [
        "# DIAG 224 Table 10 Full Compare",
        "",
        f"Generated at {dt.datetime.now().isoformat(timespec='seconds')}.",
        "",
        "## Commands",
        "",
        f"- `baseline`: `{' '.join(shlex.quote(part) for part in baseline_cmd)}`",
        f"- `{variant_name}`: `{' '.join(shlex.quote(part) for part in variant_cmd)}`",
        "",
        "## Summary",
        "",
        f"- total_points: `{len(rows)}`",
        f"- improved_primary: `{improved_primary}`",
        f"- improved_secondary: `{improved_secondary}`",
        f"- changed_selected_eig_index: `{changed_eig}`",
        f"- worsened_primary: `{worsened_primary}`",
        f"- worsened_secondary: `{worsened_secondary}`",
        f"- max_primary_err: `{max_base_primary:.3f}% -> {max_var_primary:.3f}%`",
        f"- max_secondary_err: `{max_base_secondary:.3f}% -> {max_var_secondary:.3f}%`",
        "",
        "## Worst Variant Errors",
        "",
        "| Mode | Baseline A% | Variant A% | Delta A% | Baseline H% | Variant H% | Delta gap Ez->Et | eig changed |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in top_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["mode"],
                    f"{float(row['baseline_err_primary_pct']):.3f}",
                    f"{float(row['variant_err_primary_pct']):.3f}",
                    f"{float(row['delta_err_primary_pct']):.3f}",
                    f"{float(row['baseline_err_secondary_pct']):.3f}",
                    f"{float(row['variant_err_secondary_pct']):.3f}",
                    f"{float(row['delta_ez_to_et_gap']):.3f}",
                    row["selected_eig_changed"],
                ]
            )
            + " |"
        )

    if regressions:
        lines.extend(
            [
                "",
                "## Primary Regressions",
                "",
                "| Mode | Baseline A% | Variant A% | Delta A% |",
                "| --- | ---: | ---: | ---: |",
            ]
        )
        for row in regressions:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["mode"],
                        f"{float(row['baseline_err_primary_pct']):.3f}",
                        f"{float(row['variant_err_primary_pct']):.3f}",
                        f"{float(row['delta_err_primary_pct']):.3f}",
                    ]
                )
                + " |"
            )

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=DEFAULT_OUT_ROOT / "validation" / "diag_224_table10_full_compare.csv",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=DEFAULT_OUT_ROOT / "validation" / "DIAG_224_TABLE10_FULL_COMPARE.md",
    )
    parser.add_argument("--build-bin", type=Path, default=DEFAULT_BUILD_BIN)
    parser.add_argument("--variant", default=DEFAULT_VARIANT)
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=5)
    args = parser.parse_args()

    build_bin = _resolve(args.build_bin)
    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)

    baseline_table, baseline_raw, baseline_cmd, baseline_root = _variant_data(
        build_bin, out_root, "baseline", args.nx, args.ny
    )
    variant_table, variant_raw, variant_cmd, variant_root = _variant_data(
        build_bin, out_root, args.variant, args.nx, args.ny
    )
    rows = _build_rows(
        baseline_table,
        baseline_raw,
        variant_table,
        variant_raw,
        baseline_root,
        variant_root,
        args.variant,
        args.nx,
        args.ny,
    )
    _write_csv(out_csv, rows)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(_build_markdown(rows, baseline_cmd, variant_cmd, args.variant), encoding="utf-8")

    print(f"[diag_224_table10_full_compare] wrote {out_csv}")
    print(f"[diag_224_table10_full_compare] wrote {out_md}")


if __name__ == "__main__":
    main()
