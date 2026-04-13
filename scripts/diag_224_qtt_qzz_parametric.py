#!/usr/bin/env python3
"""
Run a small parametric diagnostic campaign over Q_tt blend and Q_zz scale
for HELMVEC3 Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_qtt_qzz_parametric.csv
- out/validation/DIAG_224_QTT_QZZ_PARAMETRIC.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_ROOT = ROOT / "out"
DEFAULT_BUILD_BIN = ROOT / "build" / "helmvec3_fig13_rect"
SECTION = "2.2.4"
CASE = "Figure13_Table10"
DEFAULT_BACKEND = "closed-form"
RUNS_SUBDIR = ".diag_224_qtt_qzz_parametric"
DELTA_TOL = 1.0e-12
DEFAULT_QTT_BLEND_VALUES = (0.0, 0.5, 1.0)
DEFAULT_QZZ_SCALE_VALUES = (1.0, 0.75, 0.5)


@dataclass(frozen=True)
class RepresentativePoint:
    key: str
    regime: str
    d_over_a: float
    br_over_lambda0: float

    @property
    def mode(self) -> str:
        return f"d/a={self.d_over_a},br/lambda0={self.br_over_lambda0}"


REPRESENTATIVE_POINTS = (
    RepresentativePoint("P1", "low", 0.167, 0.5),
    RepresentativePoint("P2", "low", 0.286, 0.5),
    RepresentativePoint("P3", "mid", 0.5, 0.4),
    RepresentativePoint("H1", "high", 0.5, 1.0),
    RepresentativePoint("H2", "high", 0.6, 1.0),
    RepresentativePoint("H3", "high", 0.8, 0.8),
    RepresentativePoint("H4", "high", 0.8, 1.0),
)


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


def _run_root(out_root: Path, label: str) -> Path:
    return out_root / RUNS_SUBDIR / label / DEFAULT_BACKEND.replace("-", "_")


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


def _load_table10(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    rows = _read_csv_rows(path)
    if len(rows) != 45:
        raise SystemExit(f"Expected 45 Table 10 rows in {_rel_from_root(path)}, found {len(rows)}.")
    lookup: dict[tuple[float, float], dict[str, str]] = {}
    for row in rows:
        key = _point_key(
            _to_float(row["d_over_a"], field="d_over_a", path=path),
            _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path),
        )
        if key in lookup:
            raise SystemExit(f"Duplicate Table 10 point {key} in {_rel_from_root(path)}.")
        lookup[key] = row
    return lookup


def _raw_summary(path: Path) -> dict[str, float | int]:
    rows = _read_csv_rows(path)
    useful = [row for row in rows if row["kept_after_dedup"] == "1"]
    if not useful:
        raise SystemExit(f"No useful spectrum rows in {_rel_from_root(path)}")

    def beta(row: dict[str, str]) -> float:
        return _to_float(row["beta_ratio_if_real_positive"], field="beta_ratio_if_real_positive", path=path)

    def ez_ratio(row: dict[str, str]) -> float:
        return _to_float(row["ez_ratio"], field="ez_ratio", path=path)

    ez_like = [row for row in useful if ez_ratio(row) >= 0.5]
    et_like = [row for row in useful if ez_ratio(row) < 0.5]
    highest_ez_like_beta = max((beta(row) for row in ez_like), default=float("nan"))
    lowest_et_like_beta = min((beta(row) for row in et_like), default=float("nan"))
    gap = (
        lowest_et_like_beta - highest_ez_like_beta if ez_like and et_like else float("nan")
    )
    return {
        "physical_kept_count": sum(row["kept_after_physical_filter"] == "1" for row in rows),
        "dedup_kept_count": len(useful),
        "highest_ez_like_beta": highest_ez_like_beta,
        "lowest_et_like_beta": lowest_et_like_beta,
        "ez_to_et_gap": gap,
    }


def _run_case(
    build_bin: Path,
    out_root: Path,
    label: str,
    variant: str,
    nx: int,
    ny: int,
    *,
    qtt_blend: float | None = None,
    qzz_scale: float | None = None,
) -> tuple[list[str], Path]:
    run_root = _run_root(out_root, label)
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(run_root)
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM"] = "1"
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM_TABLE10_ALL"] = "1"
    env["TP3485_HELMVEC3_DIAG_BETA_VARIANT"] = variant
    if qtt_blend is not None:
        env["TP3485_HELMVEC3_DIAG_QTT_BLEND"] = f"{qtt_blend:.16g}"
    else:
        env.pop("TP3485_HELMVEC3_DIAG_QTT_BLEND", None)
    if qzz_scale is not None:
        env["TP3485_HELMVEC3_DIAG_QZZ_SCALE"] = f"{qzz_scale:.16g}"
    else:
        env.pop("TP3485_HELMVEC3_DIAG_QZZ_SCALE", None)

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


def _load_run_data(
    build_bin: Path,
    out_root: Path,
    label: str,
    variant: str,
    nx: int,
    ny: int,
    *,
    qtt_blend: float | None = None,
    qzz_scale: float | None = None,
) -> tuple[dict[tuple[float, float], dict[str, str]], dict[tuple[float, float], dict[str, float | int]], list[str], Path]:
    cmd, run_root = _run_case(
        build_bin,
        out_root,
        label,
        variant,
        nx,
        ny,
        qtt_blend=qtt_blend,
        qzz_scale=qzz_scale,
    )
    table = _load_table10(_table10_path(run_root))
    raw = {
        _point_key(point.d_over_a, point.br_over_lambda0): _raw_summary(
            _raw_spectrum_path(run_root, point.d_over_a, point.br_over_lambda0)
        )
        for point in REPRESENTATIVE_POINTS
    }
    return table, raw, cmd, run_root


def _combo_summary(
    combo_table: dict[tuple[float, float], dict[str, str]],
    baseline_table: dict[tuple[float, float], dict[str, str]],
) -> dict[str, float | int]:
    selected_keys = {_point_key(p.d_over_a, p.br_over_lambda0) for p in REPRESENTATIVE_POINTS}

    def err_primary(row: dict[str, str]) -> float:
        return float(row["error_percent_analytic"])

    def err_secondary(row: dict[str, str]) -> float:
        return float(row["error_percent_helmvec3"])

    selected_rows = [combo_table[key] for key in selected_keys]
    all_rows = list(combo_table.values())

    selected_max_primary = max(err_primary(row) for row in selected_rows)
    selected_max_secondary = max(err_secondary(row) for row in selected_rows)
    full_max_primary = max(err_primary(row) for row in all_rows)
    full_max_secondary = max(err_secondary(row) for row in all_rows)

    improved_selected_primary = 0
    improved_selected_secondary = 0
    worsened_selected_primary = 0
    worsened_selected_secondary = 0
    improved_full_primary = 0
    improved_full_secondary = 0
    worsened_full_primary = 0
    worsened_full_secondary = 0

    for key, base in baseline_table.items():
        row = combo_table[key]
        delta_primary = err_primary(row) - err_primary(base)
        delta_secondary = err_secondary(row) - err_secondary(base)
        if key in selected_keys:
            if delta_primary < -DELTA_TOL:
                improved_selected_primary += 1
            elif delta_primary > DELTA_TOL:
                worsened_selected_primary += 1
            if delta_secondary < -DELTA_TOL:
                improved_selected_secondary += 1
            elif delta_secondary > DELTA_TOL:
                worsened_selected_secondary += 1
        if delta_primary < -DELTA_TOL:
            improved_full_primary += 1
        elif delta_primary > DELTA_TOL:
            worsened_full_primary += 1
        if delta_secondary < -DELTA_TOL:
            improved_full_secondary += 1
        elif delta_secondary > DELTA_TOL:
            worsened_full_secondary += 1

    return {
        "selected_max_primary_pct": selected_max_primary,
        "selected_max_secondary_pct": selected_max_secondary,
        "selected_score_max_pair": selected_max_primary + selected_max_secondary,
        "full_max_primary_pct": full_max_primary,
        "full_max_secondary_pct": full_max_secondary,
        "full_score_max_pair": full_max_primary + full_max_secondary,
        "selected_improved_primary_count": improved_selected_primary,
        "selected_improved_secondary_count": improved_selected_secondary,
        "selected_worsened_primary_count": worsened_selected_primary,
        "selected_worsened_secondary_count": worsened_selected_secondary,
        "full_improved_primary_count": improved_full_primary,
        "full_improved_secondary_count": improved_full_secondary,
        "full_worsened_primary_count": worsened_full_primary,
        "full_worsened_secondary_count": worsened_full_secondary,
    }


def _build_rows(
    baseline_table: dict[tuple[float, float], dict[str, str]],
    baseline_raw: dict[tuple[float, float], dict[str, float | int]],
    combo_table: dict[tuple[float, float], dict[str, str]],
    combo_raw: dict[tuple[float, float], dict[str, float | int]],
    combo_summary: dict[str, float | int],
    baseline_root: Path,
    combo_root: Path,
    qtt_blend: float,
    qzz_scale: float,
    nx: int,
    ny: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for point in REPRESENTATIVE_POINTS:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        base = baseline_table[key]
        cand = combo_table[key]
        base_raw = baseline_raw[key]
        cand_raw = combo_raw[key]
        rows.append(
            {
                "backend": DEFAULT_BACKEND,
                "section": SECTION,
                "case": CASE,
                "regime": point.regime,
                "priority_point": point.key,
                "mode": point.mode,
                "mesh_label": f"{nx}x{ny}",
                "nx": str(nx),
                "ny": str(ny),
                "qtt_blend": f"{qtt_blend:.16g}",
                "qzz_scale": f"{qzz_scale:.16g}",
                "ref_primary": base["beta_over_k0_analytic"],
                "ref_secondary": base["beta_over_k0_helmvec3"],
                "baseline_fem": base["beta_over_k0_fem_matched"],
                "variant_fem": cand["beta_over_k0_fem_matched"],
                "baseline_err_primary_pct": base["error_percent_analytic"],
                "variant_err_primary_pct": cand["error_percent_analytic"],
                "delta_err_primary_pct": f"{float(cand['error_percent_analytic']) - float(base['error_percent_analytic']):.16g}",
                "baseline_err_secondary_pct": base["error_percent_helmvec3"],
                "variant_err_secondary_pct": cand["error_percent_helmvec3"],
                "delta_err_secondary_pct": f"{float(cand['error_percent_helmvec3']) - float(base['error_percent_helmvec3']):.16g}",
                "baseline_ez_ratio": base["ez_ratio"],
                "variant_ez_ratio": cand["ez_ratio"],
                "baseline_selected_rank": base["selected_candidate_rank"],
                "variant_selected_rank": cand["selected_candidate_rank"],
                "baseline_selected_eig_index": base["selected_eig_index"],
                "variant_selected_eig_index": cand["selected_eig_index"],
                "selected_eig_changed": "1" if base["selected_eig_index"] != cand["selected_eig_index"] else "0",
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
                    f"{float(cand_raw['ez_to_et_gap']) - float(base_raw['ez_to_et_gap']):.16g}"
                    if math.isfinite(float(base_raw["ez_to_et_gap"])) and math.isfinite(float(cand_raw["ez_to_et_gap"]))
                    else "nan"
                ),
                "selected_score_max_pair": f"{float(combo_summary['selected_score_max_pair']):.16g}",
                "full_score_max_pair": f"{float(combo_summary['full_score_max_pair']):.16g}",
                "selected_max_primary_pct": f"{float(combo_summary['selected_max_primary_pct']):.16g}",
                "selected_max_secondary_pct": f"{float(combo_summary['selected_max_secondary_pct']):.16g}",
                "full_max_primary_pct": f"{float(combo_summary['full_max_primary_pct']):.16g}",
                "full_max_secondary_pct": f"{float(combo_summary['full_max_secondary_pct']):.16g}",
                "selected_improved_primary_count": str(combo_summary["selected_improved_primary_count"]),
                "selected_improved_secondary_count": str(combo_summary["selected_improved_secondary_count"]),
                "selected_worsened_primary_count": str(combo_summary["selected_worsened_primary_count"]),
                "selected_worsened_secondary_count": str(combo_summary["selected_worsened_secondary_count"]),
                "full_improved_primary_count": str(combo_summary["full_improved_primary_count"]),
                "full_improved_secondary_count": str(combo_summary["full_improved_secondary_count"]),
                "full_worsened_primary_count": str(combo_summary["full_worsened_primary_count"]),
                "full_worsened_secondary_count": str(combo_summary["full_worsened_secondary_count"]),
                "baseline_raw_spectrum_csv": _rel_from_root(_raw_spectrum_path(baseline_root, point.d_over_a, point.br_over_lambda0)),
                "variant_raw_spectrum_csv": _rel_from_root(_raw_spectrum_path(combo_root, point.d_over_a, point.br_over_lambda0)),
                "baseline_table10_csv": _rel_from_root(_table10_path(baseline_root)),
                "variant_table10_csv": _rel_from_root(_table10_path(combo_root)),
                "baseline_run_log": _rel_from_root(_run_log_path(baseline_root)),
                "variant_run_log": _rel_from_root(_run_log_path(combo_root)),
                "baseline_run_timing_csv": _rel_from_root(_run_timing_path(baseline_root)),
                "variant_run_timing_csv": _rel_from_root(_run_timing_path(combo_root)),
            }
        )
    return rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "regime",
        "priority_point",
        "mode",
        "mesh_label",
        "nx",
        "ny",
        "qtt_blend",
        "qzz_scale",
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
        "selected_score_max_pair",
        "full_score_max_pair",
        "selected_max_primary_pct",
        "selected_max_secondary_pct",
        "full_max_primary_pct",
        "full_max_secondary_pct",
        "selected_improved_primary_count",
        "selected_improved_secondary_count",
        "selected_worsened_primary_count",
        "selected_worsened_secondary_count",
        "full_improved_primary_count",
        "full_improved_secondary_count",
        "full_worsened_primary_count",
        "full_worsened_secondary_count",
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
) -> str:
    combo_rows: dict[tuple[float, float], dict[str, str]] = {}
    for row in rows:
        combo_rows.setdefault((float(row["qtt_blend"]), float(row["qzz_scale"])), row)

    ordered = sorted(
        combo_rows.values(),
        key=lambda row: (
            float(row["full_score_max_pair"]),
            float(row["selected_score_max_pair"]),
            float(row["qtt_blend"]),
            float(row["qzz_scale"]),
        ),
    )

    lines = [
        "# DIAG 224 Parametric Qtt/Qzz",
        "",
        f"Generated at {dt.datetime.now().isoformat(timespec='seconds')}.",
        "",
        "## Commands",
        "",
        f"- `baseline`: `{' '.join(shlex.quote(part) for part in baseline_cmd)}`",
        "- parametric variant: `TP3485_HELMVEC3_DIAG_BETA_VARIANT=diag_parametric_qtt_qzz "
        "TP3485_HELMVEC3_DIAG_QTT_BLEND=<w> TP3485_HELMVEC3_DIAG_QZZ_SCALE=<s> build/helmvec3_fig13_rect ...`",
        "",
        "## Combo Summary",
        "",
        "| qtt_blend | qzz_scale | selected_score | full_score | selected impA/impH | selected regA/regH | full impA/impH | full regA/regH | full max A | full max H |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    for row in ordered:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{float(row['qtt_blend']):.2f}",
                    f"{float(row['qzz_scale']):.2f}",
                    f"{float(row['selected_score_max_pair']):.3f}",
                    f"{float(row['full_score_max_pair']):.3f}",
                    f"{row['selected_improved_primary_count']}/{row['selected_improved_secondary_count']}",
                    f"{row['selected_worsened_primary_count']}/{row['selected_worsened_secondary_count']}",
                    f"{row['full_improved_primary_count']}/{row['full_improved_secondary_count']}",
                    f"{row['full_worsened_primary_count']}/{row['full_worsened_secondary_count']}",
                    f"{float(row['full_max_primary_pct']):.3f}",
                    f"{float(row['full_max_secondary_pct']):.3f}",
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Best Full-Score Combo Points",
            "",
            "| qtt_blend | qzz_scale | point | regime | delta A% | delta H% | baseline A% | variant A% | baseline H% | variant H% | delta gap | eig changed |",
            "| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )

    best_full = ordered[0]
    best_key = (float(best_full["qtt_blend"]), float(best_full["qzz_scale"]))
    for row in [r for r in rows if (float(r["qtt_blend"]), float(r["qzz_scale"])) == best_key]:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{float(row['qtt_blend']):.2f}",
                    f"{float(row['qzz_scale']):.2f}",
                    row["priority_point"],
                    row["regime"],
                    f"{float(row['delta_err_primary_pct']):.3f}",
                    f"{float(row['delta_err_secondary_pct']):.3f}",
                    f"{float(row['baseline_err_primary_pct']):.3f}",
                    f"{float(row['variant_err_primary_pct']):.3f}",
                    f"{float(row['baseline_err_secondary_pct']):.3f}",
                    f"{float(row['variant_err_secondary_pct']):.3f}",
                    f"{float(row['delta_ez_to_et_gap']):.3f}" if row["delta_ez_to_et_gap"] != "nan" else "nan",
                    row["selected_eig_changed"],
                ]
            )
            + " |"
        )

    return "\n".join(lines) + "\n"


def _parse_float_list(raw: str | None, fallback: tuple[float, ...]) -> tuple[float, ...]:
    if raw is None or raw.strip() == "":
        return fallback
    return tuple(float(part.strip()) for part in raw.split(",") if part.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=DEFAULT_OUT_ROOT / "validation" / "diag_224_qtt_qzz_parametric.csv",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=DEFAULT_OUT_ROOT / "validation" / "DIAG_224_QTT_QZZ_PARAMETRIC.md",
    )
    parser.add_argument("--build-bin", type=Path, default=DEFAULT_BUILD_BIN)
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=5)
    parser.add_argument("--qtt-blends", help="Comma-separated list, default 0,0.5,1.0")
    parser.add_argument("--qzz-scales", help="Comma-separated list, default 1.0,0.75,0.5")
    args = parser.parse_args()

    build_bin = _resolve(args.build_bin)
    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)
    qtt_blends = _parse_float_list(args.qtt_blends, DEFAULT_QTT_BLEND_VALUES)
    qzz_scales = _parse_float_list(args.qzz_scales, DEFAULT_QZZ_SCALE_VALUES)

    baseline_table, baseline_raw, baseline_cmd, baseline_root = _load_run_data(
        build_bin,
        out_root,
        "baseline",
        "baseline",
        args.nx,
        args.ny,
    )

    rows: list[dict[str, str]] = []
    for qtt_blend in qtt_blends:
        for qzz_scale in qzz_scales:
            label = f"qtt{_label_number(qtt_blend)}_qzz{_label_number(qzz_scale)}"
            combo_table, combo_raw, _cmd, combo_root = _load_run_data(
                build_bin,
                out_root,
                label,
                "diag_parametric_qtt_qzz",
                args.nx,
                args.ny,
                qtt_blend=qtt_blend,
                qzz_scale=qzz_scale,
            )
            summary = _combo_summary(combo_table, baseline_table)
            rows.extend(
                _build_rows(
                    baseline_table,
                    baseline_raw,
                    combo_table,
                    combo_raw,
                    summary,
                    baseline_root,
                    combo_root,
                    qtt_blend,
                    qzz_scale,
                    args.nx,
                    args.ny,
                )
            )

    rows.sort(
        key=lambda row: (
            float(row["qtt_blend"]),
            float(row["qzz_scale"]),
            row["priority_point"],
        )
    )
    _write_csv(out_csv, rows)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(_build_markdown(rows, baseline_cmd), encoding="utf-8")

    print(f"[diag_224_qtt_qzz_parametric] wrote {out_csv}")
    print(f"[diag_224_qtt_qzz_parametric] wrote {out_md}")


if __name__ == "__main__":
    main()
