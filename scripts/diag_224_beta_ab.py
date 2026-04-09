#!/usr/bin/env python3
"""
Run a controlled A/B diagnostic campaign for HELMVEC3 Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_beta_ab.csv
- out/validation/DIAG_224_BETA_AB.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
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
DEFAULT_BACKENDS = ("closed-form", "gauss")
DEFAULT_VARIANTS = ("baseline", "diag_eq141_eps_mass_qtt", "diag_eq142_doc_qzz")
RUNS_SUBDIR = ".diag_224_beta_ab"


@dataclass(frozen=True)
class PriorityPoint:
    key: str
    d_over_a: float
    br_over_lambda0: float

    @property
    def mode(self) -> str:
        return f"d/a={self.d_over_a},br/lambda0={self.br_over_lambda0}"


PRIORITY_POINTS = (
    PriorityPoint("P1", 0.167, 0.5),
    PriorityPoint("P2", 0.286, 0.5),
    PriorityPoint("P3", 0.5, 0.4),
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


def _backend_tag(backend: str) -> str:
    return backend.replace("-", "_")


def _run_root(out_root: Path, variant: str, backend: str) -> Path:
    return out_root / RUNS_SUBDIR / variant / _backend_tag(backend)


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _table10_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "csv" / "helmvec3_fig13_rect_table10.csv"


def _raw_spectrum_path(out_root: Path, point: PriorityPoint) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(point.d_over_a)
        + "_br"
        + _label_number(point.br_over_lambda0)
        + "_raw_spectrum.csv"
    )
    return out_root / "helmvec3" / "fig13_rect" / "csv" / filename


def _matrix_audit_path(out_root: Path, point: PriorityPoint) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(point.d_over_a)
        + "_br"
        + _label_number(point.br_over_lambda0)
        + "_matrix_audit.csv"
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


def _to_int(raw: str, *, field: str, path: Path) -> int:
    try:
        return int(raw)
    except (TypeError, ValueError) as exc:
        raise SystemExit(
            f"Invalid integer value for {field!r} in {_rel_from_root(path)}: {raw!r}"
        ) from exc


def _load_table10_lookup(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    rows = _read_csv_rows(path)
    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in rows:
        key = _point_key(
            _to_float(row["d_over_a"], field="d_over_a", path=path),
            _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path),
        )
        grouped.setdefault(key, []).append(row)

    out: dict[tuple[float, float], dict[str, str]] = {}
    for point in PRIORITY_POINTS:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        matches = grouped.get(key, [])
        if len(matches) != 1:
            raise SystemExit(
                f"Expected exactly one Table 10 row for {point.mode} in {_rel_from_root(path)}, found {len(matches)}."
            )
        out[key] = matches[0]
    return out


def _run_variant(
    build_bin: Path,
    out_root: Path,
    backend: str,
    variant: str,
    nx: int,
    ny: int,
) -> tuple[list[str], Path]:
    run_root = _run_root(out_root, variant, backend)
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(run_root)
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM"] = "1"
    env["TP3485_HELMVEC3_DIAG_MATRIX_AUDIT"] = "1"
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
        backend,
    ]
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    return cmd, run_root


def _best_useful_candidates(
    raw_rows: list[dict[str, str]],
    table_row: dict[str, str],
    raw_path: Path,
) -> dict[str, float | int]:
    useful = [row for row in raw_rows if row["kept_after_dedup"] == "1"]
    if not useful:
        raise SystemExit(f"No useful spectrum rows in {_rel_from_root(raw_path)}")

    analytic = _to_float(table_row["beta_over_k0_analytic"], field="beta_over_k0_analytic", path=raw_path)
    helmvec3_ref = _to_float(table_row["beta_over_k0_helmvec3"], field="beta_over_k0_helmvec3", path=raw_path)

    def beta(row: dict[str, str]) -> float:
        return _to_float(row["beta_ratio_if_real_positive"], field="beta_ratio_if_real_positive", path=raw_path)

    def ez_ratio(row: dict[str, str]) -> float:
        return _to_float(row["ez_ratio"], field="ez_ratio", path=raw_path)

    closest_analytic = min(useful, key=lambda row: abs(beta(row) - analytic))
    closest_helmvec3 = min(useful, key=lambda row: abs(beta(row) - helmvec3_ref))
    ez_like = [row for row in useful if ez_ratio(row) >= 0.5]
    et_like = [row for row in useful if ez_ratio(row) < 0.5]

    out: dict[str, float | int] = {
        "physical_kept_count": sum(row["kept_after_physical_filter"] == "1" for row in raw_rows),
        "dedup_kept_count": len(useful),
        "closest_useful_beta_to_analytic": beta(closest_analytic),
        "closest_useful_err_abs_analytic": abs(beta(closest_analytic) - analytic),
        "closest_useful_rank_to_analytic": _to_int(
            closest_analytic["candidate_rank_after_dedup"], field="candidate_rank_after_dedup", path=raw_path
        ),
        "closest_useful_beta_to_helmvec3": beta(closest_helmvec3),
        "closest_useful_err_abs_helmvec3": abs(beta(closest_helmvec3) - helmvec3_ref),
        "closest_useful_rank_to_helmvec3": _to_int(
            closest_helmvec3["candidate_rank_after_dedup"], field="candidate_rank_after_dedup", path=raw_path
        ),
        "highest_ez_like_beta": max((beta(row) for row in ez_like), default=float("nan")),
        "lowest_et_like_beta": min((beta(row) for row in et_like), default=float("nan")),
    }
    if ez_like and et_like:
        out["ez_to_et_gap"] = out["lowest_et_like_beta"] - out["highest_ez_like_beta"]  # type: ignore[operator]
    else:
        out["ez_to_et_gap"] = float("nan")
    return out


def _load_summary_rows(
    build_bin: Path,
    out_root: Path,
    backend: str,
    variant: str,
    nx: int,
    ny: int,
) -> tuple[list[dict[str, str]], list[str]]:
    cmd, run_root = _run_variant(build_bin, out_root, backend, variant, nx, ny)
    table10_path = _table10_path(run_root)
    table10_lookup = _load_table10_lookup(table10_path)

    rows: list[dict[str, str]] = []
    for point in PRIORITY_POINTS:
        raw_path = _raw_spectrum_path(run_root, point)
        matrix_path = _matrix_audit_path(run_root, point)
        raw_rows = _read_csv_rows(raw_path)
        matrix_rows = _read_csv_rows(matrix_path)
        selected_matrix = [row for row in matrix_rows if row["selected_after_matching"] == "1"]
        if len(selected_matrix) != 1:
            raise SystemExit(
                f"Expected exactly one selected matrix-audit row for {point.mode} in {_rel_from_root(matrix_path)}, "
                f"found {len(selected_matrix)}."
            )
        table_row = table10_lookup[_point_key(point.d_over_a, point.br_over_lambda0)]
        summary = _best_useful_candidates(raw_rows, table_row, raw_path)
        sel = selected_matrix[0]
        rows.append(
            {
                "backend": backend,
                "variant": variant,
                "section": SECTION,
                "case": CASE,
                "priority_point": point.key,
                "mode": point.mode,
                "mesh_label": f"{nx}x{ny}",
                "nx": str(nx),
                "ny": str(ny),
                "fem": table_row["beta_over_k0_fem_matched"],
                "ref_primary": table_row["beta_over_k0_analytic"],
                "ref_secondary": table_row["beta_over_k0_helmvec3"],
                "err_primary_pct": table_row["error_percent_analytic"],
                "err_secondary_pct": table_row["error_percent_helmvec3"],
                "ez_ratio": sel["ez_ratio"],
                "selected_rank": table_row["selected_candidate_rank"],
                "selected_eig_index": table_row["selected_eig_index"],
                "match_status": table_row["match_status"],
                "field_status": table_row["field_status"],
                "residual_rel": sel["residual_rel"],
                "P_tt_scale_rel": sel["P_tt_scale_rel"],
                "P_zz_scale_rel": sel["P_zz_scale_rel"],
                "Q_tt_scale_rel": sel["Q_tt_scale_rel"],
                "Q_tz_scale_rel": sel["Q_tz_scale_rel"],
                "Q_zt_scale_rel": sel["Q_zt_scale_rel"],
                "Q_zz_scale_rel": sel["Q_zz_scale_rel"],
                "Q_tz_Q_zt_transpose_mismatch_rel": sel["Q_tz_Q_zt_transpose_mismatch_rel"],
                "physical_kept_count": str(summary["physical_kept_count"]),
                "dedup_kept_count": str(summary["dedup_kept_count"]),
                "closest_useful_beta_to_analytic": str(summary["closest_useful_beta_to_analytic"]),
                "closest_useful_err_abs_analytic": str(summary["closest_useful_err_abs_analytic"]),
                "closest_useful_rank_to_analytic": str(summary["closest_useful_rank_to_analytic"]),
                "closest_useful_beta_to_helmvec3": str(summary["closest_useful_beta_to_helmvec3"]),
                "closest_useful_err_abs_helmvec3": str(summary["closest_useful_err_abs_helmvec3"]),
                "closest_useful_rank_to_helmvec3": str(summary["closest_useful_rank_to_helmvec3"]),
                "highest_ez_like_beta": str(summary["highest_ez_like_beta"]),
                "lowest_et_like_beta": str(summary["lowest_et_like_beta"]),
                "ez_to_et_gap": str(summary["ez_to_et_gap"]),
                "raw_spectrum_csv": _rel_from_root(raw_path),
                "matrix_audit_csv": _rel_from_root(matrix_path),
                "table10_csv": _rel_from_root(table10_path),
                "run_log": _rel_from_root(_run_log_path(run_root)),
                "run_timing_csv": _rel_from_root(_run_timing_path(run_root)),
            }
        )

    rows.sort(key=lambda row: (row["variant"], row["backend"], row["priority_point"]))
    return rows, cmd


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "variant",
        "section",
        "case",
        "priority_point",
        "mode",
        "mesh_label",
        "nx",
        "ny",
        "fem",
        "ref_primary",
        "ref_secondary",
        "err_primary_pct",
        "err_secondary_pct",
        "ez_ratio",
        "selected_rank",
        "selected_eig_index",
        "match_status",
        "field_status",
        "residual_rel",
        "P_tt_scale_rel",
        "P_zz_scale_rel",
        "Q_tt_scale_rel",
        "Q_tz_scale_rel",
        "Q_zt_scale_rel",
        "Q_zz_scale_rel",
        "Q_tz_Q_zt_transpose_mismatch_rel",
        "physical_kept_count",
        "dedup_kept_count",
        "closest_useful_beta_to_analytic",
        "closest_useful_err_abs_analytic",
        "closest_useful_rank_to_analytic",
        "closest_useful_beta_to_helmvec3",
        "closest_useful_err_abs_helmvec3",
        "closest_useful_rank_to_helmvec3",
        "highest_ez_like_beta",
        "lowest_et_like_beta",
        "ez_to_et_gap",
        "raw_spectrum_csv",
        "matrix_audit_csv",
        "table10_csv",
        "run_log",
        "run_timing_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _baseline_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if row["variant"] != "baseline":
            continue
        out[(row["backend"], row["priority_point"])] = row
    return out


def _build_markdown(rows: list[dict[str, str]], commands: dict[tuple[str, str], list[str]]) -> str:
    baseline = _baseline_lookup(rows)
    lines = [
        "# DIAG 224 Beta A/B",
        "",
        f"Generated at {dt.datetime.now().isoformat(timespec='seconds')}.",
        "",
        "## Commands",
        "",
    ]
    for variant in DEFAULT_VARIANTS:
        for backend in DEFAULT_BACKENDS:
            cmd = commands.get((variant, backend))
            if cmd:
                lines.append(
                    f"- `{variant}` / `{backend}`: `{' '.join(shlex.quote(part) for part in cmd)}`"
                )

    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Variant | Backend | Point | FEM | Analytic | HELMVEC3 | Err A% | Err H% | ez_ratio | residual_rel | Qtt rel | Qzz rel | gap Ez->Et | delta err A vs baseline |",
            "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )

    for row in rows:
        base = baseline.get((row["backend"], row["priority_point"]))
        delta = ""
        if base is not None:
            delta = f"{float(row['err_primary_pct']) - float(base['err_primary_pct']):.3f}"
        lines.append(
            "| "
            + " | ".join(
                [
                    row["variant"],
                    row["backend"],
                    row["priority_point"],
                    f"{float(row['fem']):.6f}",
                    f"{float(row['ref_primary']):.6f}",
                    f"{float(row['ref_secondary']):.6f}",
                    f"{float(row['err_primary_pct']):.3f}",
                    f"{float(row['err_secondary_pct']):.3f}",
                    f"{float(row['ez_ratio']):.6f}",
                    f"{float(row['residual_rel']):.3e}",
                    f"{float(row['Q_tt_scale_rel']):.6f}",
                    f"{float(row['Q_zz_scale_rel']):.6f}",
                    f"{float(row['ez_to_et_gap']):.6f}",
                    delta,
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
        default=Path("out/validation/diag_224_beta_ab.csv"),
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/DIAG_224_BETA_AB.md"),
    )
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=5)
    parser.add_argument("--build-bin", type=Path, default=DEFAULT_BUILD_BIN)
    args = parser.parse_args()

    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)
    build_bin = _resolve(args.build_bin)
    if not build_bin.exists():
        raise SystemExit(f"Missing executable: {_rel_from_root(build_bin)}")

    all_rows: list[dict[str, str]] = []
    commands: dict[tuple[str, str], list[str]] = {}
    for variant in DEFAULT_VARIANTS:
        for backend in DEFAULT_BACKENDS:
            rows, cmd = _load_summary_rows(build_bin, out_root, backend, variant, args.nx, args.ny)
            all_rows.extend(rows)
            commands[(variant, backend)] = cmd

    _write_csv(out_csv, all_rows)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(_build_markdown(all_rows, commands), encoding="utf-8")
    print(f"[ok] wrote {_rel_from_root(out_csv)}")
    print(f"[ok] wrote {_rel_from_root(out_md)}")


if __name__ == "__main__":
    main()
