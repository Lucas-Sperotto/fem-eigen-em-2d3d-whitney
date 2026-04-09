#!/usr/bin/env python3
"""
Run and consolidate the HELMVEC3 matrix audit diagnostic for Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_matrix_audit.csv
- out/validation/DIAG_224_MATRIX_AUDIT.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import shlex
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_ROOT = ROOT / "out"
DEFAULT_BUILD_BIN = ROOT / "build" / "helmvec3_fig13_rect"
DEFAULT_RUNS_SUBDIR = ".diag_224_matrix_audit"
SECTION = "2.2.4"
CASE = "Figure13_Table10"
DEFAULT_BACKENDS = ("closed-form", "gauss")


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


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _backend_tag(backend: str) -> str:
    return backend.replace("-", "_")


def _run_root(out_root: Path, backend: str) -> Path:
    return out_root / DEFAULT_RUNS_SUBDIR / _backend_tag(backend)


def _matrix_audit_path(out_root: Path, point: PriorityPoint) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(point.d_over_a)
        + "_br"
        + _label_number(point.br_over_lambda0)
        + "_matrix_audit.csv"
    )
    return out_root / "helmvec3" / "fig13_rect" / "csv" / filename


def _table10_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "csv" / "helmvec3_fig13_rect_table10.csv"


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


def _table10_lookup(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    rows = _read_csv_rows(path)
    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in rows:
        key = _point_key(
            _to_float(row["d_over_a"], field="d_over_a", path=path),
            _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path),
        )
        grouped.setdefault(key, []).append(row)

    selected: dict[tuple[float, float], dict[str, str]] = {}
    for point in PRIORITY_POINTS:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        matches = grouped.get(key, [])
        if len(matches) != 1:
            raise SystemExit(
                f"Expected exactly one Table 10 row for {point.mode} in {_rel_from_root(path)}, "
                f"found {len(matches)}."
            )
        selected[key] = matches[0]
    return selected


def _run_backend(build_bin: Path, out_root: Path, backend: str, nx: int, ny: int) -> tuple[list[str], Path]:
    run_root = _run_root(out_root, backend)
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(run_root)
    env["TP3485_HELMVEC3_DIAG_MATRIX_AUDIT"] = "1"
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


def _load_backend_rows(
    build_bin: Path,
    out_root: Path,
    backend: str,
    nx: int,
    ny: int,
) -> tuple[list[dict[str, str]], list[str]]:
    build_cmd, run_root = _run_backend(build_bin, out_root, backend, nx, ny)
    table10_path = _table10_path(run_root)
    table10_rows = _table10_lookup(table10_path)

    rows: list[dict[str, str]] = []
    for point in PRIORITY_POINTS:
        matrix_path = _matrix_audit_path(run_root, point)
        if not matrix_path.exists():
            raise SystemExit(
                f"Missing matrix audit CSV for {point.mode}: {_rel_from_root(matrix_path)}. "
                f"Re-run with TP3485_HELMVEC3_DIAG_MATRIX_AUDIT=1."
            )
        matrix_rows = _read_csv_rows(matrix_path)
        table_row = table10_rows[_point_key(point.d_over_a, point.br_over_lambda0)]
        for row in matrix_rows:
            rows.append(
                {
                    "backend": backend,
                    "section": SECTION,
                    "case": CASE,
                    "priority_point": point.key,
                    "mode": point.mode,
                    "mesh_label": f"{nx}x{ny}",
                    "nx": str(nx),
                    "ny": str(ny),
                    "d_over_a": row["d_over_a"],
                    "br_over_lambda0": row["br_over_lambda0"],
                    "fem": table_row["beta_over_k0_fem_matched"],
                    "ref_primary": table_row["beta_over_k0_analytic"],
                    "ref_secondary": table_row["beta_over_k0_helmvec3"],
                    "err_primary_pct": table_row["error_percent_analytic"],
                    "err_secondary_pct": table_row["error_percent_helmvec3"],
                    "selected_rank": table_row["selected_candidate_rank"],
                    "selected_eig_index": table_row["selected_eig_index"],
                    "match_status": table_row["match_status"],
                    "field_status": table_row["field_status"],
                    "ordered_rank": row["ordered_rank"],
                    "solver_index": row["solver_index"],
                    "lambda_real": row["lambda_real"],
                    "lambda_imag": row["lambda_imag"],
                    "beta_ratio_if_real_positive": row["beta_ratio_if_real_positive"],
                    "filter_reason": row["filter_reason"],
                    "kept_after_physical_filter": row["kept_after_physical_filter"],
                    "kept_after_dedup": row["kept_after_dedup"],
                    "candidate_rank_after_dedup": row["candidate_rank_after_dedup"],
                    "selected_after_matching": row["selected_after_matching"],
                    "et_energy": row["et_energy"],
                    "ez_energy": row["ez_energy"],
                    "ez_ratio": row["ez_ratio"],
                    "residual_abs": row["residual_abs"],
                    "residual_rel": row["residual_rel"],
                    "P_fro": row["P_fro"],
                    "Q_fro": row["Q_fro"],
                    "P_tt_fro": row["P_tt_fro"],
                    "P_zz_fro": row["P_zz_fro"],
                    "Q_tt_fro": row["Q_tt_fro"],
                    "Q_tz_fro": row["Q_tz_fro"],
                    "Q_zt_fro": row["Q_zt_fro"],
                    "Q_zz_fro": row["Q_zz_fro"],
                    "P_tt_asym_rel": row["P_tt_asym_rel"],
                    "P_zz_asym_rel": row["P_zz_asym_rel"],
                    "Q_tt_asym_rel": row["Q_tt_asym_rel"],
                    "Q_zz_asym_rel": row["Q_zz_asym_rel"],
                    "Q_tz_Q_zt_transpose_mismatch_rel": row["Q_tz_Q_zt_transpose_mismatch_rel"],
                    "block_norm_max": row["block_norm_max"],
                    "P_tt_scale_rel": row["P_tt_scale_rel"],
                    "P_zz_scale_rel": row["P_zz_scale_rel"],
                    "Q_tt_scale_rel": row["Q_tt_scale_rel"],
                    "Q_tz_scale_rel": row["Q_tz_scale_rel"],
                    "Q_zt_scale_rel": row["Q_zt_scale_rel"],
                    "Q_zz_scale_rel": row["Q_zz_scale_rel"],
                    "matrix_audit_csv": _rel_from_root(matrix_path),
                    "table10_csv": _rel_from_root(table10_path),
                    "run_log": _rel_from_root(_run_log_path(run_root)),
                    "run_timing_csv": _rel_from_root(_run_timing_path(run_root)),
                }
            )

    rows.sort(
        key=lambda row: (
            row["backend"],
            row["priority_point"],
            _to_int(row["ordered_rank"], field="ordered_rank", path=ROOT / row["matrix_audit_csv"]),
        )
    )
    return rows, build_cmd


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "priority_point",
        "mode",
        "mesh_label",
        "nx",
        "ny",
        "d_over_a",
        "br_over_lambda0",
        "fem",
        "ref_primary",
        "ref_secondary",
        "err_primary_pct",
        "err_secondary_pct",
        "selected_rank",
        "selected_eig_index",
        "match_status",
        "field_status",
        "ordered_rank",
        "solver_index",
        "lambda_real",
        "lambda_imag",
        "beta_ratio_if_real_positive",
        "filter_reason",
        "kept_after_physical_filter",
        "kept_after_dedup",
        "candidate_rank_after_dedup",
        "selected_after_matching",
        "et_energy",
        "ez_energy",
        "ez_ratio",
        "residual_abs",
        "residual_rel",
        "P_fro",
        "Q_fro",
        "P_tt_fro",
        "P_zz_fro",
        "Q_tt_fro",
        "Q_tz_fro",
        "Q_zt_fro",
        "Q_zz_fro",
        "P_tt_asym_rel",
        "P_zz_asym_rel",
        "Q_tt_asym_rel",
        "Q_zz_asym_rel",
        "Q_tz_Q_zt_transpose_mismatch_rel",
        "block_norm_max",
        "P_tt_scale_rel",
        "P_zz_scale_rel",
        "Q_tt_scale_rel",
        "Q_tz_scale_rel",
        "Q_zt_scale_rel",
        "Q_zz_scale_rel",
        "matrix_audit_csv",
        "table10_csv",
        "run_log",
        "run_timing_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _selected_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    selected: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if row["selected_after_matching"] != "1":
            continue
        key = (row["backend"], row["priority_point"])
        selected[key] = row
    return selected


def _build_markdown(rows: list[dict[str, str]], commands: dict[str, list[str]]) -> str:
    selected = _selected_rows(rows)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["backend"], row["priority_point"]), []).append(row)

    lines = [
        "# DIAG 224 Matrix Audit",
        "",
        f"Generated at {dt.datetime.now().isoformat(timespec='seconds')}.",
        "",
        "## Commands",
        "",
    ]
    for backend in DEFAULT_BACKENDS:
        cmd = commands.get(backend)
        if cmd:
            lines.append(f"- `{backend}`: `{' '.join(shlex.quote(part) for part in cmd)}`")
    lines.extend(
        [
            "",
            "## Selected Candidates",
            "",
            "| Backend | Point | FEM | Analytic | HELMVEC3 | Err A% | Err H% | Sel rank | Sel eig | ez_ratio | residual_rel | Physical kept | Dedup kept | Q_tz/Q_zt mismatch |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )

    for backend in DEFAULT_BACKENDS:
        for point in PRIORITY_POINTS:
            key = (backend, point.key)
            point_rows = grouped.get(key, [])
            if not point_rows or key not in selected:
                continue
            chosen = selected[key]
            physical_kept = sum(row["kept_after_physical_filter"] == "1" for row in point_rows)
            dedup_kept = sum(row["kept_after_dedup"] == "1" for row in point_rows)
            lines.append(
                "| "
                + " | ".join(
                    [
                        backend,
                        point.key,
                        f"{float(chosen['fem']):.6f}",
                        f"{float(chosen['ref_primary']):.6f}",
                        f"{float(chosen['ref_secondary']):.6f}",
                        f"{float(chosen['err_primary_pct']):.3f}",
                        f"{float(chosen['err_secondary_pct']):.3f}",
                        chosen["selected_rank"],
                        chosen["selected_eig_index"],
                        f"{float(chosen['ez_ratio']):.6f}",
                        f"{float(chosen['residual_rel']):.3e}",
                        str(physical_kept),
                        str(dedup_kept),
                        f"{float(chosen['Q_tz_Q_zt_transpose_mismatch_rel']):.3e}",
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Backend Comparison",
            "",
            "| Point | Selected beta closed-form | Selected beta gauss | abs delta | rel residual closed-form | rel residual gauss |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for point in PRIORITY_POINTS:
        cf = selected.get(("closed-form", point.key))
        gq = selected.get(("gauss", point.key))
        if cf is None or gq is None:
            continue
        beta_cf = float(cf["beta_ratio_if_real_positive"])
        beta_gq = float(gq["beta_ratio_if_real_positive"])
        lines.append(
            "| "
            + " | ".join(
                [
                    point.key,
                    f"{beta_cf:.6f}",
                    f"{beta_gq:.6f}",
                    f"{abs(beta_cf - beta_gq):.6f}",
                    f"{float(cf['residual_rel']):.3e}",
                    f"{float(gq['residual_rel']):.3e}",
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Filter Summary",
            "",
        ]
    )
    for backend in DEFAULT_BACKENDS:
        backend_rows = [row for row in rows if row["backend"] == backend]
        counts = Counter(row["filter_reason"] for row in backend_rows)
        lines.append(f"- `{backend}`: {dict(counts)}")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_OUT_ROOT,
        help="Base output root used for temporary runs and final validation outputs.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/diag_224_matrix_audit.csv"),
        help="Final consolidated CSV path.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/DIAG_224_MATRIX_AUDIT.md"),
        help="Final markdown summary path.",
    )
    parser.add_argument("--nx", type=int, default=10, help="Mesh divisions in x for both backends.")
    parser.add_argument("--ny", type=int, default=5, help="Mesh divisions in y for both backends.")
    parser.add_argument(
        "--build-bin",
        type=Path,
        default=DEFAULT_BUILD_BIN,
        help="Path to helmvec3_fig13_rect executable.",
    )
    args = parser.parse_args()

    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)
    build_bin = _resolve(args.build_bin)

    if not build_bin.exists():
        raise SystemExit(f"Missing executable: {_rel_from_root(build_bin)}")

    all_rows: list[dict[str, str]] = []
    commands: dict[str, list[str]] = {}
    for backend in DEFAULT_BACKENDS:
        backend_rows, cmd = _load_backend_rows(build_bin, out_root, backend, args.nx, args.ny)
        commands[backend] = cmd
        all_rows.extend(backend_rows)

    _write_csv(out_csv, all_rows)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(_build_markdown(all_rows, commands), encoding="utf-8")
    print(f"[ok] wrote {_rel_from_root(out_csv)}")
    print(f"[ok] wrote {_rel_from_root(out_md)}")


if __name__ == "__main__":
    main()
