#!/usr/bin/env python3
"""
Aggregate CSV-based scientific regression for Section 2.2.4 Figure 13 / Table 10.

Default output:
- out/validation/validation_2d_224_table10.csv
"""

from __future__ import annotations

import argparse
import csv
import shlex
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SECTION = "2.2.4"
CASE = "Figure13_Table10"


@dataclass(frozen=True)
class RequiredPoint:
    d_over_a: float
    br_over_lambda0: float


@dataclass(frozen=True)
class CaseSpec:
    case: str
    source_rel: Path
    rerun_cmd: tuple[str, ...]
    required_points: tuple[RequiredPoint, ...]
    expected_block_counts: tuple[tuple[float, int], ...]


REQUIRED_POINTS = (
    RequiredPoint(0.0, 0.5),
    RequiredPoint(0.0, 0.6),
    RequiredPoint(0.0, 0.7),
    RequiredPoint(0.0, 0.8),
    RequiredPoint(0.0, 0.9),
    RequiredPoint(0.0, 1.0),
    RequiredPoint(0.167, 0.5),
    RequiredPoint(0.167, 0.6),
    RequiredPoint(0.167, 0.7),
    RequiredPoint(0.167, 0.8),
    RequiredPoint(0.167, 0.9),
    RequiredPoint(0.167, 1.0),
    RequiredPoint(0.286, 0.5),
    RequiredPoint(0.286, 0.6),
    RequiredPoint(0.286, 0.7),
    RequiredPoint(0.286, 0.8),
    RequiredPoint(0.286, 0.9),
    RequiredPoint(0.286, 1.0),
    RequiredPoint(0.375, 0.5),
    RequiredPoint(0.375, 0.6),
    RequiredPoint(0.375, 0.7),
    RequiredPoint(0.375, 0.8),
    RequiredPoint(0.375, 0.9),
    RequiredPoint(0.375, 1.0),
    RequiredPoint(0.5, 0.4),
    RequiredPoint(0.5, 0.5),
    RequiredPoint(0.5, 0.6),
    RequiredPoint(0.5, 0.7),
    RequiredPoint(0.5, 0.8),
    RequiredPoint(0.5, 0.9),
    RequiredPoint(0.5, 1.0),
    RequiredPoint(0.6, 0.4),
    RequiredPoint(0.6, 0.5),
    RequiredPoint(0.6, 0.6),
    RequiredPoint(0.6, 0.7),
    RequiredPoint(0.6, 0.8),
    RequiredPoint(0.6, 0.9),
    RequiredPoint(0.6, 1.0),
    RequiredPoint(0.8, 0.4),
    RequiredPoint(0.8, 0.5),
    RequiredPoint(0.8, 0.6),
    RequiredPoint(0.8, 0.7),
    RequiredPoint(0.8, 0.8),
    RequiredPoint(0.8, 0.9),
    RequiredPoint(0.8, 1.0),
)


CASE_SPEC = CaseSpec(
    case=CASE,
    source_rel=Path("helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv"),
    rerun_cmd=(
        "./build/helmvec3_fig13_rect",
        "--d-over-a-preview",
        "0.20",
        "--nx",
        "10",
        "--ny",
        "5",
    ),
    required_points=REQUIRED_POINTS,
    expected_block_counts=(
        (0.0, 6),
        (0.167, 6),
        (0.286, 6),
        (0.375, 6),
        (0.5, 7),
        (0.6, 7),
        (0.8, 7),
    ),
)

TOTAL_EXPECTED_ROWS = len(REQUIRED_POINTS)


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


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


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _point_key_from_row(row: dict[str, str], *, path: Path) -> tuple[float, float]:
    return _point_key(
        _to_float(row["d_over_a"], field="d_over_a", path=path),
        _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path),
    )


def _mode_label(point: RequiredPoint) -> str:
    return f"d/a={point.d_over_a},br/lambda0={point.br_over_lambda0}"


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required_columns = {
        "d_over_a",
        "br_over_lambda0",
        "beta_over_k0_fem_matched",
        "beta_over_k0_analytic",
        "beta_over_k0_helmvec3",
        "selected_candidate_rank",
        "selected_eig_index",
        "ez_ratio",
        "error_percent_analytic",
        "error_percent_helmvec3",
        "match_status",
        "field_status",
        "et_fields_csv_file",
        "ez_fields_csv_file",
        "et_vtk_file",
        "ez_vtk_file",
    }
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    missing_columns = sorted(required_columns - set(rows[0].keys()))
    if missing_columns:
        raise SystemExit(
            f"Missing required columns in {_rel_from_root(path)}: {', '.join(missing_columns)}"
        )
    return rows


def _group_rows(rows: list[dict[str, str]], *, path: Path) -> dict[tuple[float, float], list[dict[str, str]]]:
    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(_point_key_from_row(row, path=path), []).append(row)
    return grouped


def _select_row(candidates: list[dict[str, str]], *, path: Path) -> dict[str, str]:
    return min(
        candidates,
        key=lambda row: (
            _to_float(row["error_percent_analytic"], field="error_percent_analytic", path=path),
            _to_float(row["error_percent_helmvec3"], field="error_percent_helmvec3", path=path),
            _to_int(row["selected_candidate_rank"], field="selected_candidate_rank", path=path),
            _to_int(row["selected_eig_index"], field="selected_eig_index", path=path),
        ),
    )


def _rerun_hint(out_root: Path, backend: str) -> str:
    parts: list[str] = []
    default_out_root = (ROOT / "out").resolve()
    out_root_resolved = out_root.resolve()
    if out_root_resolved != default_out_root:
        parts.append(f"TP3485_OUT_DIR={shlex.quote(str(out_root_resolved))}")
    parts.extend(shlex.quote(part) for part in CASE_SPEC.rerun_cmd)
    parts.extend(["--backend", shlex.quote(backend)])
    return " ".join(parts)


def _fail_with_missing_points(
    *,
    source_csv: Path,
    missing_points: list[RequiredPoint],
    out_root: Path,
    backend: str,
) -> None:
    lines = [f"Missing required Table 10 point(s) in {_rel_from_root(source_csv)}:"]
    for point in missing_points:
        lines.append(f"- {_mode_label(point)}")
    lines.append("Re-run helmvec3_fig13_rect to regenerate the canonical CSV, for example:")
    lines.append(f"  {_rerun_hint(out_root, backend)}")
    raise SystemExit("\n".join(lines))


def _build_output_row(
    *,
    backend: str,
    point: RequiredPoint,
    selected: dict[str, str],
    duplicate_count: int,
    source_csv: Path,
) -> dict[str, str]:
    return {
        "backend": backend,
        "section": SECTION,
        "case": CASE_SPEC.case,
        "mode": _mode_label(point),
        "fem": selected["beta_over_k0_fem_matched"],
        "ref_primary": selected["beta_over_k0_analytic"],
        "ref_secondary": selected["beta_over_k0_helmvec3"],
        "err_primary_pct": selected["error_percent_analytic"],
        "err_secondary_pct": selected["error_percent_helmvec3"],
        "d_over_a": selected["d_over_a"],
        "br_over_lambda0": selected["br_over_lambda0"],
        "selected_rank": selected["selected_candidate_rank"],
        "selected_eig_index": selected["selected_eig_index"],
        "ez_ratio": selected["ez_ratio"],
        "match_status": selected["match_status"],
        "field_status": selected["field_status"],
        "duplicate_count": str(duplicate_count),
        "et_fields_csv_file": selected["et_fields_csv_file"],
        "ez_fields_csv_file": selected["ez_fields_csv_file"],
        "et_vtk_file": selected["et_vtk_file"],
        "ez_vtk_file": selected["ez_vtk_file"],
        "source_csv": _rel_from_root(source_csv),
    }


def _validate_expected_counts(rows: list[dict[str, str]]) -> None:
    mismatches: list[str] = []
    block_counts = Counter(round(float(row["d_over_a"]), 12) for row in rows)

    for d_over_a, expected_count in CASE_SPEC.expected_block_counts:
        actual = block_counts.get(round(d_over_a, 12), 0)
        if actual != expected_count:
            mismatches.append(
                f"- d_over_a={d_over_a}: expected {expected_count} row(s), got {actual}"
            )

    total_rows = len(rows)
    if total_rows != TOTAL_EXPECTED_ROWS:
        mismatches.append(f"- total: expected {TOTAL_EXPECTED_ROWS} row(s), got {total_rows}")

    unexpected_blocks = sorted(set(block_counts) - {round(d_over_a, 12) for d_over_a, _ in CASE_SPEC.expected_block_counts})
    for block in unexpected_blocks:
        mismatches.append(f"- unexpected d_over_a block in aggregated rows: {block} ({block_counts[block]} row(s))")

    if mismatches:
        lines = ["Aggregated row count mismatch for Section 2.2.4 Figure 13 / Table 10:"]
        lines.extend(mismatches)
        raise SystemExit("\n".join(lines))


def aggregate_case(out_root: Path, backend: str) -> list[dict[str, str]]:
    source_csv = out_root / CASE_SPEC.source_rel
    try:
        rows = _read_rows(source_csv)
    except FileNotFoundError:
        lines = ["Missing required CSV input for Section 2.2.4 Figure 13 / Table 10 aggregation:"]
        lines.append(f"- {_rel_from_root(source_csv)}")
        lines.append(f"  Re-run: {_rerun_hint(out_root, backend)}")
        raise SystemExit("\n".join(lines))

    grouped = _group_rows(rows, path=source_csv)
    expected_keys = {
        _point_key(point.d_over_a, point.br_over_lambda0): point for point in CASE_SPEC.required_points
    }
    unexpected_keys = sorted(set(grouped) - set(expected_keys), key=lambda item: (item[0], item[1]))
    if unexpected_keys:
        lines = [f"Unexpected Table 10 point(s) found in {_rel_from_root(source_csv)}:"]
        for d_over_a, br_over_lambda0 in unexpected_keys:
            lines.append(f"- d/a={d_over_a},br/lambda0={br_over_lambda0}")
        raise SystemExit("\n".join(lines))

    output_rows: list[dict[str, str]] = []
    missing_points: list[RequiredPoint] = []
    for point in CASE_SPEC.required_points:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        candidates = grouped.get(key, [])
        if not candidates:
            missing_points.append(point)
            continue
        selected = _select_row(candidates, path=source_csv)
        output_rows.append(
            _build_output_row(
                backend=backend,
                point=point,
                selected=selected,
                duplicate_count=len(candidates),
                source_csv=source_csv,
            )
        )

    if missing_points:
        _fail_with_missing_points(
            source_csv=source_csv,
            missing_points=missing_points,
            out_root=out_root,
            backend=backend,
        )

    _validate_expected_counts(output_rows)
    return output_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate Section 2.2.4 Figure 13 / Table 10 CSV regressions into one CSV."
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path("out"),
        help="Root output directory containing helmvec3 case folders.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_2d_224_table10.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--backend",
        default="closed-form",
        help="Backend label to record in the aggregated CSV (default: closed-form).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)

    rows = aggregate_case(out_root=out_root, backend=args.backend)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "mode",
        "fem",
        "ref_primary",
        "ref_secondary",
        "err_primary_pct",
        "err_secondary_pct",
        "d_over_a",
        "br_over_lambda0",
        "selected_rank",
        "selected_eig_index",
        "ez_ratio",
        "match_status",
        "field_status",
        "duplicate_count",
        "et_fields_csv_file",
        "ez_fields_csv_file",
        "et_vtk_file",
        "ez_vtk_file",
        "source_csv",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
