#!/usr/bin/env python3
"""
Consolidate the raw HELMVEC3 spectrum diagnostic for Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_raw_spectrum.csv
- out/validation/DIAG_224_RAW_SPECTRUM.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import shlex
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_ROOT = ROOT / "out"
SECTION = "2.2.4"
CASE = "Figure13_Table10"


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


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _label_number(value: float) -> str:
    text = f"{value:.12g}"
    text = text.replace(".", "_")
    return text


def _raw_csv_path(out_root: Path, point: PriorityPoint) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(point.d_over_a)
        + "_br"
        + _label_number(point.br_over_lambda0)
        + "_raw_spectrum.csv"
    )
    return out_root / "helmvec3" / "fig13_rect" / "csv" / filename


def _table10_csv_path(out_root: Path) -> Path:
    return out_root / "helmvec3" / "fig13_rect" / "csv" / "helmvec3_fig13_rect_table10.csv"


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


def _load_table10_rows(path: Path) -> dict[tuple[float, float], dict[str, str]]:
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


def _missing_raw_message(out_root: Path, point: PriorityPoint) -> str:
    raw_path = _raw_csv_path(out_root, point)
    cmd = (
        f"TP3485_OUT_DIR={shlex.quote(str(out_root))} "
        "TP3485_HELMVEC3_DIAG_RAW_SPECTRUM=1 "
        "./build/helmvec3_fig13_rect --d-over-a-preview 0.20 --nx 10 --ny 5 --backend closed-form"
    )
    return (
        f"Missing raw spectrum CSV for {point.mode}: {_rel_from_root(raw_path)}. "
        f"Re-run the diagnostic export with: {cmd}"
    )


def _build_rows(out_root: Path) -> list[dict[str, str]]:
    table10_path = _table10_csv_path(out_root)
    table10_rows = _load_table10_rows(table10_path)

    rows: list[dict[str, str]] = []
    for point in PRIORITY_POINTS:
        raw_path = _raw_csv_path(out_root, point)
        if not raw_path.exists():
            raise SystemExit(_missing_raw_message(out_root, point))
        raw_rows = _read_csv_rows(raw_path)
        table_row = table10_rows[_point_key(point.d_over_a, point.br_over_lambda0)]

        for row in raw_rows:
            rows.append(
                {
                    "backend": "closed-form",
                    "section": SECTION,
                    "case": CASE,
                    "priority_point": point.key,
                    "mode": point.mode,
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
                    "raw_spectrum_csv": _rel_from_root(raw_path),
                    "table10_csv": _rel_from_root(table10_path),
                }
            )

    rows.sort(
        key=lambda row: (
            row["priority_point"],
            _to_int(row["ordered_rank"], field="ordered_rank", path=ROOT / row["raw_spectrum_csv"]),
        )
    )
    return rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "priority_point",
        "mode",
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
        "raw_spectrum_csv",
        "table10_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _build_markdown(rows: list[dict[str, str]]) -> str:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["priority_point"], []).append(row)

    lines = [
        "# DIAG 224 Raw Spectrum",
        "",
        f"Generated at {dt.datetime.now().isoformat(timespec='seconds')}.",
        "",
        "| Point | Total eigs | Physical kept | Dedup kept | Selected rank | Selected eig | Nearest physical to analytic | Nearest dedup to analytic |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]

    for point in PRIORITY_POINTS:
        point_rows = grouped.get(point.key, [])
        if not point_rows:
            continue
        total = len(point_rows)
        physical = sum(row["kept_after_physical_filter"] == "1" for row in point_rows)
        dedup = sum(row["kept_after_dedup"] == "1" for row in point_rows)
        selected = next((row for row in point_rows if row["selected_after_matching"] == "1"), None)
        analytic = float(point_rows[0]["ref_primary"])

        def _closest(target_rows: list[dict[str, str]]) -> str:
            if not target_rows:
                return "n/a"
            best = min(
                target_rows,
                key=lambda row: abs(float(row["beta_ratio_if_real_positive"]) - analytic),
            )
            return (
                f"{float(best['beta_ratio_if_real_positive']):.6f}"
                f" (|delta|={abs(float(best['beta_ratio_if_real_positive']) - analytic):.6f})"
            )

        physical_rows = [row for row in point_rows if row["kept_after_physical_filter"] == "1"]
        dedup_rows = [row for row in point_rows if row["kept_after_dedup"] == "1"]
        lines.append(
            "| "
            + point.key
            + " | "
            + str(total)
            + " | "
            + str(physical)
            + " | "
            + str(dedup)
            + " | "
            + (selected["selected_rank"] if selected is not None else "0")
            + " | "
            + (selected["selected_eig_index"] if selected is not None else "-1")
            + " | "
            + _closest(physical_rows)
            + " | "
            + _closest(dedup_rows)
            + " |"
        )

    lines.extend(
        [
            "",
            "Interpretation:",
            "- If a target-like branch appears among `kept_after_physical_filter=1` rows but disappears from `kept_after_dedup=1`, the loss is happening in deduplication.",
            "- If it never appears among `kept_after_physical_filter=1` rows but positive `lambda_real` values exist nearby, the physical filter is the first suspect.",
            "- If even the raw positive-real spectrum never approaches the target branch, the branch is already missing in the solved discretized problem.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Consolidate the raw HELMVEC3 spectrum diagnostic for Figure 13 / Table 10."
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_OUT_ROOT,
        help="TP3485 output root containing out/helmvec3/fig13_rect/csv.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/diag_224_raw_spectrum.csv"),
        help="Consolidated diagnostic CSV output.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/DIAG_224_RAW_SPECTRUM.md"),
        help="Optional Markdown summary output.",
    )
    args = parser.parse_args()

    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)

    rows = _build_rows(out_root)
    _write_csv(out_csv, rows)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(_build_markdown(rows), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
