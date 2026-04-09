#!/usr/bin/env python3
"""
Aggregate CSV-based scientific regression for Section 2.2.2 Table 6 (mixed_rect).

Default output:
- out/validation/validation_2d_222_table6.csv
"""

from __future__ import annotations

import argparse
import csv
import shlex
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SECTION = "2.2.2"


@dataclass(frozen=True)
class RequiredMode:
    mode: str
    formulation: str
    dominant_block: str
    family: str
    mode_label: str


@dataclass(frozen=True)
class CaseSpec:
    case: str
    source_rel: Path
    rerun_cmd: tuple[str, ...]
    required_modes: tuple[RequiredMode, ...]
    required_columns: tuple[str, ...]
    expected_count: int


TABLE6_E_TE_MODES = (
    RequiredMode("TE10", "E", "edge", "TE", "TE_m1_n0"),
    RequiredMode("TE20", "E", "edge", "TE", "TE_m2_n0"),
    RequiredMode("TE30", "E", "edge", "TE", "TE_m3_n0"),
    RequiredMode("TE01", "E", "edge", "TE", "TE_m0_n1"),
    RequiredMode("TE11", "E", "edge", "TE", "TE_m1_n1"),
    RequiredMode("TE21", "E", "edge", "TE", "TE_m2_n1"),
    RequiredMode("TE31", "E", "edge", "TE", "TE_m3_n1"),
    RequiredMode("TE02", "E", "edge", "TE", "TE_m0_n2"),
)

TABLE6_H_TE_MODES = tuple(
    RequiredMode(mode.mode, "H", "scalar", "TE", mode.mode_label) for mode in TABLE6_E_TE_MODES
)

TABLE6_E_TM_MODES = (
    RequiredMode("TM11", "E", "scalar", "TM", "TM_m1_n1"),
    RequiredMode("TM21", "E", "scalar", "TM", "TM_m2_n1"),
    RequiredMode("TM31", "E", "scalar", "TM", "TM_m3_n1"),
    # The paper prints TM02 with the same analytical cutoff as TE02.
    # We preserve the paper label in `mode` and record the matched repository
    # row explicitly in `mode_label`/`family`.
    RequiredMode("TM02", "E", "edge", "TE", "TE_m0_n2"),
)

TABLE6_H_TM_MODES = (
    RequiredMode("TM11", "H", "edge", "TM", "TM_m1_n1"),
    RequiredMode("TM21", "H", "edge", "TM", "TM_m2_n1"),
    RequiredMode("TM31", "H", "edge", "TM", "TM_m3_n1"),
    RequiredMode("TM02", "H", "scalar", "TE", "TE_m0_n2"),
)

CASE_SPECS = (
    CaseSpec(
        case="mixed_rect_table6_E_TE",
        source_rel=Path("helmvec1/rect/csv/mixed_rect_modes.csv"),
        rerun_cmd=("./build/mixed_rect", "--nx", "12", "--ny", "6"),
        required_modes=TABLE6_E_TE_MODES,
        required_columns=("n",),
        expected_count=8,
    ),
    CaseSpec(
        case="mixed_rect_table6_H_TE",
        source_rel=Path("helmvec1/rect/csv/mixed_rect_modes.csv"),
        rerun_cmd=("./build/mixed_rect", "--nx", "12", "--ny", "6"),
        required_modes=TABLE6_H_TE_MODES,
        required_columns=("n",),
        expected_count=8,
    ),
    CaseSpec(
        case="mixed_rect_table6_E_TM",
        source_rel=Path("helmvec1/rect/csv/mixed_rect_modes.csv"),
        rerun_cmd=("./build/mixed_rect", "--nx", "12", "--ny", "6"),
        required_modes=TABLE6_E_TM_MODES,
        required_columns=("n",),
        expected_count=4,
    ),
    CaseSpec(
        case="mixed_rect_table6_H_TM",
        source_rel=Path("helmvec1/rect/csv/mixed_rect_modes.csv"),
        rerun_cmd=("./build/mixed_rect", "--nx", "12", "--ny", "6"),
        required_modes=TABLE6_H_TM_MODES,
        required_columns=("n",),
        expected_count=4,
    ),
)

TOTAL_EXPECTED_ROWS = sum(case.expected_count for case in CASE_SPECS)


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
        raise SystemExit(f"Invalid numeric value for {field!r} in {_rel_from_root(path)}: {raw!r}") from exc


def _to_rank(raw: str, *, path: Path) -> int:
    try:
        return int(raw)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"Invalid positive_rank in {_rel_from_root(path)}: {raw!r}") from exc


def _read_rows(path: Path, *, required_columns: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    base_required_columns = {
        "formulation",
        "dominant_block",
        "family",
        "mode_label",
        "positive_rank",
        "m",
        "kc_fem",
        "kc_ana",
        "error_percent",
        "rho_abs",
    }
    required = base_required_columns | set(required_columns)
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    missing_columns = sorted(required - set(rows[0].keys()))
    if missing_columns:
        missing = ", ".join(missing_columns)
        raise SystemExit(f"Missing required columns in {_rel_from_root(path)}: {missing}")
    return rows


def _group_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (
            row["formulation"],
            row["dominant_block"],
            row["family"],
            row["mode_label"],
        )
        grouped.setdefault(key, []).append(row)
    return grouped


def _select_row(candidates: list[dict[str, str]], *, path: Path) -> dict[str, str]:
    return min(
        candidates,
        key=lambda row: (
            _to_float(row["error_percent"], field="error_percent", path=path),
            -_to_float(row["rho_abs"], field="rho_abs", path=path),
            _to_rank(row["positive_rank"], path=path),
        ),
    )


def _rerun_hint(case_spec: CaseSpec, out_root: Path, backend: str) -> str:
    parts: list[str] = []
    default_out_root = (ROOT / "out").resolve()
    out_root_resolved = out_root.resolve()
    if out_root_resolved != default_out_root:
        parts.append(f"TP3485_OUT_DIR={shlex.quote(str(out_root_resolved))}")
    parts.extend(shlex.quote(part) for part in case_spec.rerun_cmd)
    parts.extend(["--backend", shlex.quote(backend)])
    return " ".join(parts)


def _build_output_row(
    *,
    backend: str,
    case_spec: CaseSpec,
    mode: RequiredMode,
    selected: dict[str, str],
    duplicate_count: int,
    source_csv: Path,
) -> dict[str, str]:
    return {
        "backend": backend,
        "section": SECTION,
        "case": case_spec.case,
        "mode": mode.mode,
        "fem": selected["kc_fem"],
        "ref_primary": selected["kc_ana"],
        "ref_secondary": "",
        "err_primary_pct": selected["error_percent"],
        "err_secondary_pct": "",
        "mode_label": selected["mode_label"],
        "formulation": selected["formulation"],
        "dominant_block": selected["dominant_block"],
        "family": selected["family"],
        "m": selected.get("m", ""),
        "n": selected.get("n", ""),
        "p": selected.get("p", ""),
        "rho_abs": selected["rho_abs"],
        "selected_rank": selected["positive_rank"],
        "duplicate_count": str(duplicate_count),
        "source_csv": _rel_from_root(source_csv),
    }


def _fail_with_missing_modes(
    *,
    case_spec: CaseSpec,
    source_csv: Path,
    missing_modes: list[RequiredMode],
    out_root: Path,
    backend: str,
) -> None:
    lines = [
        f"Missing required mode_label(s) in {_rel_from_root(source_csv)} for {case_spec.case}:",
    ]
    for mode in missing_modes:
        lines.append(
            "- "
            f"mode={mode.mode} formulation={mode.formulation} dominant_block={mode.dominant_block} "
            f"family={mode.family} mode_label={mode.mode_label}"
        )
    lines.append("Re-run mixed_rect to regenerate the canonical CSV, for example:")
    lines.append(f"  {_rerun_hint(case_spec, out_root, backend)}")
    raise SystemExit("\n".join(lines))


def _validate_expected_counts(rows: list[dict[str, str]]) -> None:
    case_counts = Counter(row["case"] for row in rows)
    mismatches: list[str] = []

    for case_spec in CASE_SPECS:
        actual = case_counts.get(case_spec.case, 0)
        if actual != case_spec.expected_count:
            mismatches.append(
                f"- {case_spec.case}: expected {case_spec.expected_count} row(s), got {actual}"
            )

    total_rows = len(rows)
    if total_rows != TOTAL_EXPECTED_ROWS:
        mismatches.append(f"- total: expected {TOTAL_EXPECTED_ROWS} row(s), got {total_rows}")

    unexpected_cases = sorted(set(case_counts) - {case.case for case in CASE_SPECS})
    for case in unexpected_cases:
        mismatches.append(f"- unexpected case in aggregated rows: {case} ({case_counts[case]} row(s))")

    if mismatches:
        lines = ["Aggregated row count mismatch for Section 2.2.2 Table 6:"]
        lines.extend(mismatches)
        raise SystemExit("\n".join(lines))


def aggregate_cases(out_root: Path, backend: str) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    missing_files: list[tuple[CaseSpec, Path]] = []
    cached_rows: dict[Path, list[dict[str, str]]] = {}
    cached_groups: dict[Path, dict[tuple[str, str, str, str], list[dict[str, str]]]] = {}

    for case_spec in CASE_SPECS:
        source_csv = out_root / case_spec.source_rel
        if source_csv not in cached_rows:
            try:
                rows = _read_rows(source_csv, required_columns=case_spec.required_columns)
            except FileNotFoundError:
                missing_files.append((case_spec, source_csv))
                continue
            cached_rows[source_csv] = rows
            cached_groups[source_csv] = _group_rows(rows)

        grouped = cached_groups[source_csv]
        missing_modes: list[RequiredMode] = []
        for mode in case_spec.required_modes:
            key = (mode.formulation, mode.dominant_block, mode.family, mode.mode_label)
            candidates = grouped.get(key, [])
            if not candidates:
                missing_modes.append(mode)
                continue
            selected = _select_row(candidates, path=source_csv)
            output_rows.append(
                _build_output_row(
                    backend=backend,
                    case_spec=case_spec,
                    mode=mode,
                    selected=selected,
                    duplicate_count=len(candidates),
                    source_csv=source_csv,
                )
            )

        if missing_modes:
            _fail_with_missing_modes(
                case_spec=case_spec,
                source_csv=source_csv,
                missing_modes=missing_modes,
                out_root=out_root,
                backend=backend,
            )

    if missing_files:
        lines = ["Missing required CSV input(s) for Section 2.2.2 Table 6 aggregation:"]
        for case_spec, source_csv in missing_files:
            lines.append(f"- {_rel_from_root(source_csv)}")
            lines.append(f"  Re-run: {_rerun_hint(case_spec, out_root, backend)}")
        raise SystemExit("\n".join(lines))

    _validate_expected_counts(output_rows)
    return output_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate Section 2.2.2 Table 6 CSV regressions into one CSV.")
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path("out"),
        help="Root output directory containing helmvec1 case folders.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_2d_222_table6.csv"),
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

    rows = aggregate_cases(out_root=out_root, backend=args.backend)
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
        "mode_label",
        "formulation",
        "dominant_block",
        "family",
        "m",
        "n",
        "p",
        "rho_abs",
        "selected_rank",
        "duplicate_count",
        "source_csv",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {_rel_from_root(out_csv)}")


if __name__ == "__main__":
    main()
