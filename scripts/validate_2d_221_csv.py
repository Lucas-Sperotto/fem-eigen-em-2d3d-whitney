#!/usr/bin/env python3
"""
Aggregate CSV-based scientific regression for Section 2.2.1 (Tables 4 and 5).

Default output:
- out/validation/validation_2d_221.csv
"""

from __future__ import annotations

import argparse
import csv
import shlex
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SECTION = "2.2.1"


@dataclass(frozen=True)
class RequiredMode:
    mode: str
    family: str
    mode_label: str


@dataclass(frozen=True)
class CaseSpec:
    case: str
    source_rel: Path
    rerun_cmd: tuple[str, ...]
    required_modes: tuple[RequiredMode, ...]


CASE_SPECS = (
    CaseSpec(
        case="edge_rect_table4",
        source_rel=Path("helmvec/rect/csv/edge_rect_modes.csv"),
        rerun_cmd=("./build/edge_rect", "--nx", "10", "--ny", "20", "--nmodos", "10"),
        required_modes=(
            RequiredMode("TE10", "TE", "TE_m1_n0"),
            RequiredMode("TE20", "TE", "TE_m2_n0"),
            RequiredMode("TE30", "TE", "TE_m3_n0"),
            RequiredMode("TE01", "TE", "TE_m0_n1"),
            RequiredMode("TE11", "TE", "TE_m1_n1"),
            RequiredMode("TE21", "TE", "TE_m2_n1"),
            RequiredMode("TE31", "TE", "TE_m3_n1"),
            RequiredMode("TE02", "TE", "TE_m0_n2"),
            # The paper/translation labels are preserved in `mode`; the matched
            # repository row is recorded separately in `mode_label`.
            RequiredMode("TM11", "TM", "TM_m1_n1"),
            RequiredMode("TM12", "TM", "TM_m2_n1"),
            RequiredMode("TM13", "TM", "TM_m3_n1"),
            RequiredMode("TM21", "TM", "TM_m1_n2"),
            RequiredMode("TM22", "TM", "TM_m2_n2"),
            RequiredMode("TM23", "TM", "TM_m3_n2"),
            RequiredMode("TM31", "TM", "TM_m1_n3"),
            RequiredMode("TM32", "TM", "TM_m2_n3"),
        ),
    ),
    CaseSpec(
        case="edge_circle_table5",
        source_rel=Path("helmvec/circle/csv/edge_circle_modes.csv"),
        rerun_cmd=("./build/edge_circle", "--nr", "8", "--nt", "15", "--nmodos", "10"),
        required_modes=(
            RequiredMode("TE01", "TE", "TE_m0_p1"),
            RequiredMode("TE11", "TE", "TE_m1_p1"),
            RequiredMode("TE12", "TE", "TE_m2_p1"),
            RequiredMode("TE13", "TE", "TE_m3_p1"),
            RequiredMode("TE20", "TE", "TE_m0_p2"),
            RequiredMode("TE21", "TE", "TE_m1_p2"),
            RequiredMode("TE22", "TE", "TE_m2_p2"),
            RequiredMode("TE23", "TE", "TE_m3_p2"),
            RequiredMode("TM10", "TM", "TM_m0_p1"),
            RequiredMode("TM11", "TM", "TM_m1_p1"),
            RequiredMode("TM12", "TM", "TM_m2_p1"),
            RequiredMode("TM13", "TM", "TM_m3_p1"),
            RequiredMode("TM20", "TM", "TM_m0_p2"),
            RequiredMode("TM21", "TM", "TM_m1_p2"),
            RequiredMode("TM22", "TM", "TM_m2_p2"),
            RequiredMode("TM23", "TM", "TM_m3_p2"),
        ),
    ),
)


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


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required_columns = {
        "family",
        "mode_label",
        "positive_rank",
        "m",
        "kc_fem",
        "kc_ana",
        "error_percent",
        "rho_abs",
    }
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    missing_columns = sorted(required_columns - set(rows[0].keys()))
    if missing_columns:
        missing = ", ".join(missing_columns)
        raise SystemExit(f"Missing required columns in {_rel_from_root(path)}: {missing}")
    return rows


def _group_rows(rows: list[dict[str, str]]) -> dict[tuple[str, str], list[dict[str, str]]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (row["family"], row["mode_label"])
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
        lines.append(f"- mode={mode.mode} family={mode.family} mode_label={mode.mode_label}")
    lines.append("Re-run the corresponding executable with sufficient mode export, for example:")
    lines.append(f"  {_rerun_hint(case_spec, out_root, backend)}")
    raise SystemExit("\n".join(lines))


def aggregate_cases(out_root: Path, backend: str) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []
    missing_files: list[tuple[CaseSpec, Path]] = []

    for case_spec in CASE_SPECS:
        source_csv = out_root / case_spec.source_rel
        try:
            rows = _read_rows(source_csv)
        except FileNotFoundError:
            missing_files.append((case_spec, source_csv))
            continue

        grouped = _group_rows(rows)
        missing_modes: list[RequiredMode] = []
        for mode in case_spec.required_modes:
            candidates = grouped.get((mode.family, mode.mode_label), [])
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
        lines = ["Missing required CSV input(s) for Section 2.2.1 aggregation:"]
        for case_spec, source_csv in missing_files:
            lines.append(f"- {_rel_from_root(source_csv)}")
            lines.append(f"  Re-run: {_rerun_hint(case_spec, out_root, backend)}")
        raise SystemExit("\n".join(lines))

    return output_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate Section 2.2.1 CSV regressions into one CSV.")
    parser.add_argument("--out-root", type=Path, default=Path("out"), help="Root output directory containing helmvec case folders.")
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_2d_221.csv"),
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
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
