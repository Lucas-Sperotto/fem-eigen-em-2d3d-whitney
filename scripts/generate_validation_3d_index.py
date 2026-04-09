#!/usr/bin/env python3
"""
Generate a consolidated 3D scientific regression index from existing validation CSVs.

Default outputs:
- out/validation/validation_3d_index.csv
- out/validation/VALIDATION_3D_INDEX.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class IndexSpec:
    section: str
    article_ref: str
    family: str
    executables: str
    case: str
    modes_rel: Path
    summary_rel: Path
    primary_metric: str
    primary_error_column: str
    secondary_metric: str
    secondary_error_column: str
    notes: str


SPECS = (
    IndexSpec(
        section="3.1",
        article_ref="Figura 15 / Tabela 12",
        family="FEM3D0/FEM3D1",
        executables="fem3d0_air|fem3d1_air",
        case="air",
        modes_rel=Path("validation_3d_31_modes.csv"),
        summary_rel=Path("validation_3d_31_summary.csv"),
        primary_metric="k0_fem vs k0_ana",
        primary_error_column="err_ana_pct",
        secondary_metric="k0_fem vs ref_paper",
        secondary_error_column="err_ref_pct",
        notes="Cavidade retangular preenchida com ar.",
    ),
    IndexSpec(
        section="3.1",
        article_ref="Figura 16 / Tabela 13",
        family="FEM3D0/FEM3D1",
        executables="fem3d0_half|fem3d1_half",
        case="half",
        modes_rel=Path("validation_3d_31_modes.csv"),
        summary_rel=Path("validation_3d_31_summary.csv"),
        primary_metric="k0_fem vs k0_ana",
        primary_error_column="err_ana_pct",
        secondary_metric="k0_fem vs ref_paper",
        secondary_error_column="err_ref_pct",
        notes="Cavidade retangular semi-preenchida.",
    ),
    IndexSpec(
        section="3.1",
        article_ref="Figura 17 / Tabela 14",
        family="FEM3D0/FEM3D1",
        executables="fem3d0_cyl|fem3d1_cyl",
        case="cyl",
        modes_rel=Path("validation_3d_31_modes.csv"),
        summary_rel=Path("validation_3d_31_summary.csv"),
        primary_metric="k0_fem vs k0_ana",
        primary_error_column="err_ana_pct",
        secondary_metric="k0_fem vs ref_paper",
        secondary_error_column="err_ref_pct",
        notes="Cavidade cilindrica circular com ar.",
    ),
    IndexSpec(
        section="3.1",
        article_ref="Tabela 15",
        family="FEM3D0/FEM3D1",
        executables="fem3d0_sphere|fem3d1_sphere",
        case="sphere",
        modes_rel=Path("validation_3d_31_modes.csv"),
        summary_rel=Path("validation_3d_31_summary.csv"),
        primary_metric="k0_fem vs k0_ana",
        primary_error_column="err_ana_pct",
        secondary_metric="k0_fem vs ref_paper",
        secondary_error_column="err_ref_pct",
        notes="Cavidade esferica com raio de 1 cm.",
    ),
)


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _to_float(raw: str) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _max_abs(rows: list[dict[str, str]], column: str) -> float | None:
    if not column:
        return None
    values = [abs(val) for row in rows if (val := _to_float(row.get(column, ""))) is not None]
    if not values:
        return None
    return max(values)


def _join_unique(rows: list[dict[str, str]], column: str) -> str:
    values = sorted({str(row.get(column, "")).strip() for row in rows if str(row.get(column, "")).strip()})
    return "|".join(values)


def _status_for(
    modes_exists: bool,
    summary_exists: bool,
    mode_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
) -> str:
    if not modes_exists and not summary_exists:
        return "missing_files"
    if not modes_exists:
        return "missing_modes_file"
    if not summary_exists:
        return "missing_summary_file"
    if not mode_rows and not summary_rows:
        return "case_missing"
    if not mode_rows:
        return "case_missing_modes"
    if not summary_rows:
        return "case_missing_summary"
    return "present"


def build_index(validation_dir: Path) -> list[dict[str, str]]:
    cached_rows: dict[Path, list[dict[str, str]]] = {}
    output_rows: list[dict[str, str]] = []

    for spec in SPECS:
        modes_path = validation_dir / spec.modes_rel
        summary_path = validation_dir / spec.summary_rel
        if modes_path not in cached_rows:
            cached_rows[modes_path] = _read_rows(modes_path)
        if summary_path not in cached_rows:
            cached_rows[summary_path] = _read_rows(summary_path)

        all_mode_rows = cached_rows[modes_path]
        all_summary_rows = cached_rows[summary_path]
        mode_rows = [row for row in all_mode_rows if row.get("case", "") == spec.case]
        summary_rows = [row for row in all_summary_rows if row.get("case", "") == spec.case]
        max_err_primary = _max_abs(mode_rows, spec.primary_error_column)
        max_err_secondary = _max_abs(mode_rows, spec.secondary_error_column)

        output_rows.append(
            {
                "section": spec.section,
                "article_ref": spec.article_ref,
                "family": spec.family,
                "executables": spec.executables,
                "case": spec.case,
                "modes_file": _rel_from_root(modes_path),
                "summary_file": _rel_from_root(summary_path),
                "primary_metric": spec.primary_metric,
                "primary_error_column": spec.primary_error_column,
                "secondary_metric": spec.secondary_metric,
                "secondary_error_column": spec.secondary_error_column,
                "mode_rows": str(len(mode_rows)),
                "summary_rows": str(len(summary_rows)),
                "solvers_present": _join_unique(mode_rows, "solver"),
                "backends_present": _join_unique(mode_rows, "backend"),
                "max_err_primary_pct": "" if max_err_primary is None else f"{max_err_primary:.12g}",
                "max_err_secondary_pct": "" if max_err_secondary is None else f"{max_err_secondary:.12g}",
                "status": _status_for(modes_path.exists(), summary_path.exists(), mode_rows, summary_rows),
                "notes": spec.notes,
            }
        )

    return output_rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "section",
        "article_ref",
        "family",
        "executables",
        "case",
        "modes_file",
        "summary_file",
        "primary_metric",
        "primary_error_column",
        "secondary_metric",
        "secondary_error_column",
        "mode_rows",
        "summary_rows",
        "solvers_present",
        "backends_present",
        "max_err_primary_pct",
        "max_err_secondary_pct",
        "status",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_md(path: Path, rows: list[dict[str, str]]) -> None:
    now = dt.datetime.now().isoformat(timespec="seconds")
    lines: list[str] = []
    lines.append("# Indice Consolidado de Validacao 3D")
    lines.append("")
    lines.append(f"Gerado em: `{now}`")
    lines.append("")
    lines.append("Este indice consolida os agregados cientificos 3D ja produzidos em `out/validation/`.")
    lines.append("Ele aponta o item do artigo, os executaveis envolvidos, os agregados canônicos correspondentes e as metricas principais de erro.")
    lines.append("")
    lines.append("| Secao | Tabela/Figura | Familia | Executaveis | Agregados | Mode rows | Summary rows | Max err ana (%) | Max err paper (%) | Status |")
    lines.append("|---|---|---|---|---|---:|---:|---:|---:|---|")
    for row in rows:
        modes_name = Path(row["modes_file"]).name
        summary_name = Path(row["summary_file"]).name
        aggregations = (
            f"[{modes_name}]({modes_name}) / "
            f"[{summary_name}]({summary_name})"
        )
        lines.append(
            f"| `{row['section']}` | `{row['article_ref']}` | `{row['family']}` | `{row['executables']}` | "
            f"{aggregations} | `{row['mode_rows']}` | `{row['summary_rows']}` | "
            f"`{row['max_err_primary_pct']}` | `{row['max_err_secondary_pct']}` | `{row['status']}` |"
        )
    lines.append("")
    lines.append("## Observacoes")
    lines.append("")
    for row in rows:
        if row["status"] != "present":
            lines.append(
                f"- `{row['article_ref']}` esta `{row['status']}` para o caso `{row['case']}`."
            )
    if all(row["status"] == "present" for row in rows):
        lines.append("- Todos os itens catalogados neste indice estao presentes nesta rodada.")
    lines.append("")
    lines.append("- CSV maquina-legivel: `[validation_3d_index.csv](validation_3d_index.csv)`")

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a consolidated 3D validation index from existing aggregates.")
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("out/validation"),
        help="Directory containing validation aggregate CSVs.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_3d_index.csv"),
        help="Output CSV path for the consolidated 3D index.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/VALIDATION_3D_INDEX.md"),
        help="Output Markdown path for the consolidated 3D index.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation_dir = _resolve(args.validation_dir)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)

    rows = build_index(validation_dir)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_csv, rows)
    _write_md(out_md, rows)

    print(f"Saved: {_rel_from_root(out_csv)}")
    print(f"Saved: {_rel_from_root(out_md)}")


if __name__ == "__main__":
    main()
