#!/usr/bin/env python3
"""
Generate a consolidated 2D scientific regression index from existing validation CSVs.

Default outputs:
- out/validation/validation_2d_index.csv
- out/validation/VALIDATION_2D_INDEX.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class IndexSpec:
    section: str
    article_ref: str
    family: str
    executable: str
    aggregate_rel: Path
    case_filters: tuple[str, ...]
    primary_metric: str
    primary_error_column: str
    secondary_metric: str
    secondary_error_column: str
    notes: str


SPECS = (
    IndexSpec(
        section="2.1",
        article_ref="Tabela 1",
        family="HELM10",
        executable="helm10_rect",
        aggregate_rel=Path("validation_2d_21.csv"),
        case_filters=("helm10_rect_table1",),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Guia retangular homogeneo.",
    ),
    IndexSpec(
        section="2.1",
        article_ref="Tabela 2",
        family="HELM10",
        executable="helm10_circle",
        aggregate_rel=Path("validation_2d_21.csv"),
        case_filters=("helm10_circle_table2",),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Guia circular homogeneo.",
    ),
    IndexSpec(
        section="2.1",
        article_ref="Tabela 3",
        family="HELM10",
        executable="helm10_coax",
        aggregate_rel=Path("validation_2d_21.csv"),
        case_filters=("helm10_coax_table3",),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Linha coaxial homogenea.",
    ),
    IndexSpec(
        section="2.2.1",
        article_ref="Tabela 4",
        family="HELMVEC",
        executable="edge_rect",
        aggregate_rel=Path("validation_2d_221.csv"),
        case_filters=("edge_rect_table4",),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Guia retangular vetorial.",
    ),
    IndexSpec(
        section="2.2.1",
        article_ref="Tabela 5",
        family="HELMVEC",
        executable="edge_circle",
        aggregate_rel=Path("validation_2d_221.csv"),
        case_filters=("edge_circle_table5",),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Guia circular vetorial.",
    ),
    IndexSpec(
        section="2.2.2",
        article_ref="Tabela 6",
        family="HELMVEC1",
        executable="mixed_rect",
        aggregate_rel=Path("validation_2d_222_table6.csv"),
        case_filters=(
            "mixed_rect_table6_E_TE",
            "mixed_rect_table6_H_TE",
            "mixed_rect_table6_E_TM",
            "mixed_rect_table6_H_TM",
        ),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Tabela mista retangular consolidando formulacoes E/H.",
    ),
    IndexSpec(
        section="2.2.2",
        article_ref="Tabela 7",
        family="HELMVEC1",
        executable="mixed_circle",
        aggregate_rel=Path("validation_2d_222_table7.csv"),
        case_filters=(
            "mixed_circle_table7_E_TE",
            "mixed_circle_table7_H_TE",
            "mixed_circle_table7_E_TM",
            "mixed_circle_table7_H_TM",
        ),
        primary_metric="kc_fem vs kc_ana",
        primary_error_column="err_primary_pct",
        secondary_metric="",
        secondary_error_column="",
        notes="Tabela mista circular consolidando formulacoes E/H.",
    ),
    IndexSpec(
        section="2.2.3",
        article_ref="Figura 11 / Tabela 8",
        family="HELMVEC2",
        executable="helmvec2_rect",
        aggregate_rel=Path("validation_2d_22.csv"),
        case_filters=("Figure11_Table8",),
        primary_metric="k0L_fem vs HELMVEC2 ref",
        primary_error_column="err_primary_pct",
        secondary_metric="k0L_fem vs Hayata ref",
        secondary_error_column="err_secondary_pct",
        notes="Guia parcialmente preenchido com beta fixo.",
    ),
    IndexSpec(
        section="2.2.4",
        article_ref="Figura 12 / Tabela 9",
        family="HELMVEC3",
        executable="helmvec3_fig12_rect",
        aggregate_rel=Path("validation_2d_22.csv"),
        case_filters=("Figure12_Table9",),
        primary_metric="beta_over_k0 FEM vs analitico",
        primary_error_column="err_primary_pct",
        secondary_metric="beta_over_k0 FEM vs HELMVEC3 ref",
        secondary_error_column="err_secondary_pct",
        notes="Dispersao para d/a = 0.5.",
    ),
    IndexSpec(
        section="2.2.4",
        article_ref="Figura 13 / Tabela 10",
        family="HELMVEC3",
        executable="helmvec3_fig13_rect",
        aggregate_rel=Path("validation_2d_22.csv"),
        case_filters=("Figure13_Table10",),
        primary_metric="beta_over_k0 FEM vs analitico",
        primary_error_column="err_primary_pct",
        secondary_metric="beta_over_k0 FEM vs HELMVEC3 ref",
        secondary_error_column="err_secondary_pct",
        notes="Varredura em d/a para guia parcialmente preenchido.",
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


def _max_abs(rows: Iterable[dict[str, str]], column: str) -> float | None:
    if not column:
        return None
    values = [abs(val) for row in rows if (val := _to_float(row.get(column, ""))) is not None]
    if not values:
        return None
    return max(values)


def _status_for(file_exists: bool, rows: list[dict[str, str]]) -> str:
    if not file_exists:
        return "missing_file"
    if not rows:
        return "case_missing"
    return "present"


def build_index(validation_dir: Path) -> list[dict[str, str]]:
    cached_rows: dict[Path, list[dict[str, str]]] = {}
    output_rows: list[dict[str, str]] = []

    for spec in SPECS:
        aggregate_path = validation_dir / spec.aggregate_rel
        if aggregate_path not in cached_rows:
            cached_rows[aggregate_path] = _read_rows(aggregate_path)
        all_rows = cached_rows[aggregate_path]
        case_set = set(spec.case_filters)
        matched_rows = [row for row in all_rows if row.get("case", "") in case_set]
        file_exists = aggregate_path.exists()
        max_err_primary = _max_abs(matched_rows, spec.primary_error_column)
        max_err_secondary = _max_abs(matched_rows, spec.secondary_error_column)

        output_rows.append(
            {
                "section": spec.section,
                "article_ref": spec.article_ref,
                "family": spec.family,
                "executable": spec.executable,
                "aggregate_file": _rel_from_root(aggregate_path),
                "aggregate_cases": "|".join(spec.case_filters),
                "primary_metric": spec.primary_metric,
                "primary_error_column": spec.primary_error_column,
                "secondary_metric": spec.secondary_metric,
                "secondary_error_column": spec.secondary_error_column,
                "rows": str(len(matched_rows)),
                "max_err_primary_pct": "" if max_err_primary is None else f"{max_err_primary:.12g}",
                "max_err_secondary_pct": "" if max_err_secondary is None else f"{max_err_secondary:.12g}",
                "status": _status_for(file_exists, matched_rows),
                "notes": spec.notes,
            }
        )

    return output_rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "section",
        "article_ref",
        "family",
        "executable",
        "aggregate_file",
        "aggregate_cases",
        "primary_metric",
        "primary_error_column",
        "secondary_metric",
        "secondary_error_column",
        "rows",
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
    lines.append("# Indice Consolidado de Validacao 2D")
    lines.append("")
    lines.append(f"Gerado em: `{now}`")
    lines.append("")
    lines.append("Este indice consolida os agregados cientificos 2D ja produzidos em `out/validation/`.")
    lines.append("Ele aponta o item do artigo, o executavel/familia responsavel, o agregado correspondente e a metrica principal de comparacao.")
    lines.append("")
    lines.append("| Secao | Tabela/Figura | Familia | Executavel | Agregado | Linhas | Metrica principal | Max err prim. (%) | Status |")
    lines.append("|---|---|---|---|---|---:|---|---:|---|")
    for row in rows:
        aggregate = Path(row["aggregate_file"])
        label = aggregate.name
        max_err = row["max_err_primary_pct"] or ""
        lines.append(
            f"| `{row['section']}` | `{row['article_ref']}` | `{row['family']}` | `{row['executable']}` | "
            f"[{label}]({label}) | `{row['rows']}` | `{row['primary_metric']}` | `{max_err}` | `{row['status']}` |"
        )
    lines.append("")
    lines.append("## Observacoes")
    lines.append("")
    for row in rows:
        if row["status"] != "present":
            lines.append(
                f"- `{row['article_ref']}` ({row['section']}) esta `{row['status']}` em "
                f"`{Path(row['aggregate_file']).name}`."
            )
    if all(row["status"] == "present" for row in rows):
        lines.append("- Todos os itens catalogados neste indice estao presentes nesta rodada.")
    lines.append("")
    lines.append(f"- CSV maquina-legivel: `[validation_2d_index.csv](validation_2d_index.csv)`")

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a consolidated 2D validation index from existing aggregates.")
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("out/validation"),
        help="Directory containing validation aggregate CSVs.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_2d_index.csv"),
        help="Output CSV path for the consolidated 2D index.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/VALIDATION_2D_INDEX.md"),
        help="Output Markdown path for the consolidated 2D index.",
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
