#!/usr/bin/env python3
"""
Generate a master scientific validation index by combining the consolidated 2D and 3D indexes.

Default outputs:
- out/validation/validation_master_index.csv
- out/validation/VALIDATION_MASTER_INDEX.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _coalesce(*values: str) -> str:
    for value in values:
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _normalize_2d(rows: list[dict[str, str]], *, source_index: Path) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        out.append(
            {
                "dimension": "2D",
                "section": row.get("section", ""),
                "article_ref": row.get("article_ref", ""),
                "family": row.get("family", ""),
                "executables": row.get("executable", ""),
                "source_index_file": _rel_from_root(source_index),
                "aggregate_files": row.get("aggregate_file", ""),
                "coverage_key": row.get("aggregate_cases", ""),
                "primary_metric": row.get("primary_metric", ""),
                "secondary_metric": row.get("secondary_metric", ""),
                "status": row.get("status", ""),
                "data_rows": row.get("rows", ""),
                "summary_rows": "",
                "solvers_present": "",
                "backends_present": "",
                "max_err_primary_pct": row.get("max_err_primary_pct", ""),
                "max_err_secondary_pct": row.get("max_err_secondary_pct", ""),
                "notes": row.get("notes", ""),
            }
        )
    return out


def _normalize_3d(rows: list[dict[str, str]], *, source_index: Path) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in rows:
        aggregate_files = " | ".join(
            part for part in (row.get("modes_file", "").strip(), row.get("summary_file", "").strip()) if part
        )
        out.append(
            {
                "dimension": "3D",
                "section": row.get("section", ""),
                "article_ref": row.get("article_ref", ""),
                "family": row.get("family", ""),
                "executables": row.get("executables", ""),
                "source_index_file": _rel_from_root(source_index),
                "aggregate_files": aggregate_files,
                "coverage_key": row.get("case", ""),
                "primary_metric": row.get("primary_metric", ""),
                "secondary_metric": row.get("secondary_metric", ""),
                "status": row.get("status", ""),
                "data_rows": row.get("mode_rows", ""),
                "summary_rows": row.get("summary_rows", ""),
                "solvers_present": row.get("solvers_present", ""),
                "backends_present": row.get("backends_present", ""),
                "max_err_primary_pct": row.get("max_err_primary_pct", ""),
                "max_err_secondary_pct": row.get("max_err_secondary_pct", ""),
                "notes": row.get("notes", ""),
            }
        )
    return out


def build_master_index(validation_dir: Path) -> tuple[list[dict[str, str]], list[str]]:
    validation_dir = _resolve(validation_dir)
    source_notes: list[str] = []
    out_rows: list[dict[str, str]] = []

    index_2d = validation_dir / "validation_2d_index.csv"
    index_3d = validation_dir / "validation_3d_index.csv"

    rows_2d = _read_rows(index_2d)
    rows_3d = _read_rows(index_3d)

    if not index_2d.exists():
        source_notes.append(f"Fonte ausente: {_rel_from_root(index_2d)}")
    if not index_3d.exists():
        source_notes.append(f"Fonte ausente: {_rel_from_root(index_3d)}")

    out_rows.extend(_normalize_2d(rows_2d, source_index=index_2d))
    out_rows.extend(_normalize_3d(rows_3d, source_index=index_3d))
    out_rows.sort(key=lambda row: (_coalesce(row["section"]), row["dimension"], row["article_ref"]))
    return out_rows, source_notes


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "dimension",
        "section",
        "article_ref",
        "family",
        "executables",
        "source_index_file",
        "aggregate_files",
        "coverage_key",
        "primary_metric",
        "secondary_metric",
        "status",
        "data_rows",
        "summary_rows",
        "solvers_present",
        "backends_present",
        "max_err_primary_pct",
        "max_err_secondary_pct",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_md(path: Path, rows: list[dict[str, str]], source_notes: list[str]) -> None:
    now = dt.datetime.now().isoformat(timespec="seconds")
    lines: list[str] = []
    lines.append("# Indice Mestre de Validacao Cientifica")
    lines.append("")
    lines.append(f"Gerado em: `{now}`")
    lines.append("")
    lines.append("Este indice unifica os indices consolidados 2D e 3D do TP-3485.")
    lines.append("Ele aponta o item do artigo, a familia/executavel, os agregados correspondentes, as metricas principais e o status da rodada.")
    lines.append("")
    lines.append("| Dim | Secao | Tabela/Figura | Familia | Executaveis | Agregados | Dados | Summary | Max err prim. (%) | Status |")
    lines.append("|---|---|---|---|---|---|---:|---:|---:|---|")
    for row in rows:
        lines.append(
            f"| `{row['dimension']}` | `{row['section']}` | `{row['article_ref']}` | `{row['family']}` | "
            f"`{row['executables']}` | `{row['aggregate_files']}` | `{row['data_rows']}` | "
            f"`{row['summary_rows']}` | `{row['max_err_primary_pct']}` | `{row['status']}` |"
        )
    lines.append("")
    lines.append("## Fontes")
    lines.append("")
    seen_sources: list[str] = []
    for row in rows:
        source = row["source_index_file"]
        if source and source not in seen_sources:
            seen_sources.append(source)
            lines.append(f"- `{source}`")
    for note in source_notes:
        lines.append(f"- {note}")
    lines.append("")
    lines.append("## Observacoes")
    lines.append("")
    missing_rows = [row for row in rows if row["status"] != "present"]
    if missing_rows:
        for row in missing_rows:
            lines.append(
                f"- `{row['dimension']} {row['article_ref']}` esta `{row['status']}` "
                f"com agregados `{row['aggregate_files']}`."
            )
    else:
        lines.append("- Todos os itens presentes nos indices 2D e 3D consolidados estao rastreados nesta rodada.")
    lines.append("")
    lines.append("- CSV maquina-legivel: `[validation_master_index.csv](validation_master_index.csv)`")

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a master validation index from the consolidated 2D and 3D indexes.")
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("out/validation"),
        help="Directory containing validation_2d_index.csv and validation_3d_index.csv.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_master_index.csv"),
        help="Output CSV path for the master validation index.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/VALIDATION_MASTER_INDEX.md"),
        help="Output Markdown path for the master validation index.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    validation_dir = _resolve(args.validation_dir)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)

    rows, source_notes = build_master_index(validation_dir)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_csv, rows)
    _write_md(out_md, rows, source_notes)

    print(f"Saved: {_rel_from_root(out_csv)}")
    print(f"Saved: {_rel_from_root(out_md)}")


if __name__ == "__main__":
    main()
