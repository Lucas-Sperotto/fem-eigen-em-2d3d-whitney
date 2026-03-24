#!/usr/bin/env python3
"""
Exporta uma curadoria leve da campanha completa de sweep de malhas,
partindo de `out/sweeps/...` para `docs/results/sweeps/...`.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SRC = ROOT / "out" / "sweeps" / "full_mesh_compare"
DEFAULT_DST = ROOT / "docs" / "results" / "sweeps" / "full_mesh_compare"


ROOT_FILES = [
    "README.md",
    "sweep_runs.csv",
    "sweep_modes.csv",
    "sweep_levels.csv",
]

OPTIONAL_LEVEL_FILES = [
    "validation/validation_2d_22.csv",
    "validation/validation_3d_31_summary.csv",
]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def copy_if_exists(src: Path, dst: Path) -> bool:
    if not src.exists():
        return False
    ensure_parent(dst)
    shutil.copy2(src, dst)
    return True


def summarize_runs(rows: list[dict[str, str]]) -> list[tuple[str, str, int, float | None, float | None]]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row.get("case_id", ""), row.get("backend", ""))].append(row)

    out: list[tuple[str, str, int, float | None, float | None]] = []
    for (case_id, backend), items in sorted(grouped.items()):
        total_vals = [float(r["total_ms"]) for r in items if r.get("total_ms") not in {None, "", "None"}]
        wall_vals = [float(r["wall_ms"]) for r in items if r.get("wall_ms") not in {None, "", "None"}]
        total_mean = sum(total_vals) / len(total_vals) if total_vals else None
        wall_mean = sum(wall_vals) / len(wall_vals) if wall_vals else None
        out.append((case_id, backend, len(items), total_mean, wall_mean))
    return out


def summarize_modes(rows: list[dict[str, str]]) -> list[tuple[str, str, int, float | None, float | None]]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        raw = row.get("err_primary_pct")
        if raw in {None, "", "None"}:
            continue
        grouped[(row.get("case", ""), row.get("backend", ""))].append(abs(float(raw)))

    out: list[tuple[str, str, int, float | None, float | None]] = []
    for (case_id, backend), vals in sorted(grouped.items()):
        if vals:
            out.append((case_id, backend, len(vals), sum(vals) / len(vals), max(vals)))
    return out


def select_levels(level_rows: list[dict[str, str]]) -> dict[str, tuple[str, str]]:
    selected: dict[str, tuple[str, str]] = {}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in level_rows:
        grouped[row.get("backend", "")].append(row)
    for backend, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda r: int(r.get("level", "0") or 0))
        if rows_sorted:
            best = rows_sorted[-1]
            selected[backend] = (best.get("level", ""), best.get("path", ""))
    return selected


def render_readme(dst_root: Path) -> str:
    run_rows = read_csv_rows(dst_root / "sweep_runs.csv")
    mode_rows = read_csv_rows(dst_root / "sweep_modes.csv")
    level_rows = read_csv_rows(dst_root / "sweep_levels.csv")
    selected = select_levels(level_rows)

    lines: list[str] = []
    lines.append("# Sweep Curado")
    lines.append("")
    lines.append("Esta pasta preserva a parte leve e mais comunicativa da campanha completa de malhas.")
    lines.append("O objetivo e registrar tempos, erros e configuracoes principais sem subir toda a massa de VTK, logs e imagens repetidas.")
    lines.append("")
    lines.append("## Como atualizar")
    lines.append("")
    lines.append("Depois de executar a campanha completa em `out/sweeps/full_mesh_compare`, rode:")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/export_curated_sweep_results.py")
    lines.append("```")
    lines.append("")
    lines.append("## Arquivos preservados")
    lines.append("")
    for rel in ROOT_FILES:
        if (dst_root / rel).exists():
            lines.append(f"- [{rel}]({rel})")
    for backend, (level, path_rel) in sorted(selected.items()):
        level_tag = f"level_{int(level):02d}" if level else "ultimo_nivel"
        base = Path(backend) / level_tag / "validation"
        for rel in OPTIONAL_LEVEL_FILES:
            if (dst_root / base / Path(rel).name).exists():
                lines.append(f"- [{base / Path(rel).name}]({base / Path(rel).name})")
    lines.append("")
    lines.append("## Resumo de tempos")
    lines.append("")
    lines.append("| Caso | Backend | Niveis | media total_ms | media wall_ms |")
    lines.append("|---|---|---:|---:|---:|")
    for case_id, backend, n, total_mean, wall_mean in summarize_runs(run_rows):
        total_txt = "-" if total_mean is None else f"{total_mean:.3f}"
        wall_txt = "-" if wall_mean is None else f"{wall_mean:.3f}"
        lines.append(f"| `{case_id}` | `{backend}` | {n} | {total_txt} | {wall_txt} |")
    lines.append("")
    lines.append("## Resumo de erros 2D")
    lines.append("")
    lines.append("| Caso | Backend | amostras | media |err| (%) | pior |err| (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for case_id, backend, n, mean_err, max_err in summarize_modes(mode_rows):
        mean_txt = "-" if mean_err is None else f"{mean_err:.3f}"
        max_txt = "-" if max_err is None else f"{max_err:.3f}"
        lines.append(f"| `{case_id}` | `{backend}` | {n} | {mean_txt} | {max_txt} |")
    lines.append("")
    lines.append("## Niveis preservados por backend")
    lines.append("")
    lines.append("Para manter a curadoria leve, copiamos o `README` e os CSVs globais da campanha inteira, e tambem os CSVs de validacao do ultimo nivel disponivel de cada backend.")
    lines.append("")
    lines.append("| Backend | Nivel curado | Pasta original |")
    lines.append("|---|---:|---|")
    for backend, (level, path_rel) in sorted(selected.items()):
        lines.append(f"| `{backend}` | {level or '-'} | `{path_rel or '-'}` |")
    lines.append("")
    lines.append("## O que ficou de fora")
    lines.append("")
    lines.append("- VTKs de todos os niveis")
    lines.append("- logs completos stdout/stderr")
    lines.append("- imagens repetidas por nivel/modo")
    lines.append("- toda a arvore pesada de `out/sweeps/...`")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Exporta uma curadoria leve da campanha de sweep.")
    parser.add_argument("--src-root", type=Path, default=DEFAULT_SRC)
    parser.add_argument("--dst-root", type=Path, default=DEFAULT_DST)
    args = parser.parse_args()

    src_root = args.src_root.resolve()
    dst_root = args.dst_root.resolve()

    if not src_root.exists():
        raise SystemExit(f"Diretorio de sweep nao encontrado: {src_root}")

    dst_root.mkdir(parents=True, exist_ok=True)

    for rel in ROOT_FILES:
        copy_if_exists(src_root / rel, dst_root / rel)

    level_rows = read_csv_rows(src_root / "sweep_levels.csv")
    selected = select_levels(level_rows)
    for backend, (level, path_rel) in selected.items():
        if not level or not path_rel:
            continue
        level_tag = f"level_{int(level):02d}"
        src_level = src_root / path_rel / "validation"
        dst_level = dst_root / backend / level_tag / "validation"
        for rel in OPTIONAL_LEVEL_FILES:
            src = src_level / Path(rel).name
            dst = dst_level / Path(rel).name
            copy_if_exists(src, dst)

    (dst_root / "README.md").write_text(render_readme(dst_root), encoding="utf-8")
    print(f"Curadoria do sweep salva em: {dst_root}")


if __name__ == "__main__":
    main()
