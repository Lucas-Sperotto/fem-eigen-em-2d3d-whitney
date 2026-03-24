#!/usr/bin/env python3
"""
Exporta uma curadoria leve e versionavel dos resultados gerados em `out/`
para `docs/results/`, evitando subir todo o volume de VTK/PNG/logs.
"""

from __future__ import annotations

import csv
import math
import shutil
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "out"
DEST_DIR = ROOT / "docs" / "results"
CSV_DIR = DEST_DIR / "csv"
FIG_DIR = DEST_DIR / "figs"


@dataclass(frozen=True)
class CopyItem:
    src_rel: str
    dst_rel: str
    label: str


COPY_ITEMS = [
    CopyItem("out/validation/validation_2d_22.csv", "csv/validation_2d_22_legacy.csv",
             "Validacao 2D legado"),
    CopyItem("out/validation/validation_2d_22_closed_form.csv", "csv/validation_2d_22_closed_form.csv",
             "Validacao 2D closed-form"),
    CopyItem("out/validation/validation_3d_31_summary.csv", "csv/validation_3d_31_summary_legacy.csv",
             "Validacao 3D legado"),
    CopyItem("out/validation/validation_3d_31_air_closed_form_summary.csv",
             "csv/validation_3d_31_air_closed_form_summary.csv",
             "Validacao 3D closed-form (caso air)"),
    CopyItem("out/benchmark_test/backend_benchmark_summary.csv", "csv/backend_benchmark_summary.csv",
             "Benchmark de backends"),
    CopyItem("out/img_all/mode_summary.csv", "csv/mode_summary.csv",
             "Resumo de modos 2D"),
    CopyItem("out/img_all/2.1_scalar/rect/te10_rect_sv.png", "figs/2d_scalar_rect_te10.png",
             "Campo escalar 2D retangular TE10"),
    CopyItem("out/img_all/2.1_scalar/circle/te_circle_sv.png", "figs/2d_scalar_circle_te.png",
             "Campo escalar 2D circular TE"),
    CopyItem("out/img_all/2.1_scalar/coax/te_coax_sv.png", "figs/2d_scalar_coax_te.png",
             "Campo escalar 2D coaxial TE"),
    CopyItem("out/img_all/2.2.1_edge/rect/edge_rect_Et.png", "figs/2d_edge_rect_Et.png",
             "Campo edge 2D retangular Et"),
    CopyItem("out/img_all/2.2.1_edge/circle/edge_circle_Et.png", "figs/2d_edge_circle_Et.png",
             "Campo edge 2D circular Et"),
    CopyItem("out/img_all/2.2.1_edge/coax/edge_coax_Et.png", "figs/2d_edge_coax_Et.png",
             "Campo edge 2D coaxial Et"),
    CopyItem("out/img_all/validation_2d_22/2.2.2/2_2_2_mixed_rect_cutoff.png",
             "figs/2d_22_mixed_rect_cutoff.png",
             "Validacao 2.2.2 retangular"),
    CopyItem("out/img_all/validation_2d_22/2.2.2/2_2_2_mixed_circle_coax_snapshots.png",
             "figs/2d_22_mixed_circle_coax.png",
             "Validacao 2.2.2 circular/coaxial"),
    CopyItem("out/img_all/validation_2d_22/2.2.3/2_2_3_table8_k0L.png",
             "figs/2d_23_table8_k0L.png",
             "Validacao 2.2.3 tabela 8"),
    CopyItem("out/img_all/validation_2d_22/2.2.3/2_2_3_table8_error_pct.png",
             "figs/2d_23_table8_error_pct.png",
             "Erro relativo 2.2.3 tabela 8"),
    CopyItem("out/img_all/validation_2d_22/2.2.4/2_2_4_table9_beta_over_k0.png",
             "figs/2d_24_table9_beta_over_k0.png",
             "Validacao 2.2.4 tabela 9"),
    CopyItem("out/img_all/validation_2d_22/2.2.4/2_2_4_table10_fem_branches.png",
             "figs/2d_24_table10_fem_branches.png",
             "Validacao 2.2.4 tabela 10"),
]


def ensure_dirs() -> None:
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def copy_items() -> list[CopyItem]:
    copied: list[CopyItem] = []
    for item in COPY_ITEMS:
        src = ROOT / item.src_rel
        dst = DEST_DIR / item.dst_rel
        if not src.exists():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied.append(item)
    return copied


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt_pct(value: float | None) -> str:
    if value is None or math.isnan(value):
        return "-"
    return f"{value:.3f}"


def float_or_none(raw: str | None) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def summarize_2d(rows: list[dict[str, str]]) -> list[tuple[str, str, int, float | None, float | None]]:
    agg: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        section = row.get("section", "")
        case = row.get("case", "")
        key = (section, case)
        stats = agg.setdefault(key, {"n": 0, "max_primary": None, "max_secondary": None})
        stats["n"] = int(stats["n"]) + 1
        p = float_or_none(row.get("err_primary_pct"))
        s = float_or_none(row.get("err_secondary_pct"))
        if p is not None:
            cur = stats["max_primary"]
            stats["max_primary"] = max(abs(p), cur if isinstance(cur, float) else 0.0)
        if s is not None:
            cur = stats["max_secondary"]
            stats["max_secondary"] = max(abs(s), cur if isinstance(cur, float) else 0.0)
    out: list[tuple[str, str, int, float | None, float | None]] = []
    for (section, case), stats in sorted(agg.items()):
        out.append((
            section,
            case,
            int(stats["n"]),
            stats["max_primary"] if isinstance(stats["max_primary"], float) else None,
            stats["max_secondary"] if isinstance(stats["max_secondary"], float) else None,
        ))
    return out


def summarize_3d(rows: list[dict[str, str]]) -> list[tuple[str, str, int, float | None, float | None]]:
    agg: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        solver = row.get("solver", "")
        case = row.get("case", "")
        key = (solver, case)
        stats = agg.setdefault(key, {"n": 0, "best_mean_ana": None, "best_mean_ref": None})
        stats["n"] = int(stats["n"]) + 1
        mean_ana = float_or_none(row.get("mean_err_ana_pct"))
        mean_ref = float_or_none(row.get("mean_err_ref_pct"))
        if mean_ana is not None:
            cur = stats["best_mean_ana"]
            if cur is None or mean_ana < cur:
                stats["best_mean_ana"] = mean_ana
        if mean_ref is not None:
            cur = stats["best_mean_ref"]
            if cur is None or mean_ref < cur:
                stats["best_mean_ref"] = mean_ref
    out: list[tuple[str, str, int, float | None, float | None]] = []
    for (solver, case), stats in sorted(agg.items()):
        out.append((
            solver,
            case,
            int(stats["n"]),
            stats["best_mean_ana"] if isinstance(stats["best_mean_ana"], float) else None,
            stats["best_mean_ref"] if isinstance(stats["best_mean_ref"], float) else None,
        ))
    return out


def summarize_benchmark(rows: list[dict[str, str]]) -> list[tuple[str, float | None, float | None, str]]:
    grouped: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        case_id = row.get("case_id", "")
        backend = row.get("backend", "")
        total = float_or_none(row.get("total_ms_mean"))
        if case_id and backend and total is not None:
            grouped[case_id][backend] = total
    out: list[tuple[str, float | None, float | None, str]] = []
    for case_id in sorted(grouped):
        gauss = grouped[case_id].get("gauss")
        closed = grouped[case_id].get("closed-form")
        faster = "-"
        if gauss is not None and closed is not None:
            faster = "closed-form" if closed < gauss else "gauss"
        out.append((case_id, gauss, closed, faster))
    return out


def render_markdown(copied: list[CopyItem]) -> str:
    rows_2d_legacy = read_csv_rows(DEST_DIR / "csv" / "validation_2d_22_legacy.csv")
    rows_2d_cf = read_csv_rows(DEST_DIR / "csv" / "validation_2d_22_closed_form.csv")
    rows_3d_legacy = read_csv_rows(DEST_DIR / "csv" / "validation_3d_31_summary_legacy.csv")
    rows_3d_cf = read_csv_rows(DEST_DIR / "csv" / "validation_3d_31_air_closed_form_summary.csv")
    bench_rows = read_csv_rows(DEST_DIR / "csv" / "backend_benchmark_summary.csv")

    lines: list[str] = []
    lines.append("# Resultados Curados")
    lines.append("")
    lines.append("Esta pasta contem apenas os resultados leves e estaveis que valem a pena versionar no Git.")
    lines.append("Os artefatos completos de `out/` continuam uteis localmente, mas nao entram todos no repositorio.")
    lines.append("")
    lines.append("## O que foi preservado")
    lines.append("")
    lines.append("- CSVs-resumo de validacao 2D, 3D e benchmark")
    lines.append("- resumo de modos 2D")
    lines.append("- figuras representativas dos campos 2D")
    lines.append("- figuras principais de validacao das secoes 2.2.2, 2.2.3 e 2.2.4")
    lines.append("")
    lines.append("## Arquivos incluidos")
    lines.append("")
    for item in copied:
        lines.append(f"- [{item.dst_rel}]({item.dst_rel}): {item.label}")
    lines.append("")

    lines.append("## Resumo 2D")
    lines.append("")
    lines.append("### Validacao 2.2.x - legado")
    lines.append("")
    lines.append("| Secao | Caso | Linhas | Max |err primary| (%) | Max |err secondary| (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for section, case, n, p, s in summarize_2d(rows_2d_legacy):
        lines.append(f"| `{section}` | `{case}` | {n} | {fmt_pct(p)} | {fmt_pct(s)} |")
    lines.append("")

    lines.append("### Validacao 2.2.x - closed-form")
    lines.append("")
    lines.append("| Secao | Caso | Linhas | Max |err primary| (%) | Max |err secondary| (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for section, case, n, p, s in summarize_2d(rows_2d_cf):
        lines.append(f"| `{section}` | `{case}` | {n} | {fmt_pct(p)} | {fmt_pct(s)} |")
    lines.append("")

    lines.append("### Figuras 2D escolhidas")
    lines.append("")
    lines.append("- Escalar retangular: [figs/2d_scalar_rect_te10.png](figs/2d_scalar_rect_te10.png)")
    lines.append("- Escalar circular: [figs/2d_scalar_circle_te.png](figs/2d_scalar_circle_te.png)")
    lines.append("- Escalar coaxial: [figs/2d_scalar_coax_te.png](figs/2d_scalar_coax_te.png)")
    lines.append("- Edge retangular: [figs/2d_edge_rect_Et.png](figs/2d_edge_rect_Et.png)")
    lines.append("- Edge circular: [figs/2d_edge_circle_Et.png](figs/2d_edge_circle_Et.png)")
    lines.append("- Edge coaxial: [figs/2d_edge_coax_Et.png](figs/2d_edge_coax_Et.png)")
    lines.append("- Secao 2.2.2: [figs/2d_22_mixed_rect_cutoff.png](figs/2d_22_mixed_rect_cutoff.png), [figs/2d_22_mixed_circle_coax.png](figs/2d_22_mixed_circle_coax.png)")
    lines.append("- Secao 2.2.3: [figs/2d_23_table8_k0L.png](figs/2d_23_table8_k0L.png), [figs/2d_23_table8_error_pct.png](figs/2d_23_table8_error_pct.png)")
    lines.append("- Secao 2.2.4: [figs/2d_24_table9_beta_over_k0.png](figs/2d_24_table9_beta_over_k0.png), [figs/2d_24_table10_fem_branches.png](figs/2d_24_table10_fem_branches.png)")
    lines.append("")

    lines.append("## Resumo 3D")
    lines.append("")
    lines.append("### Melhor erro medio por solver/caso - legado")
    lines.append("")
    lines.append("| Solver | Caso | Configs | Melhor mean err ana (%) | Melhor mean err ref (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for solver, case, n, mean_ana, mean_ref in summarize_3d(rows_3d_legacy):
        lines.append(f"| `{solver}` | `{case}` | {n} | {fmt_pct(mean_ana)} | {fmt_pct(mean_ref)} |")
    lines.append("")

    if rows_3d_cf:
        lines.append("### Closed-form 3D preservado")
        lines.append("")
        lines.append("Nesta curadoria, foi preservado o resumo `closed-form` do caso `air`, suficiente para documentar a trilha de comparacao sem subir todos os artefatos intermediarios.")
        lines.append("")
        lines.append("- [csv/validation_3d_31_air_closed_form_summary.csv](csv/validation_3d_31_air_closed_form_summary.csv)")
        lines.append("")

    lines.append("## Benchmark de backends")
    lines.append("")
    lines.append("| Caso | total_ms_mean gauss | total_ms_mean closed-form | Mais rapido |")
    lines.append("|---|---:|---:|---|")
    for case_id, gauss, closed, faster in summarize_benchmark(bench_rows):
        g = "-" if gauss is None else f"{gauss:.3f}"
        c = "-" if closed is None else f"{closed:.3f}"
        lines.append(f"| `{case_id}` | {g} | {c} | `{faster}` |")
    lines.append("")

    lines.append("## O que ficou de fora do Git")
    lines.append("")
    lines.append("- pasta `out/` completa")
    lines.append("- todos os `.vtk`")
    lines.append("- logs integrais de execucao")
    lines.append("- imagens repetidas por rank/modo")
    lines.append("- CSVs detalhados demais para leitura humana")
    lines.append("")
    lines.append("A ideia aqui e simples: o repositório guarda o que comunica bem os resultados; os artefatos pesados continuam regeneraveis localmente.")
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    ensure_dirs()
    copied = copy_items()
    readme = render_markdown(copied)
    (DEST_DIR / "README.md").write_text(readme + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
