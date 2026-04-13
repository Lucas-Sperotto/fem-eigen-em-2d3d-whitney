#!/usr/bin/env python3
"""
Generate a consolidated Markdown report for outputs under out/.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ROOT = Path(__file__).resolve().parents[1]


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _md_link(path: Path, label: Optional[str] = None) -> str:
    lbl = label or _rel_from_root(path)
    return f"[{lbl}]({_rel_from_root(path)})"


def _to_float(s: str) -> Optional[float]:
    if s is None:
        return None
    t = str(s).strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def _read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _count_by_case_dir(base: Path, suffix: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not base.exists():
        return out
    for p in sorted(base.rglob(f"*{suffix}")):
        key = str(p.parent.relative_to(base))
        out[key] = out.get(key, 0) + 1
    return out


def _count_2d_vtks(out_dir: Path) -> Dict[str, int]:
    legacy_root = out_dir / "2d"
    if legacy_root.exists():
        return _count_by_case_dir(legacy_root, ".vtk")

    out: Dict[str, int] = {}
    for family_name, case_name in (("helm10", "2.1_scalar"), ("helmvec", "2.2.1_edge")):
        family_root = out_dir / family_name
        if not family_root.exists():
            continue
        for vtk_path in sorted(family_root.glob("*/vtk/*.vtk")):
            geometry = vtk_path.parent.parent.name
            key = f"{case_name}/{geometry}"
            out[key] = out.get(key, 0) + 1
    return out


def _analyze_log(path: Path) -> Dict[str, object]:
    info: Dict[str, object] = {
        "exists": path.exists(),
        "size": 0,
        "line_count": 0,
        "keywords": {},
        "hits": [],
        "config": {},
        "tail": [],
    }
    if not path.exists():
        return info

    txt = path.read_text(encoding="utf-8", errors="ignore")
    lines = txt.splitlines()
    info["size"] = path.stat().st_size
    info["line_count"] = len(lines)

    # Use regex with word boundaries to avoid false positives from filenames
    # like "error_pct.png" in generated artifacts.
    keyword_patterns = {
        "error": re.compile(r"\berror\b", flags=re.IGNORECASE),
        "failed": re.compile(r"\bfailed\b", flags=re.IGNORECASE),
        "exception": re.compile(r"\bexception\b", flags=re.IGNORECASE),
        "traceback": re.compile(r"\btraceback\b", flags=re.IGNORECASE),
        "segmentation fault": re.compile(r"\bsegmentation\s+fault\b", flags=re.IGNORECASE),
    }

    counts = {k: 0 for k in keyword_patterns}
    hits = []
    for i, ln in enumerate(lines, start=1):
        matched = False
        for key, pat in keyword_patterns.items():
            if pat.search(ln):
                counts[key] += 1
                matched = True
        if matched:
            hits.append((i, ln))
    info["keywords"] = counts
    info["hits"] = hits

    cfg = {}
    for ln in lines:
        m = re.match(r"^\[run_all\]\s+([A-Z_]+)=(.+)$", ln.strip())
        if m:
            cfg[m.group(1)] = m.group(2)
    info["config"] = cfg
    info["tail"] = lines[-20:]
    return info


def _mode_paths(out_dir: Path, row: Dict[str, str]) -> Tuple[Optional[Path], Optional[Path]]:
    geometry = row.get("geometry", "").strip()
    formulation = row.get("formulation", "").strip()
    pol = row.get("polarization", "").strip().lower()
    field = row.get("field", "").strip()
    rank = int(float(row.get("mode_rank", "0") or 0))
    m = int(float(row.get("m", "0") or 0))
    idx2_name = row.get("index2_name", "").strip()
    idx2 = row.get("n", "").strip() if idx2_name == "n" else row.get("p", "").strip()
    idx2 = str(int(float(idx2))) if idx2 else "0"

    if formulation == "scalar":
        case = "2.1_scalar"
        fname = f"{pol}_{geometry}_m{m}_{idx2_name}{idx2}_rank{rank:02d}_sv"
        vtk_candidates = [
            out_dir / "helm10" / geometry / "vtk" / f"{fname}.vtk",
            out_dir / "2d" / case / geometry / f"{fname}.vtk",
        ]
    elif formulation == "edge":
        case = "2.2.1_edge"
        fname = f"edge_{geometry}_{pol}_m{m}_{idx2_name}{idx2}_rank{rank:02d}_{field}"
        vtk_candidates = [
            out_dir / "helmvec" / geometry / "vtk" / f"{fname}.vtk",
            out_dir / "2d" / case / geometry / f"{fname}.vtk",
        ]
    else:
        return None, None

    vtk = next((path for path in vtk_candidates if path.exists()), None)
    png = out_dir / "img_all" / case / geometry / f"{fname}.png"
    return vtk, (png if png.exists() else None)


def _top_mode_errors(mode_rows: List[Dict[str, str]], top_n: int = 12) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for r in mode_rows:
        err = _to_float(r.get("rel_error_pct", ""))
        if err is None:
            continue
        rr = dict(r)
        rr["_abs_err"] = abs(err)
        out.append(rr)
    out.sort(key=lambda r: float(r["_abs_err"]), reverse=True)
    return out[:top_n]


def _group_2d22_case_stats(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        groups[(r.get("backend", "legacy"), r.get("section", ""), r.get("case", ""))].append(r)

    out = []
    for (backend, sec, case), gr in sorted(groups.items()):
        e1 = [abs(_to_float(r.get("err_primary_pct", "")) or 0.0) for r in gr if _to_float(r.get("err_primary_pct", "")) is not None]
        e2 = [abs(_to_float(r.get("err_secondary_pct", "")) or 0.0) for r in gr if _to_float(r.get("err_secondary_pct", "")) is not None]
        out.append(
            {
                "backend": backend,
                "section": sec,
                "case": case,
                "rows": len(gr),
                "max_err_primary_pct": max(e1) if e1 else None,
                "max_err_secondary_pct": max(e2) if e2 else None,
            }
        )
    return out


def _max_3d_modes(mode_rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in mode_rows:
        groups[(r.get("solver", ""), r.get("backend", "legacy"), r.get("case", ""))].append(r)

    out = []
    for (solver, backend, case), gr in sorted(groups.items()):
        e_ana = [abs(_to_float(r.get("err_ana_pct", "")) or 0.0) for r in gr]
        e_ref = [abs(_to_float(r.get("err_ref_pct", "")) or 0.0) for r in gr]
        out.append(
            {
                "solver": solver,
                "backend": backend,
                "case": case,
                "rows": len(gr),
                "max_err_ana_pct": max(e_ana) if e_ana else None,
                "max_err_ref_pct": max(e_ref) if e_ref else None,
            }
        )
    return out


def _fmt(x: Optional[float], nd: int = 4) -> str:
    if x is None:
        return "-"
    if not math.isfinite(x):
        return "nan"
    return f"{x:.{nd}f}"


def _benchmark_speedups(summary_rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    grouped: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
    for row in summary_rows:
        grouped[row.get("case_id", "")][row.get("backend", "")] = row

    out: List[Dict[str, object]] = []
    for case_id, backends in sorted(grouped.items()):
        if "gauss" not in backends or "closed-form" not in backends:
            continue
        gauss = backends["gauss"]
        closed = backends["closed-form"]
        wall_gauss = _to_float(gauss.get("wall_ms_mean", ""))
        wall_closed = _to_float(closed.get("wall_ms_mean", ""))
        speedup = None
        if wall_gauss is not None and wall_closed is not None and wall_closed > 0.0:
            speedup = wall_gauss / wall_closed
        out.append(
            {
                "section": gauss.get("section", ""),
                "case_id": case_id,
                "wall_gauss": wall_gauss,
                "wall_closed": wall_closed,
                "speedup": speedup,
                "assembly_gauss": _to_float(gauss.get("assembly_ms_mean", "")),
                "assembly_closed": _to_float(closed.get("assembly_ms_mean", "")),
            }
        )
    return out


def _equation_trace_rows() -> List[Dict[str, str]]:
    return [
        {
            "program": "HELM10",
            "section": "2.1",
            "local_eqs": "30, 33",
            "global_eq": "43",
            "module": "src/helm10",
            "explicit": "src/explicit/tri2d_scalar_explicit.hpp",
            "assembly": "src/core/helm10_scalar_system.cpp",
            "function": "tp3485::build_eq43_helm10_system(...)",
        },
        {
            "program": "HELMVEC",
            "section": "2.2.1",
            "local_eqs": "66, 67",
            "global_eq": "65",
            "module": "src/helmvec",
            "explicit": "src/explicit/tri2d_edge_explicit.hpp",
            "assembly": "src/edge/edge_assembly.cpp",
            "function": "tp3485::build_eq65_helmvec_system(...)",
        },
        {
            "program": "HELMVEC1",
            "section": "2.2.2",
            "local_eqs": "66, 67, 30, 33",
            "global_eq": "92",
            "module": "src/helmvec1",
            "explicit": "src/explicit/tri2d_edge_explicit.hpp + src/explicit/tri2d_scalar_explicit.hpp",
            "assembly": "src/helmvec1/helmvec1_mixed_system.cpp",
            "function": "tp3485::build_eq92_helmvec1_system_E/H(...)",
        },
        {
            "program": "HELMVEC2",
            "section": "2.2.3",
            "local_eqs": "120-125",
            "global_eq": "119",
            "module": "src/helmvec2",
            "explicit": "src/explicit/tri2d_coupled_explicit.hpp",
            "assembly": "src/helmvec2/helmvec2_coupled_system.cpp",
            "function": "tp3485::build_eq119_helmvec2_system_E(...)",
        },
        {
            "program": "HELMVEC3",
            "section": "2.2.4",
            "local_eqs": "137-142",
            "global_eq": "136",
            "module": "src/helmvec3",
            "explicit": "src/explicit/tri2d_coupled_explicit.hpp",
            "assembly": "src/helmvec2/helmvec2_coupled_system.cpp",
            "function": "tp3485::build_eq136_helmvec3_system_E(...)",
        },
        {
            "program": "FEM3D0/FEM3D1",
            "section": "3.1",
            "local_eqs": "181, 182",
            "global_eq": "178",
            "module": "src/fem3d0 + src/fem3d1",
            "explicit": "src/explicit/tet3d_edge_explicit.hpp",
            "assembly": "src/edge3d/edge3d_assembly.cpp (build_eq178_local_tet_blocks -> assemble_eq178_global_generic -> dense/sparse)",
            "function": "tp3485::build_eq178_fem3d_system_dense/sparse(...)",
        },
    ]


def _write_report(
    report_path: Path,
    out_dir: Path,
    mode_rows: List[Dict[str, str]],
    v2d_rows: List[Dict[str, str]],
    v3m_rows: List[Dict[str, str]],
    v3s_rows: List[Dict[str, str]],
    bench_detail_rows: List[Dict[str, str]],
    bench_summary_rows: List[Dict[str, str]],
    log_info: Dict[str, object],
) -> None:
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    vtk_counts = _count_2d_vtks(out_dir)
    png_counts = _count_by_case_dir(out_dir / "img_all", ".png")

    mode_cov = Counter((r["geometry"], r["formulation"], r["polarization"]) for r in mode_rows)
    top_mode = _top_mode_errors(mode_rows, top_n=12)
    c2 = _group_2d22_case_stats(v2d_rows)
    c3 = _max_3d_modes(v3m_rows)
    cbench = _benchmark_speedups(bench_summary_rows)

    v22_img = out_dir / "img_all" / "validation_2d_22"
    v22_imgs = sorted(v22_img.rglob("*.png")) if v22_img.exists() else []

    lines: List[str] = []
    lines.append("# Relatorio Consolidado de Resultados")
    lines.append("")
    lines.append(f"Gerado em: `{now}`")
    lines.append("")
    lines.append("## Links Principais")
    lines.append("")
    lines.append(f"- Log principal: {_md_link(out_dir / 'run_all.log')}")
    lines.append(f"- CSV modos 2D (campos): {_md_link(out_dir / 'img_all' / 'mode_summary.csv')}")
    lines.append(f"- CSV validacao 2D (2.1): {_md_link(out_dir / 'validation' / 'validation_2d_21.csv')}")
    lines.append(f"- CSV validacao 2D (2.2.1): {_md_link(out_dir / 'validation' / 'validation_2d_221.csv')}")
    lines.append(f"- CSV validacao 2D (2.2.2 Tabela 6): {_md_link(out_dir / 'validation' / 'validation_2d_222_table6.csv')}")
    lines.append(f"- CSV validacao 2D (2.2.2 Tabela 7): {_md_link(out_dir / 'validation' / 'validation_2d_222_table7.csv')}")
    lines.append(f"- CSV validacao 2D (2.2.4 Tabela 10): {_md_link(out_dir / 'validation' / 'validation_2d_224_table10.csv')}")
    lines.append(f"- CSV validacao 2D (2.2.x): {_md_link(out_dir / 'validation' / 'validation_2d_22.csv')}")
    lines.append(f"- CSV indice consolidado 2D: {_md_link(out_dir / 'validation' / 'validation_2d_index.csv')}")
    lines.append(f"- Markdown indice consolidado 2D: {_md_link(out_dir / 'validation' / 'VALIDATION_2D_INDEX.md')}")
    lines.append(f"- CSV validacao 3D (modes): {_md_link(out_dir / 'validation' / 'validation_3d_31_modes.csv')}")
    lines.append(f"- CSV validacao 3D (summary): {_md_link(out_dir / 'validation' / 'validation_3d_31_summary.csv')}")
    lines.append(f"- CSV indice consolidado 3D: {_md_link(out_dir / 'validation' / 'validation_3d_index.csv')}")
    lines.append(f"- Markdown indice consolidado 3D: {_md_link(out_dir / 'validation' / 'VALIDATION_3D_INDEX.md')}")
    lines.append(f"- CSV indice mestre: {_md_link(out_dir / 'validation' / 'validation_master_index.csv')}")
    lines.append(f"- Markdown indice mestre: {_md_link(out_dir / 'validation' / 'VALIDATION_MASTER_INDEX.md')}")
    lines.append(f"- Politica de thresholds: {_md_link(ROOT / 'docs' / 'validation_thresholds.csv')}")
    lines.append(f"- CSV veredito cientifico: {_md_link(out_dir / 'validation' / 'validation_verdict.csv')}")
    lines.append(f"- Markdown veredito cientifico: {_md_link(out_dir / 'validation' / 'VALIDATION_VERDICT.md')}")
    lines.append(f"- Imagens de campos 2D: {_md_link(out_dir / 'img_all')}")
    lines.append(f"- Imagens de validacao 2.2.x: {_md_link(v22_img)}")
    bench_detail = out_dir / "benchmark" / "backend_benchmark_detail.csv"
    bench_summary = out_dir / "benchmark" / "backend_benchmark_summary.csv"
    bench_report = out_dir / "benchmark" / "BACKEND_BENCHMARK.md"
    if bench_detail.exists():
        lines.append(f"- Benchmark detalhado: {_md_link(bench_detail)}")
    if bench_summary.exists():
        lines.append(f"- Benchmark agregado: {_md_link(bench_summary)}")
    if bench_report.exists():
        lines.append(f"- Benchmark Markdown: {_md_link(bench_report)}")
    lines.append("")

    lines.append("## Rastreabilidade das Equacoes")
    lines.append("")
    lines.append("Esta secao resume a trilha didatica principal do repositorio:")
    lines.append("equacoes locais closed-form -> entrada didatica -> montagem interna.")
    lines.append("")
    lines.append("| Programa | Secao | Eq. locais | Eq. global | Modulo | Closed-form | Entrada didatica | Montagem interna |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for row in _equation_trace_rows():
        lines.append(
            f"| `{row['program']}` | `{row['section']}` | `{row['local_eqs']}` | `{row['global_eq']}` | "
            f"`{row['module']}` | `{row['explicit']}` | `{row['function']}` | `{row['assembly']}` |"
        )
    lines.append("")

    lines.append("## Depuracao, Backends e Reproducao")
    lines.append("")
    lines.append("Fluxos recomendados para inspecao didatica e comparacao numerica:")
    lines.append("")
    lines.append("- `--backend closed-form`: usa os helpers ligados diretamente as equacoes do artigo e e o fluxo principal do repositorio.")
    lines.append("- `--backend gauss`: usa a montagem por quadratura/cubatura como verificacao auxiliar.")
    lines.append("- `--debug-local-blocks`: imprime os blocos locais do primeiro elemento.")
    lines.append("- `--debug-candidates`: imprime as primeiras raizes/candidatos antes do matching final.")
    lines.append("- `--show-output` (scripts Python): ecoa a saida bruta dos executaveis durante a validacao.")
    lines.append("- `--show-validation-output` (`build_and_run_all.sh`): repassa esse comportamento para o pipeline completo.")
    lines.append("")
    lines.append("```bash")
    lines.append("./build/helm10_rect --ar-m 1.0 --nx 10 --ny 20 --nmodos 10 --backend closed-form --debug-local-blocks")
    lines.append("./build/edge_rect --nx 10 --ny 20 --nmodos 10 --backend closed-form --debug-candidates")
    lines.append("./build/mixed_rect --nx 10 --ny 20 --backend closed-form --debug-local-blocks --debug-candidates")
    lines.append("python3 scripts/validate_2d_22.py --backend closed-form --show-output --debug-local-blocks")
    lines.append("python3 scripts/validate_3d_31.py --backend closed-form --show-output --debug-candidates")
    lines.append("scripts/build_and_run_all.sh --backend closed-form --show-validation-output --debug-local-blocks")
    lines.append("```")
    lines.append("")

    lines.append("## Inventario de Arquivos")
    lines.append("")
    n_vtk = sum(vtk_counts.values())
    n_png = sum(png_counts.values())
    n_csv = len(list(out_dir.rglob("*.csv")))
    n_log = len(list(out_dir.rglob("*.log")))
    lines.append(f"- Total VTK: `{n_vtk}`")
    lines.append(f"- Total PNG: `{n_png}`")
    lines.append(f"- Total CSV: `{n_csv}`")
    lines.append(f"- Total LOG: `{n_log}`")
    lines.append("")
    lines.append("### VTK por caso")
    lines.append("")
    lines.append("| Caso | Qtde |")
    lines.append("|---|---:|")
    for k in sorted(vtk_counts):
        lines.append(f"| `{k}` | {vtk_counts[k]} |")
    lines.append("")
    lines.append("### PNG por caso")
    lines.append("")
    lines.append("| Caso | Qtde |")
    lines.append("|---|---:|")
    for k in sorted(png_counts):
        lines.append(f"| `{k}` | {png_counts[k]} |")
    lines.append("")

    lines.append("## Resumo do Log")
    lines.append("")
    if not log_info.get("exists"):
        lines.append("Log nao encontrado.")
    else:
        lines.append(f"- Arquivo: {_md_link(out_dir / 'run_all.log')}")
        lines.append(f"- Tamanho: `{log_info['size']} bytes`")
        lines.append(f"- Linhas: `{log_info['line_count']}`")
        cfg = log_info.get("config", {})
        if cfg:
            lines.append("")
            lines.append("### Configuracao capturada")
            lines.append("")
            for k in sorted(cfg.keys()):
                lines.append(f"- `{k}` = `{cfg[k]}`")
        lines.append("")
        lines.append("### Palavras-chave de erro no log")
        lines.append("")
        lines.append("| Palavra | Ocorrencias |")
        lines.append("|---|---:|")
        kw = log_info.get("keywords", {})
        for k in ["error", "failed", "exception", "traceback", "segmentation fault"]:
            lines.append(f"| `{k}` | {kw.get(k, 0)} |")

        hits = log_info.get("hits", [])
        if hits:
            lines.append("")
            lines.append("### Linhas com alerta")
            lines.append("")
            for ln, txt in hits[:40]:
                lines.append(f"- linha `{ln}`: `{txt}`")
        else:
            lines.append("")
            lines.append("Nenhuma linha suspeita detectada por palavras-chave.")
        lines.append("")
        lines.append("### Trecho final do log")
        lines.append("")
        lines.append("```text")
        for ln in log_info.get("tail", []):
            lines.append(ln)
        lines.append("```")
    lines.append("")

    lines.append("## Tabela de Modos 2D (mode_summary.csv)")
    lines.append("")
    lines.append(f"- Linhas: `{len(mode_rows)}`")
    lines.append("")
    lines.append("### Cobertura por geometria/formulacao/polarizacao")
    lines.append("")
    lines.append("| Geometria | Formulacao | Polarizacao | Linhas |")
    lines.append("|---|---|---|---:|")
    for (g, f, p), n in sorted(mode_cov.items()):
        lines.append(f"| `{g}` | `{f}` | `{p}` | {n} |")
    lines.append("")
    lines.append("### Top 12 maiores erros relativos (|rel_error_pct|)")
    lines.append("")
    lines.append("| # | geom | form | pol | rank | m | idx | kc_fem | err(%) | VTK | PNG |")
    lines.append("|---:|---|---|---|---:|---:|---|---:|---:|---|---|")
    for i, r in enumerate(top_mode, start=1):
        idx2 = f"n={r.get('n')}" if (r.get("index2_name") == "n") else f"p={r.get('p')}"
        vtk, png = _mode_paths(out_dir, r)
        lines.append(
            f"| {i} | `{r.get('geometry')}` | `{r.get('formulation')}` | `{r.get('polarization')}` | "
            f"{int(float(r.get('mode_rank','0')))} | {int(float(r.get('m','0')))} | `{idx2}` | "
            f"{_fmt(_to_float(r.get('kc_fem','')), 5)} | {_fmt(abs(_to_float(r.get('rel_error_pct','')) or 0.0), 4)} | "
            f"{_md_link(vtk, 'vtk') if vtk else '-'} | {_md_link(png, 'png') if png else '-'} |"
        )
    lines.append("")

    lines.append("## Validacao 2D (2.2.2 / 2.2.3 / 2.2.4)")
    lines.append("")
    lines.append(f"- Linhas em CSV: `{len(v2d_rows)}`")
    lines.append("")
    lines.append("### Estatisticas por caso")
    lines.append("")
    lines.append("| Backend | Secao | Caso | Linhas | Max err primary (%) | Max err secondary (%) |")
    lines.append("|---|---|---|---:|---:|---:|")
    for r in c2:
        lines.append(
            f"| `{r['backend']}` | `{r['section']}` | `{r['case']}` | {r['rows']} | {_fmt(r['max_err_primary_pct'])} | {_fmt(r['max_err_secondary_pct'])} |"
        )
    lines.append("")
    if v22_imgs:
        lines.append("### Imagens de validacao 2.2.x")
        lines.append("")
        for p in v22_imgs:
            lines.append(f"- {_md_link(p)}")
        lines.append("")

    lines.append("## Validacao 3D (Secao 3.1)")
    lines.append("")
    lines.append(f"- Linhas em modes CSV: `{len(v3m_rows)}`")
    lines.append(f"- Linhas em summary CSV: `{len(v3s_rows)}`")
    lines.append("")
    lines.append("### Maximos por solver/caso (a partir de validation_3d_31_modes.csv)")
    lines.append("")
    lines.append("| Solver | Backend | Caso | Linhas | Max err ana (%) | Max err ref (%) |")
    lines.append("|---|---|---|---:|---:|---:|")
    for r in c3:
        lines.append(
            f"| `{r['solver']}` | `{r['backend']}` | `{r['case']}` | {r['rows']} | {_fmt(r['max_err_ana_pct'])} | {_fmt(r['max_err_ref_pct'])} |"
        )
    lines.append("")

    lines.append("### Summary detalhado (validation_3d_31_summary.csv)")
    lines.append("")
    lines.append("| Solver | Backend | Caso | nx | ny | nz | n_modes | max_err_ana | mean_err_ana | max_err_ref | mean_err_ref |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in v3s_rows:
        lines.append(
            f"| `{r.get('solver','')}` | `{r.get('backend','legacy')}` | `{r.get('case','')}` | {r.get('nx','')} | {r.get('ny','')} | {r.get('nz','')} | "
            f"{r.get('n_modes','')} | {_fmt(_to_float(r.get('max_err_ana_pct','')))} | {_fmt(_to_float(r.get('mean_err_ana_pct','')))} | "
            f"{_fmt(_to_float(r.get('max_err_ref_pct','')))} | {_fmt(_to_float(r.get('mean_err_ref_pct','')))} |"
        )
    lines.append("")

    if bench_summary_rows:
        lines.append("## Benchmark de Backends")
        lines.append("")
        lines.append(f"- Linhas em benchmark detalhado: `{len(bench_detail_rows)}`")
        lines.append(f"- Linhas em benchmark agregado: `{len(bench_summary_rows)}`")
        lines.append("")
        lines.append("### Speedup closed-form vs gauss")
        lines.append("")
        lines.append("| Secao | Caso | Wall gauss (ms) | Wall closed-form (ms) | Speedup | Assembly gauss (ms) | Assembly closed-form (ms) |")
        lines.append("|---|---|---:|---:|---:|---:|---:|")
        for r in cbench:
            lines.append(
                f"| `{r['section']}` | `{r['case_id']}` | {_fmt(r['wall_gauss'], 3)} | {_fmt(r['wall_closed'], 3)} | "
                f"{_fmt(r['speedup'], 3)}x | {_fmt(r['assembly_gauss'], 3)} | {_fmt(r['assembly_closed'], 3)} |"
            )
        lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines).strip() + "\n", encoding="utf-8")
    print(f"Saved: {report_path}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Generate consolidated Markdown report for out/ artifacts.")
    ap.add_argument("--out-dir", type=Path, default=Path("out"), help="Output root directory.")
    ap.add_argument("--report", type=Path, default=Path("out/RESULTS_REPORT.md"), help="Markdown report output path.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = _resolve(args.out_dir)
    report = _resolve(args.report)

    mode_csv = out_dir / "img_all" / "mode_summary.csv"
    v2d_csv = out_dir / "validation" / "validation_2d_22.csv"
    v3m_csv = out_dir / "validation" / "validation_3d_31_modes.csv"
    v3s_csv = out_dir / "validation" / "validation_3d_31_summary.csv"
    bench_detail_csv = out_dir / "benchmark" / "backend_benchmark_detail.csv"
    bench_summary_csv = out_dir / "benchmark" / "backend_benchmark_summary.csv"
    log_file = out_dir / "run_all.log"

    mode_rows = _read_csv(mode_csv)
    v2d_rows = _read_csv(v2d_csv)
    v3m_rows = _read_csv(v3m_csv)
    v3s_rows = _read_csv(v3s_csv)
    bench_detail_rows = _read_csv(bench_detail_csv)
    bench_summary_rows = _read_csv(bench_summary_csv)
    log_info = _analyze_log(log_file)

    _write_report(
        report_path=report,
        out_dir=out_dir,
        mode_rows=mode_rows,
        v2d_rows=v2d_rows,
        v3m_rows=v3m_rows,
        v3s_rows=v3s_rows,
        bench_detail_rows=bench_detail_rows,
        bench_summary_rows=bench_summary_rows,
        log_info=log_info,
    )


if __name__ == "__main__":
    main()
