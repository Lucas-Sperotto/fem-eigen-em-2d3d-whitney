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
    elif formulation == "edge":
        case = "2.2.1_edge"
        fname = f"edge_{geometry}_{pol}_m{m}_{idx2_name}{idx2}_rank{rank:02d}_{field}"
    else:
        return None, None

    vtk = out_dir / "2d" / case / geometry / f"{fname}.vtk"
    png = out_dir / "img_all" / case / geometry / f"{fname}.png"
    return (vtk if vtk.exists() else None, png if png.exists() else None)


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
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        groups[(r.get("section", ""), r.get("case", ""))].append(r)

    out = []
    for (sec, case), gr in sorted(groups.items()):
        e1 = [abs(_to_float(r.get("err_primary_pct", "")) or 0.0) for r in gr if _to_float(r.get("err_primary_pct", "")) is not None]
        e2 = [abs(_to_float(r.get("err_secondary_pct", "")) or 0.0) for r in gr if _to_float(r.get("err_secondary_pct", "")) is not None]
        out.append(
            {
                "section": sec,
                "case": case,
                "rows": len(gr),
                "max_err_primary_pct": max(e1) if e1 else None,
                "max_err_secondary_pct": max(e2) if e2 else None,
            }
        )
    return out


def _max_3d_modes(mode_rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for r in mode_rows:
        groups[(r.get("solver", ""), r.get("case", ""))].append(r)

    out = []
    for (solver, case), gr in sorted(groups.items()):
        e_ana = [abs(_to_float(r.get("err_ana_pct", "")) or 0.0) for r in gr]
        e_ref = [abs(_to_float(r.get("err_ref_pct", "")) or 0.0) for r in gr]
        out.append(
            {
                "solver": solver,
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


def _write_report(
    report_path: Path,
    out_dir: Path,
    mode_rows: List[Dict[str, str]],
    v2d_rows: List[Dict[str, str]],
    v3m_rows: List[Dict[str, str]],
    v3s_rows: List[Dict[str, str]],
    log_info: Dict[str, object],
) -> None:
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    vtk_counts = _count_by_case_dir(out_dir / "2d", ".vtk")
    png_counts = _count_by_case_dir(out_dir / "img_all", ".png")

    mode_cov = Counter((r["geometry"], r["formulation"], r["polarization"]) for r in mode_rows)
    top_mode = _top_mode_errors(mode_rows, top_n=12)
    c2 = _group_2d22_case_stats(v2d_rows)
    c3 = _max_3d_modes(v3m_rows)

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
    lines.append(f"- CSV validacao 2D (2.2.x): {_md_link(out_dir / 'validation' / 'validation_2d_22.csv')}")
    lines.append(f"- CSV validacao 3D (modes): {_md_link(out_dir / 'validation' / 'validation_3d_31_modes.csv')}")
    lines.append(f"- CSV validacao 3D (summary): {_md_link(out_dir / 'validation' / 'validation_3d_31_summary.csv')}")
    lines.append(f"- Imagens de campos 2D: {_md_link(out_dir / 'img_all')}")
    lines.append(f"- Imagens de validacao 2.2.x: {_md_link(v22_img)}")
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
    lines.append("| Secao | Caso | Linhas | Max err primary (%) | Max err secondary (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for r in c2:
        lines.append(
            f"| `{r['section']}` | `{r['case']}` | {r['rows']} | {_fmt(r['max_err_primary_pct'])} | {_fmt(r['max_err_secondary_pct'])} |"
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
    lines.append("| Solver | Caso | Linhas | Max err ana (%) | Max err ref (%) |")
    lines.append("|---|---|---:|---:|---:|")
    for r in c3:
        lines.append(
            f"| `{r['solver']}` | `{r['case']}` | {r['rows']} | {_fmt(r['max_err_ana_pct'])} | {_fmt(r['max_err_ref_pct'])} |"
        )
    lines.append("")

    lines.append("### Summary detalhado (validation_3d_31_summary.csv)")
    lines.append("")
    lines.append("| Solver | Caso | nx | ny | nz | n_modes | max_err_ana | mean_err_ana | max_err_ref | mean_err_ref |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in v3s_rows:
        lines.append(
            f"| `{r.get('solver','')}` | `{r.get('case','')}` | {r.get('nx','')} | {r.get('ny','')} | {r.get('nz','')} | "
            f"{r.get('n_modes','')} | {_fmt(_to_float(r.get('max_err_ana_pct','')))} | {_fmt(_to_float(r.get('mean_err_ana_pct','')))} | "
            f"{_fmt(_to_float(r.get('max_err_ref_pct','')))} | {_fmt(_to_float(r.get('mean_err_ref_pct','')))} |"
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
    log_file = out_dir / "run_all.log"

    mode_rows = _read_csv(mode_csv)
    v2d_rows = _read_csv(v2d_csv)
    v3m_rows = _read_csv(v3m_csv)
    v3s_rows = _read_csv(v3s_csv)
    log_info = _analyze_log(log_file)

    _write_report(
        report_path=report,
        out_dir=out_dir,
        mode_rows=mode_rows,
        v2d_rows=v2d_rows,
        v3m_rows=v3m_rows,
        v3s_rows=v3s_rows,
        log_info=log_info,
    )


if __name__ == "__main__":
    main()
