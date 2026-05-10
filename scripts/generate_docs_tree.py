#!/usr/bin/env python3
"""
Gera a arvore de documentacao em docs/results/:

  docs/results/README.md             <- indice dos 14 casos
  docs/results/caso_NN_*.md          <- pagina por caso (14 arquivos)
  docs/results/fem_vs_efgmi.md       <- comparativo consolidado FEM x EFGMI

Cada pagina inclui:
  - descricao do problema e link para teoria
  - figura do artigo
  - condicoes de calculo (CLI, malha, backend)
  - tabela de resultados FEM (extraida do modes.csv)
  - imagens de campos FEM
  - tabela e imagens EFGMI (quando disponivel)
  - comparacao de timing FEM x EFGMI

Uso:
  python3 scripts/generate_docs_tree.py
  python3 scripts/generate_docs_tree.py --dry-run   # lista arquivos sem criar
"""

from __future__ import annotations

import argparse
import csv
import json
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DOCS_RESULTS = ROOT / "docs" / "results"
IMG_DIR = DOCS_RESULTS / "img"
MANIFEST_PATH = IMG_DIR / "manifest.json"
OUT = ROOT / "out"
EFGMI_BASE = OUT / "fem_efgmi_mode_report_base" / "base"


# ---------------------------------------------------------------------------
# Utilitarios
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _fmt(val: Any, decimals: int = 4) -> str:
    try:
        return f"{float(val):.{decimals}f}"
    except (ValueError, TypeError):
        return str(val)


def _pct(val: Any) -> str:
    try:
        return f"{float(val):.3f}"
    except (ValueError, TypeError):
        return str(val)


def _img_links(paths: list[str], prefix: str = "img/", max_imgs: int | None = None,
               cols: int = 3) -> str:
    """Gera grade de imagens em markdown."""
    if not paths:
        return "_Sem imagens disponiveis._\n"
    shown = paths if max_imgs is None else paths[:max_imgs]
    lines: list[str] = []
    for i, p in enumerate(shown):
        alt = Path(p).stem
        lines.append(f"![{alt}]({prefix}{p})")
        if (i + 1) % cols == 0:
            lines.append("")
    return "\n".join(lines) + "\n"


def _pick_summary_imgs(imgs: list[str]) -> list[str]:
    """Seleciona imagens de resumo (error_by_mode, rho, table summaries)."""
    return [p for p in imgs if any(kw in Path(p).name for kw in
                                   ["error_by_mode", "rho_by_mode", "table8", "table9",
                                    "table10", "preview", "k0L", "beta_over", "error_pct",
                                    "branches", "candidates"])]


def _pick_mode_imgs(imgs: list[str], max_imgs: int = 9) -> list[str]:
    """Seleciona imagens de campos modais (exclui resumos)."""
    mode_imgs = [p for p in imgs if not any(kw in Path(p).name for kw in
                                             ["error_by_mode", "rho_by_mode", "table8", "table9",
                                              "table10", "preview", "k0L", "beta_over",
                                              "error_pct", "branches", "candidates"])]
    return mode_imgs[:max_imgs]


def _load_manifest() -> dict[str, dict[str, list[str]]]:
    if not MANIFEST_PATH.exists():
        return {}
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _read_timing_generic(path: Path) -> dict[str, str]:
    rows = _read_csv(path)
    if not rows:
        return {}
    r = rows[0]
    return {k: r.get(k, "") for k in ["assembly_ms", "solve_ms", "post_ms", "total_ms",
                                       "assembly_ms_total", "solve_ms_total",
                                       "mesh_nodes", "mesh_tris", "mesh_tets", "nx", "ny", "nr", "nt"]}


def _timing_str(t: dict[str, str]) -> str:
    asm = t.get("assembly_ms") or t.get("assembly_ms_total", "?")
    slv = t.get("solve_ms") or t.get("solve_ms_total", "?")
    pst = t.get("post_ms", "?")
    tot = t.get("total_ms", "?")
    try:
        return (f"assembly={float(asm):.1f} ms | solve={float(slv):.1f} ms | "
                f"post={float(pst):.1f} ms | total={float(tot):.1f} ms")
    except (ValueError, TypeError):
        return f"assembly={asm} | solve={slv} | post={pst} | total={tot}"


# ---------------------------------------------------------------------------
# Leitores de CSV por familia
# ---------------------------------------------------------------------------

def _read_helm10_modes(path: Path, max_rows: int = 10) -> str:
    rows = _read_csv(path)[:max_rows]
    if not rows:
        return "_Sem dados._\n"
    header = "| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |"
    sep    = "|---|---|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r.get('mode_label','?')} | {r.get('family','?')} "
            f"| {_fmt(r.get('kc_ana','?'))} | {_fmt(r.get('kc_fem','?'))} "
            f"| {_pct(r.get('error_percent','?'))} | {_fmt(r.get('rho_abs','?'),6)} |"
        )
    return "\n".join(lines) + "\n"


def _read_helmvec_modes(path: Path, max_rows: int = 10) -> str:
    rows = _read_csv(path)[:max_rows]
    if not rows:
        return "_Sem dados._\n"
    header = "| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |"
    sep    = "|---|---|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r.get('mode_label','?')} | {r.get('family','?')} "
            f"| {_fmt(r.get('kc_ana','?'))} | {_fmt(r.get('kc_fem','?'))} "
            f"| {_pct(r.get('error_percent','?'))} | {_fmt(r.get('rho_abs','?'),6)} |"
        )
    return "\n".join(lines) + "\n"


def _read_helmvec1_modes(path: Path, max_rows: int = 12) -> str:
    rows = _read_csv(path)[:max_rows]
    if not rows:
        return "_Sem dados._\n"
    header = "| Form. | Bloco dom. | Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |"
    sep    = "|---|---|---|---|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r.get('formulation','?')} | {r.get('dominant_block','?')} "
            f"| {r.get('mode_label','?')} | {r.get('family','?')} "
            f"| {_fmt(r.get('kc_ana','?'))} | {_fmt(r.get('kc_fem','?'))} "
            f"| {_pct(r.get('error_percent','?'))} | {_fmt(r.get('rho_abs','?'),6)} |"
        )
    return "\n".join(lines) + "\n"


def _read_helmvec2_modes(path: Path) -> str:
    rows = _read_csv(path)
    if not rows:
        return "_Sem dados._\n"
    header = "| Modo | `k0L` FEM | Ref. HELMVEC2 | Ref. Hayata | Erro HELMVEC2 (%) | Erro Hayata (%) |"
    sep    = "|---|---:|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r.get('mode','?')} "
            f"| {_fmt(r.get('k0l_fem_matched','?'))} "
            f"| {_fmt(r.get('ref_helmvec2','?'))} "
            f"| {_fmt(r.get('ref_hayata','?'))} "
            f"| {_pct(r.get('error_percent_helmvec2','?'))} "
            f"| {_pct(r.get('error_percent_hayata','?'))} |"
        )
    return "\n".join(lines) + "\n"


def _read_helmvec3_fig12(path: Path) -> str:
    rows = _read_csv(path)
    if not rows:
        return "_Sem dados._\n"
    header = "| `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo | Erro analitico (%) | Erro artigo (%) |"
    sep    = "|---|---:|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {_fmt(r.get('br_over_lambda0','?'),3)} "
            f"| {_fmt(r.get('beta_over_k0_fem','?'))} "
            f"| {_fmt(r.get('beta_over_k0_analytic','?'),4)} "
            f"| {_fmt(r.get('beta_over_k0_helmvec3','?'),4)} "
            f"| {_pct(r.get('error_percent_analytic','?'))} "
            f"| {_pct(r.get('error_percent_helmvec3','?'))} |"
        )
    return "\n".join(lines) + "\n"


def _read_helmvec3_fig13_table10(path: Path) -> str:
    rows = _read_csv(path)
    if not rows:
        return "_Sem dados._\n"
    header = "| `d/a` | `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo/HELMVEC3 | Erro analítico (%) | Erro artigo/HELMVEC3 (%) |"
    sep    = "|---|---|---:|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {_fmt(r.get('d_over_a','?'),3)} "
            f"| {_fmt(r.get('br_over_lambda0','?'),3)} "
            f"| {_fmt(r.get('beta_over_k0_fem_matched','?'))} "
            f"| {_fmt(r.get('beta_over_k0_analytic','?'),4)} "
            f"| {_fmt(r.get('beta_over_k0_helmvec3','?'),4)} "
            f"| {_pct(r.get('error_percent_analytic','?'))} "
            f"| {_pct(r.get('error_percent_helmvec3','?'))} |"
        )
    return "\n".join(lines) + "\n"


def _read_fem3d_modes(path: Path, max_rows: int = 10) -> str:
    rows = _read_csv(path)[:max_rows]
    if not rows:
        return "_Sem dados._\n"
    header = "| # | Modo | `k0` analitico | `k0` FEM | Ref. artigo | Erro analítico (%) | Erro artigo (%) |"
    sep    = "|---|---|---:|---:|---:|---:|---:|"
    lines = [header, sep]
    for r in rows:
        lines.append(
            f"| {r.get('reference_index','?')} | {r.get('mode_label','?')} "
            f"| {_fmt(r.get('k0_analytic','?'))} "
            f"| {_fmt(r.get('k0_fem','?'))} "
            f"| {_fmt(r.get('ref_paper','?'))} "
            f"| {_pct(r.get('error_percent_analytic','?'))} "
            f"| {_pct(r.get('error_percent_ref_paper','?'))} |"
        )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Especificacoes dos 14 casos
# ---------------------------------------------------------------------------

@dataclass
class CaseSpec:
    num: int
    slug: str
    title: str
    section: str
    formula: str
    quantity: str
    geometry_desc: str
    theory_link: str
    article_figs: list[str]
    cli_cmd: str
    mesh_desc: str
    fem_modes_path: Path
    fem_timing_path: Path
    efgm_modes_path: Path | None
    efgm_timing_path: Path | None
    manifest_key: str
    modes_reader: str  # "helm10", "helmvec", "helmvec1", "helmvec2", "hv3fig12", "hv3fig13", "fem3d"
    extra_notes: str = ""


CASES: list[CaseSpec] = [
    CaseSpec(
        num=1, slug="tab1_helm10_rect",
        title="Tabela 1 — Guia Retangular Escalar (Sec. 2.1)",
        section="2.1", formula="Eq. (43) — S φ = kc² T φ", quantity="`kc` (cutoff escalar)",
        geometry_desc="Guia retangular homogêneo, a/b = 2, PEC",
        theory_link="../traducao/2.1_Guias de Onda Homogêneos.md",
        article_figs=["figura4.png", "figura5.png" if (ROOT / "docs/figs/figura5.png").exists() else None],
        cli_cmd="./build/helm10_rect --ar-m 1.0 --nx 10 --ny 20 --nmodos 10 --backend closed-form",
        mesh_desc="231 nós, 400 triângulos (nx=10, ny=20)",
        fem_modes_path=OUT / "helm10/rect/csv/helm10_rect_modes.csv",
        fem_timing_path=OUT / "helm10/rect/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helm10/rect/csv/helm10_rect_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helm10/rect/run_timing.csv",
        manifest_key="helm10_rect", modes_reader="helm10",
    ),
    CaseSpec(
        num=2, slug="tab2_helm10_circle",
        title="Tabela 2 — Guia Circular Escalar (Sec. 2.1)",
        section="2.1", formula="Eq. (43)", quantity="`kc`",
        geometry_desc="Guia circular homogêneo, raio unitário, PEC",
        theory_link="../traducao/2.1_Guias de Onda Homogêneos.md",
        article_figs=["figura6.png"],
        cli_cmd="./build/helm10_circle --nr 8 --nt 15 --nmodos 10 --backend closed-form",
        mesh_desc="121 nós, 225 triângulos (nr=8, nt=15)",
        fem_modes_path=OUT / "helm10/circle/csv/helm10_circle_modes.csv",
        fem_timing_path=OUT / "helm10/circle/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helm10/circle/csv/helm10_circle_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helm10/circle/run_timing.csv",
        manifest_key="helm10_circle", modes_reader="helm10",
    ),
    CaseSpec(
        num=3, slug="tab3_helm10_coax",
        title="Tabela 3 — Linha Coaxial Escalar (Sec. 2.1)",
        section="2.1", formula="Eq. (43)", quantity="`kc`",
        geometry_desc="Linha coaxial homogênea, r₂/r₁ = 4",
        theory_link="../traducao/2.1_Guias de Onda Homogêneos.md",
        article_figs=["figura8.png"],
        cli_cmd="./build/helm10_coax --nr 10 --nt 17 --nmodos 10 --backend closed-form",
        mesh_desc="187 nós, 340 triângulos (nr=10, nt=17)",
        fem_modes_path=OUT / "helm10/coax/csv/helm10_coax_modes.csv",
        fem_timing_path=OUT / "helm10/coax/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helm10/coax/csv/helm10_coax_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helm10/coax/run_timing.csv",
        manifest_key="helm10_coax", modes_reader="helm10",
    ),
    CaseSpec(
        num=4, slug="tab4_helmvec_rect",
        title="Tabela 4 — Guia Retangular Vetorial (Sec. 2.2.1)",
        section="2.2.1", formula="Eq. (65) — S e = kc² T e (edge 2D)", quantity="`kc` vetorial",
        geometry_desc="Guia retangular homogêneo, elementos de aresta Whitney",
        theory_link="../traducao/2.2.1_Guia de Onda Campos Vetoriais.md",
        article_figs=["figura4.png"],
        cli_cmd="./build/edge_rect --nx 10 --ny 20 --nmodos 10 --backend closed-form",
        mesh_desc="231 nós, 400 triângulos (nx=10, ny=20)",
        fem_modes_path=OUT / "helmvec/rect/csv/edge_rect_modes.csv",
        fem_timing_path=OUT / "helmvec/rect/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="helmvec_rect", modes_reader="helmvec",
    ),
    CaseSpec(
        num=5, slug="tab5_helmvec_circle",
        title="Tabela 5 — Guia Circular Vetorial (Sec. 2.2.1)",
        section="2.2.1", formula="Eq. (65)", quantity="`kc` vetorial",
        geometry_desc="Guia circular homogêneo, raio unitário, elementos de aresta",
        theory_link="../traducao/2.2.1_Guia de Onda Campos Vetoriais.md",
        article_figs=["figura6.png"],
        cli_cmd="./build/edge_circle --nr 8 --nt 15 --nmodos 10 --backend closed-form",
        mesh_desc="121 nós, 225 triângulos (nr=8, nt=15)",
        fem_modes_path=OUT / "helmvec/circle/csv/edge_circle_modes.csv",
        fem_timing_path=OUT / "helmvec/circle/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="helmvec_circle", modes_reader="helmvec",
    ),
    CaseSpec(
        num=6, slug="tab6_helmvec1_rect",
        title="Tabela 6 — Guia Retangular Misto 3 Comp. (Sec. 2.2.2)",
        section="2.2.2", formula="Eq. (92) — bloco diagonal [St 0; 0 Sz] e = kc² [Tt 0; 0 Tz] e",
        quantity="`kc` sistema misto",
        geometry_desc="Guia retangular homogêneo, formulações E e H",
        theory_link="../traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md",
        article_figs=["figura4.png"],
        cli_cmd="./build/mixed_rect --nx 10 --ny 20 --backend closed-form",
        mesh_desc="231 nós, 400 triângulos (nx=10, ny=20)",
        fem_modes_path=OUT / "helmvec1/rect/csv/mixed_rect_modes.csv",
        fem_timing_path=OUT / "helmvec1/rect/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helmvec1/rect/csv/mixed_rect_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helmvec1/rect/run_timing.csv",
        manifest_key="helmvec1_rect", modes_reader="helmvec1",
    ),
    CaseSpec(
        num=7, slug="tab7_helmvec1_circle",
        title="Tabela 7 — Guia Circular Misto 3 Comp. (Sec. 2.2.2)",
        section="2.2.2", formula="Eq. (92)", quantity="`kc` sistema misto",
        geometry_desc="Guia circular homogêneo, raio unitário, formulações E e H",
        theory_link="../traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md",
        article_figs=["figura6.png"],
        cli_cmd="./build/mixed_circle --nr 8 --nt 15 --backend closed-form",
        mesh_desc="121 nós, 225 triângulos (nr=8, nt=15)",
        fem_modes_path=OUT / "helmvec1/circle/csv/mixed_circle_modes.csv",
        fem_timing_path=OUT / "helmvec1/circle/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helmvec1/circle/csv/mixed_circle_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helmvec1/circle/run_timing.csv",
        manifest_key="helmvec1_circle", modes_reader="helmvec1",
    ),
    CaseSpec(
        num=8, slug="fig11_tab8_helmvec2",
        title="Figura 11 / Tabela 8 — Guia Parcialmente Preenchido, `k0` dado `β` (Sec. 2.2.3)",
        section="2.2.3", formula="Eq. (119) — A x = k0² B x, x = [Et; Ez]",
        quantity="`k0L` (número de onda × comprimento)",
        geometry_desc="Guia retangular quadrado, dielétrico inferior (εr=1.5), β=10",
        theory_link="../traducao/2.2.3_Determinação do número de onda.md",
        article_figs=["figura11.png"],
        cli_cmd="./build/helmvec2_rect --beta 10 --nx 20 --ny 20 --backend closed-form",
        mesh_desc="441 nós, 800 triângulos (nx=20, ny=20)",
        fem_modes_path=OUT / "helmvec2/rect/csv/helmvec2_rect_modes.csv",
        fem_timing_path=OUT / "helmvec2/rect/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helmvec2/rect/csv/helmvec2_rect_modes.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helmvec2/rect/run_timing.csv",
        manifest_key="helmvec2_rect", modes_reader="helmvec2",
        extra_notes="> **Nota:** A Eq. (120) impressa no artigo contém inconsistência de impressão "
                    "(falta fator β² no bloco de massa vetorial). "
                    "O código reconstrói A_tt a partir dos blocos elementares validados. "
                    "Ver [src/helmvec2/README.md](../../src/helmvec2/README.md).",
    ),
    CaseSpec(
        num=9, slug="fig12_tab9_helmvec3",
        title="Figura 12 / Tabela 9 — Dispersão, `β` dado `k0`, Exemplo 1 (Sec. 2.2.4)",
        section="2.2.4", formula="Eq. (136) — P x = β² Q x",
        quantity="`β/k0` (relação de dispersão)",
        geometry_desc="Guia retangular, b/a=0.45, d/a=0.5, εr=2.45",
        theory_link="../traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md",
        article_figs=["figura12.png"],
        cli_cmd="./build/helmvec3_fig12_rect --nx 10 --ny 5 --backend closed-form",
        mesh_desc="66 nós, 100 triângulos (nx=10, ny=5)",
        fem_modes_path=OUT / "helmvec3/fig12_rect/csv/helmvec3_fig12_rect_table9.csv",
        fem_timing_path=OUT / "helmvec3/fig12_rect/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_table9.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helmvec3/fig12_rect/run_timing.csv",
        manifest_key="helmvec3_fig12", modes_reader="hv3fig12",
    ),
    CaseSpec(
        num=10, slug="fig13_tab10_helmvec3",
        title="Figura 13 / Tabela 10 — Dispersão, `β` dado `k0`, Exemplo 2 (Sec. 2.2.4)",
        section="2.2.4", formula="Eq. (136)", quantity="`β/k0` — múltiplos ramos por d/a",
        geometry_desc="Guia retangular, b/a=0.45, εr=2.45, d/a variável",
        theory_link="../traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md",
        article_figs=["figura13.png"],
        cli_cmd="./build/helmvec3_fig13_rect --d-over-a-preview 0.20 --nx 10 --ny 5 --backend closed-form",
        mesh_desc="66 nós, 100 triângulos (nx=10, ny=5)",
        fem_modes_path=OUT / "helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv",
        fem_timing_path=OUT / "helmvec3/fig13_rect/run_timing.csv",
        efgm_modes_path=EFGMI_BASE / "efgm/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv",
        efgm_timing_path=EFGMI_BASE / "efgm/helmvec3/fig13_rect/run_timing.csv",
        manifest_key="helmvec3_fig13", modes_reader="hv3fig13",
    ),
    CaseSpec(
        num=11, slug="tab12_fem3d_air",
        title="Tabela 12 — Cavidade Retangular 3D, Ar (Sec. 3.1)",
        section="3.1", formula="Eq. (178) — S e = k0² T e (aresta tetraédrica)",
        quantity="`k0` (número de onda de ressonância)",
        geometry_desc="Cavidade retangular preenchida com ar, PEC",
        theory_link="../traducao/3_Problemas_Tridimensionais.md",
        article_figs=["figura15.png"],
        cli_cmd="./build/fem3d0_air --backend closed-form\n./build/fem3d1_air --backend closed-form",
        mesh_desc="~112 nós, ~324 tetraedros (FEM3D0 denso / FEM3D1 esparso)",
        fem_modes_path=OUT / "fem3d0/air/csv/fem3d0_air_modes.csv",
        fem_timing_path=OUT / "fem3d0/air/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="fem3d0_air", modes_reader="fem3d",
    ),
    CaseSpec(
        num=12, slug="tab13_fem3d_half",
        title="Tabela 13 — Cavidade Retangular 3D Semi-Preenchida (Sec. 3.1)",
        section="3.1", formula="Eq. (178)", quantity="`k0`",
        geometry_desc="Cavidade retangular, dielétrico εr=2 em z=[0.5, 1] cm",
        theory_link="../traducao/3_Problemas_Tridimensionais.md",
        article_figs=["figura16.png"],
        cli_cmd="./build/fem3d0_half --backend closed-form\n./build/fem3d1_half --backend closed-form",
        mesh_desc="~200 nós, ~615 tetraedros",
        fem_modes_path=OUT / "fem3d0/half/csv/fem3d0_half_modes.csv",
        fem_timing_path=OUT / "fem3d0/half/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="fem3d0_half", modes_reader="fem3d",
    ),
    CaseSpec(
        num=13, slug="tab14_fem3d_cyl",
        title="Tabela 14 — Cavidade Cilíndrica 3D, Ar (Sec. 3.1)",
        section="3.1", formula="Eq. (178)", quantity="`k0`",
        geometry_desc="Cavidade cilíndrica circular com ar",
        theory_link="../traducao/3_Problemas_Tridimensionais.md",
        article_figs=["figura17.png"],
        cli_cmd="./build/fem3d0_cyl --backend closed-form\n./build/fem3d1_cyl --backend closed-form",
        mesh_desc="~200 nós, ~633 tetraedros",
        fem_modes_path=OUT / "fem3d0/cyl/csv/fem3d0_cyl_modes.csv",
        fem_timing_path=OUT / "fem3d0/cyl/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="fem3d0_cyl", modes_reader="fem3d",
    ),
    CaseSpec(
        num=14, slug="tab15_fem3d_sphere",
        title="Tabela 15 — Cavidade Esférica 3D (Sec. 3.1)",
        section="3.1", formula="Eq. (178)", quantity="`k0`",
        geometry_desc="Cavidade esférica, raio 1 cm",
        theory_link="../traducao/3_Problemas_Tridimensionais.md",
        article_figs=[],
        cli_cmd="./build/fem3d0_sphere --backend closed-form\n./build/fem3d1_sphere --backend closed-form",
        mesh_desc="~166 nós, ~473 tetraedros",
        fem_modes_path=OUT / "fem3d0/sphere/csv/fem3d0_sphere_modes.csv",
        fem_timing_path=OUT / "fem3d0/sphere/run_timing.csv",
        efgm_modes_path=None, efgm_timing_path=None,
        manifest_key="fem3d0_sphere", modes_reader="fem3d",
    ),
]


# ---------------------------------------------------------------------------
# Geradores de .md
# ---------------------------------------------------------------------------

def _modes_table(spec: CaseSpec, efgm: bool = False) -> str:
    path = spec.efgm_modes_path if efgm else spec.fem_modes_path
    if path is None or not path.exists():
        return "_Sem dados disponíveis._\n"
    readers = {
        "helm10": _read_helm10_modes,
        "helmvec": _read_helmvec_modes,
        "helmvec1": _read_helmvec1_modes,
        "helmvec2": _read_helmvec2_modes,
        "hv3fig12": _read_helmvec3_fig12,
        "hv3fig13": _read_helmvec3_fig13_table10,
        "fem3d": _read_fem3d_modes,
    }
    return readers[spec.modes_reader](path)


def _timing_row(label: str, path: Path | None) -> str:
    if path is None or not path.exists():
        return f"| {label} | — | — | — | — |\n"
    t = _read_timing_generic(path)
    asm = t.get("assembly_ms") or t.get("assembly_ms_total", "?")
    slv = t.get("solve_ms") or t.get("solve_ms_total", "?")
    pst = t.get("post_ms", "?")
    tot = t.get("total_ms", "?")
    try:
        return (f"| {label} | {float(asm):.1f} | {float(slv):.1f} "
                f"| {float(pst):.1f} | {float(tot):.1f} |\n")
    except (ValueError, TypeError):
        return f"| {label} | {asm} | {slv} | {pst} | {tot} |\n"


def generate_case_md(spec: CaseSpec, manifest: dict, dry_run: bool = False) -> Path:
    out_path = DOCS_RESULTS / f"caso_{spec.num:02d}_{spec.slug}.md"
    if dry_run:
        print(f"  [dry-run] {out_path.relative_to(ROOT)}")
        return out_path

    case_imgs = manifest.get(spec.manifest_key, {"fem": [], "efgm": []})
    fem_imgs = case_imgs.get("fem", [])
    efgm_imgs = case_imgs.get("efgm", [])

    fem_summary = _pick_summary_imgs(fem_imgs)
    fem_modes_imgs = _pick_mode_imgs(fem_imgs, max_imgs=9)
    efgm_summary = _pick_summary_imgs(efgm_imgs)
    efgm_modes_imgs = _pick_mode_imgs(efgm_imgs, max_imgs=9)

    # Artigo figures
    artigo_block = ""
    for fig in spec.article_figs:
        if fig and (ROOT / "docs" / "figs" / fig).exists():
            label = fig.replace(".png", "").replace("figura", "Figura ")
            artigo_block += f"![{label}](../figs/{fig})\n"
    if not artigo_block:
        artigo_block = "_Figura do artigo não disponível._\n"

    lines: list[str] = [
        f"# Caso {spec.num:02d} — {spec.title}",
        "",
        "## Problema",
        "",
        f"- **Seção do artigo:** {spec.section}",
        f"- **Formulação:** {spec.formula}",
        f"- **Grandeza calculada:** {spec.quantity}",
        f"- **Geometria:** {spec.geometry_desc}",
        f"- **Teoria:** [{spec.theory_link.split('/')[-1]}]({spec.theory_link})",
        f"- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)",
        "",
        "## Figura do artigo",
        "",
        artigo_block,
        "## Condições de cálculo",
        "",
        "| Parâmetro | Valor |",
        "|---|---|",
        f"| Malha | {spec.mesh_desc} |",
        "| Backend | `closed-form` |",
        "",
        "```bash",
        spec.cli_cmd,
        "```",
        "",
    ]

    if spec.extra_notes:
        lines += [spec.extra_notes, ""]

    # FEM results
    lines += [
        "## Resultados — FEM",
        "",
        _modes_table(spec, efgm=False),
    ]
    if fem_summary:
        lines += ["### Gráficos de resumo — FEM", "", _img_links(fem_summary, cols=2)]
    if fem_modes_imgs:
        lines += ["### Campos modais — FEM", "", _img_links(fem_modes_imgs, cols=3)]

    # Timing
    lines += [
        "## Tempo de execução",
        "",
        "| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |",
        "|---|---:|---:|---:|---:|",
    ]
    lines.append(_timing_row("FEM", spec.fem_timing_path))
    if spec.efgm_timing_path:
        lines.append(_timing_row("EFGMI", spec.efgm_timing_path))

    # EFGMI results
    if spec.efgm_modes_path is not None:
        lines += [
            "",
            "## Resultados — EFGMI",
            "",
            _modes_table(spec, efgm=True),
        ]
        if efgm_summary:
            lines += ["### Gráficos de resumo — EFGMI", "", _img_links(efgm_summary, cols=2)]
        if efgm_modes_imgs:
            lines += ["### Campos modais — EFGMI", "", _img_links(efgm_modes_imgs, cols=3)]

    # Navigation
    prev_link = f"[← Caso {spec.num-1:02d}](caso_{spec.num-1:02d}_{''.join(c.slug for c in CASES if c.num == spec.num-1)}.md)" if spec.num > 1 else ""
    next_link = f"[Caso {spec.num+1:02d} →](caso_{spec.num+1:02d}_{''.join(c.slug for c in CASES if c.num == spec.num+1)}.md)" if spec.num < 14 else ""

    lines += [
        "",
        "---",
        "",
        f"{prev_link}  {next_link}".strip() + "  [↑ Índice de Resultados](README.md)",
        "",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Gerado: {out_path.relative_to(ROOT)}")
    return out_path


def generate_results_readme(manifest: dict, dry_run: bool = False) -> Path:
    out_path = DOCS_RESULTS / "README.md"
    if dry_run:
        print(f"  [dry-run] {out_path.relative_to(ROOT)}")
        return out_path

    lines: list[str] = [
        "# Resultados Numéricos — NASA TP-3485",
        "",
        "Reprodução computacional de todos os 14 casos numéricos publicados no artigo.",
        "Cada página cobre: problema, condições de cálculo, resultados FEM, campos modais e comparação com EFGMI.",
        "",
        "→ [Comparação consolidada FEM × EFGMI](fem_vs_efgmi.md)",
        "→ [Índice geral de documentação](../INDICE.md)",
        "",
        "## Casos 2D — Escalar (Sec. 2.1)",
        "",
        "| # | Caso | Geometria | Grandeza | Figura |",
        "|---|---|---|---|---|",
    ]
    for c in CASES[:3]:
        fig = f"![](../figs/{c.article_figs[0]})" if c.article_figs else "—"
        lines.append(f"| {c.num:02d} | [{c.title}](caso_{c.num:02d}_{c.slug}.md) "
                     f"| {c.geometry_desc} | {c.quantity} | {fig} |")

    lines += [
        "",
        "## Casos 2D — Vetorial Edge (Sec. 2.2.1)",
        "",
        "| # | Caso | Geometria | Grandeza | Figura |",
        "|---|---|---|---|---|",
    ]
    for c in CASES[3:5]:
        fig = f"![](../figs/{c.article_figs[0]})" if c.article_figs else "—"
        lines.append(f"| {c.num:02d} | [{c.title}](caso_{c.num:02d}_{c.slug}.md) "
                     f"| {c.geometry_desc} | {c.quantity} | {fig} |")

    lines += [
        "",
        "## Casos 2D — Sistema Misto (Sec. 2.2.2)",
        "",
        "| # | Caso | Geometria | Grandeza | Figura |",
        "|---|---|---|---|---|",
    ]
    for c in CASES[5:7]:
        fig = f"![](../figs/{c.article_figs[0]})" if c.article_figs else "—"
        lines.append(f"| {c.num:02d} | [{c.title}](caso_{c.num:02d}_{c.slug}.md) "
                     f"| {c.geometry_desc} | {c.quantity} | {fig} |")

    lines += [
        "",
        "## Casos 2D — Acoplados (Sec. 2.2.3 e 2.2.4)",
        "",
        "| # | Caso | Geometria | Grandeza | Figura |",
        "|---|---|---|---|---|",
    ]
    for c in CASES[7:10]:
        fig = f"![](../figs/{c.article_figs[0]})" if c.article_figs else "—"
        lines.append(f"| {c.num:02d} | [{c.title}](caso_{c.num:02d}_{c.slug}.md) "
                     f"| {c.geometry_desc} | {c.quantity} | {fig} |")

    lines += [
        "",
        "## Casos 3D — Cavidades Ressonantes (Sec. 3.1)",
        "",
        "| # | Caso | Geometria | Grandeza | Figura |",
        "|---|---|---|---|---|",
    ]
    for c in CASES[10:]:
        fig = f"![](../figs/{c.article_figs[0]})" if c.article_figs else "—"
        lines.append(f"| {c.num:02d} | [{c.title}](caso_{c.num:02d}_{c.slug}.md) "
                     f"| {c.geometry_desc} | {c.quantity} | {fig} |")

    lines += ["", "---", "", "→ [README do repositório](../../README.md)"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Gerado: {out_path.relative_to(ROOT)}")
    return out_path


def generate_fem_vs_efgmi_md(manifest: dict, dry_run: bool = False) -> Path:
    out_path = DOCS_RESULTS / "fem_vs_efgmi.md"
    if dry_run:
        print(f"  [dry-run] {out_path.relative_to(ROOT)}")
        return out_path

    lines: list[str] = [
        "# Comparação FEM × EFGMI",
        "",
        "Tabela consolidada de timing e erro médio por caso, comparando o backend FEM (formas fechadas) "
        "e o backend EFGMI (interpolantes de consistência, malha triangular como fundo de integração).",
        "",
        "Para detalhes por modo e campos individuais, veja as páginas de cada caso.",
        "",
        "## Timing por caso",
        "",
        "| Caso | FEM assembly (ms) | FEM solve (ms) | FEM total (ms) | EFGMI assembly (ms) | EFGMI solve (ms) | EFGMI total (ms) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for spec in CASES:
        if spec.efgm_timing_path is None:
            continue
        fem_t = _read_timing_generic(spec.fem_timing_path)
        efgm_t = _read_timing_generic(spec.efgm_timing_path)

        def _v(t: dict, key: str) -> str:
            val = t.get(key) or t.get(key + "_total", "?")
            try:
                return f"{float(val):.1f}"
            except (ValueError, TypeError):
                return "?"

        lines.append(
            f"| [{spec.num:02d}](caso_{spec.num:02d}_{spec.slug}.md) "
            f"| {_v(fem_t,'assembly_ms')} | {_v(fem_t,'solve_ms')} | {_v(fem_t,'total_ms')} "
            f"| {_v(efgm_t,'assembly_ms')} | {_v(efgm_t,'solve_ms')} | {_v(efgm_t,'total_ms')} |"
        )

    lines += [
        "",
        "> **Nota:** EFGM tem assembly mais custoso (construção de interpolantes nodais por "
        "consistência para cada elemento), mas solve comparável (mesmo LAPACK dsygv).",
        "",
        "## Imagens de resumo de erro — FEM",
        "",
    ]

    for spec in CASES[:7]:  # 2D cases com error_by_mode
        imgs = manifest.get(spec.manifest_key, {}).get("fem", [])
        summary = [p for p in imgs if "error_by_mode" in p]
        if summary:
            lines.append(f"### {spec.title}")
            lines.append("")
            lines += [f"![{Path(p).stem}](img/{p})" for p in summary[:2]]
            lines.append("")

    lines += [
        "## Imagens de resumo de erro — EFGMI",
        "",
    ]
    for spec in CASES[:7]:
        imgs = manifest.get(spec.manifest_key, {}).get("efgm", [])
        summary = [p for p in imgs if "error_by_mode" in p]
        if summary:
            lines.append(f"### {spec.title}")
            lines.append("")
            lines += [f"![{Path(p).stem}](img/{p})" for p in summary[:2]]
            lines.append("")

    lines += [
        "---",
        "",
        "→ [Índice de Resultados](README.md)  "
        "→ [Relatório completo FEM×EFGMI](../../out/fem_efgmi_mode_report_base/FEM_EFGM_MODE_REPORT.md)",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Gerado: {out_path.relative_to(ROOT)}")
    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Lista arquivos que seriam gerados sem criar nada.")
    args = parser.parse_args()

    manifest = _load_manifest()
    if not manifest:
        print("AVISO: manifest.json nao encontrado. Execute export_presentation_assets.py primeiro.")

    print("Gerando docs/results/ ...")

    generate_results_readme(manifest, dry_run=args.dry_run)
    generate_fem_vs_efgmi_md(manifest, dry_run=args.dry_run)
    for spec in CASES:
        generate_case_md(spec, manifest, dry_run=args.dry_run)

    if not args.dry_run:
        total = len(list(DOCS_RESULTS.glob("*.md")))
        print(f"\nPronto: {total} arquivos .md em docs/results/")


if __name__ == "__main__":
    main()
