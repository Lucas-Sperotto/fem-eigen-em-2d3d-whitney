#!/usr/bin/env python3
"""
Gera um relatorio Markdown comparando:
- FEM (`closed-form`)
- EFGM/hibrido (`efgmi`)

O foco e a validacao 2D com selecao explicita de malha:
- `base`: discretizacao do artigo;
- `doubled`: discretizacao dobrada.

Para cada caso, o relatorio registra:
- discretizacao usada;
- numero de nos, triangulos e arestas da malha;
- tempos reportados pelo executavel;
- tabelas por modo/ponto do artigo, comparando FEM e EFGM.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class ModeRow:
    group: str
    article_mode: str
    repo_mode_label: str
    family: str
    variable: str
    value: float
    ref_primary: float
    ref_secondary: float | None
    err_primary_pct: float
    err_secondary_pct: float | None
    rho_label: str
    rho_value: float


@dataclass(frozen=True)
class MeshPreset:
    key: str
    title: str


@dataclass(frozen=True)
class CaseSpec:
    key: str
    title: str
    section: str
    article_ref: str
    exe: str
    out_rel: Path
    base_args: tuple[str, ...]
    doubled_args: tuple[str, ...]
    primary_ref_label: str
    secondary_ref_label: str | None
    parse_modes: Callable[[Path], list[ModeRow]]


BASE = MeshPreset("base", "Discretizacao Base do Artigo")
DOUBLED = MeshPreset("doubled", "Discretizacao Dobrada")
PRESET_BY_KEY = {
    BASE.key: BASE,
    DOUBLED.key: DOUBLED,
}
BACKENDS = (
    ("FEM", "closed-form"),
    ("EFGM", "efgmi"),
)


@dataclass(frozen=True)
class ScalarRequiredMode:
    article_mode: str
    family: str
    mode_label: str


@dataclass(frozen=True)
class MixedRequiredMode:
    group: str
    article_mode: str
    formulation: str
    dominant_block: str
    family: str
    mode_label: str


SCALAR_RECT_MODES = (
    ScalarRequiredMode("TE10", "TE", "TE_m1_n0"),
    ScalarRequiredMode("TE20", "TE", "TE_m2_n0"),
    ScalarRequiredMode("TE01", "TE", "TE_m0_n1"),
    ScalarRequiredMode("TE11", "TE", "TE_m1_n1"),
    ScalarRequiredMode("TM11", "TM", "TM_m1_n1"),
    ScalarRequiredMode("TE12", "TE", "TE_m1_n2"),
    ScalarRequiredMode("TM12", "TM", "TM_m1_n2"),
    ScalarRequiredMode("TE21", "TE", "TE_m2_n1"),
    ScalarRequiredMode("TM21", "TM", "TM_m2_n1"),
)

SCALAR_CIRCLE_MODES = (
    ScalarRequiredMode("TE01", "TE", "TE_m0_p1"),
    ScalarRequiredMode("TE11", "TE", "TE_m1_p1"),
    ScalarRequiredMode("TE12", "TE", "TE_m1_p2"),
    ScalarRequiredMode("TE21", "TE", "TE_m2_p1"),
    ScalarRequiredMode("TE22", "TE", "TE_m2_p2"),
    ScalarRequiredMode("TM01", "TM", "TM_m0_p1"),
    ScalarRequiredMode("TM11", "TM", "TM_m1_p1"),
    ScalarRequiredMode("TM12", "TM", "TM_m1_p2"),
    ScalarRequiredMode("TM21", "TM", "TM_m2_p1"),
    ScalarRequiredMode("TM22", "TM", "TM_m2_p2"),
)

SCALAR_COAX_MODES = (
    ScalarRequiredMode("TE11", "TE", "TE_m1_p1"),
    ScalarRequiredMode("TE21", "TE", "TE_m2_p1"),
    ScalarRequiredMode("TE31", "TE", "TE_m3_p1"),
    ScalarRequiredMode("TM01", "TM", "TM_m0_p1"),
    ScalarRequiredMode("TM11", "TM", "TM_m1_p1"),
)

MIXED_RECT_MODES = (
    MixedRequiredMode("E / edge / TE", "TE10", "E", "edge", "TE", "TE_m1_n0"),
    MixedRequiredMode("E / edge / TE", "TE20", "E", "edge", "TE", "TE_m2_n0"),
    MixedRequiredMode("E / edge / TE", "TE30", "E", "edge", "TE", "TE_m3_n0"),
    MixedRequiredMode("E / edge / TE", "TE01", "E", "edge", "TE", "TE_m0_n1"),
    MixedRequiredMode("E / edge / TE", "TE11", "E", "edge", "TE", "TE_m1_n1"),
    MixedRequiredMode("E / edge / TE", "TE21", "E", "edge", "TE", "TE_m2_n1"),
    MixedRequiredMode("E / edge / TE", "TE31", "E", "edge", "TE", "TE_m3_n1"),
    MixedRequiredMode("E / edge / TE", "TE02", "E", "edge", "TE", "TE_m0_n2"),
    MixedRequiredMode("H / scalar / TE", "TE10", "H", "scalar", "TE", "TE_m1_n0"),
    MixedRequiredMode("H / scalar / TE", "TE20", "H", "scalar", "TE", "TE_m2_n0"),
    MixedRequiredMode("H / scalar / TE", "TE30", "H", "scalar", "TE", "TE_m3_n0"),
    MixedRequiredMode("H / scalar / TE", "TE01", "H", "scalar", "TE", "TE_m0_n1"),
    MixedRequiredMode("H / scalar / TE", "TE11", "H", "scalar", "TE", "TE_m1_n1"),
    MixedRequiredMode("H / scalar / TE", "TE21", "H", "scalar", "TE", "TE_m2_n1"),
    MixedRequiredMode("H / scalar / TE", "TE31", "H", "scalar", "TE", "TE_m3_n1"),
    MixedRequiredMode("H / scalar / TE", "TE02", "H", "scalar", "TE", "TE_m0_n2"),
    MixedRequiredMode("E / scalar / TM", "TM11", "E", "scalar", "TM", "TM_m1_n1"),
    MixedRequiredMode("E / scalar / TM", "TM21", "E", "scalar", "TM", "TM_m2_n1"),
    MixedRequiredMode("E / scalar / TM", "TM31", "E", "scalar", "TM", "TM_m3_n1"),
    MixedRequiredMode("E / scalar / TM", "TM02", "E", "edge", "TE", "TE_m0_n2"),
    MixedRequiredMode("H / edge / TM", "TM11", "H", "edge", "TM", "TM_m1_n1"),
    MixedRequiredMode("H / edge / TM", "TM21", "H", "edge", "TM", "TM_m2_n1"),
    MixedRequiredMode("H / edge / TM", "TM31", "H", "edge", "TM", "TM_m3_n1"),
    MixedRequiredMode("H / edge / TM", "TM02", "H", "scalar", "TE", "TE_m0_n2"),
)

MIXED_CIRCLE_MODES = (
    MixedRequiredMode("E / edge / TE", "TE01", "E", "edge", "TE", "TE_m0_p1"),
    MixedRequiredMode("E / edge / TE", "TE11", "E", "edge", "TE", "TE_m1_p1"),
    MixedRequiredMode("E / edge / TE", "TE12", "E", "edge", "TE", "TE_m2_p1"),
    MixedRequiredMode("E / edge / TE", "TE13", "E", "edge", "TE", "TE_m3_p1"),
    MixedRequiredMode("E / edge / TE", "TE20", "E", "edge", "TE", "TE_m0_p2"),
    MixedRequiredMode("E / edge / TE", "TE21", "E", "edge", "TE", "TE_m1_p2"),
    MixedRequiredMode("E / edge / TE", "TE22", "E", "edge", "TE", "TE_m2_p2"),
    MixedRequiredMode("E / edge / TE", "TE23", "E", "edge", "TE", "TE_m3_p2"),
    MixedRequiredMode("H / scalar / TE", "TE01", "H", "scalar", "TE", "TE_m0_p1"),
    MixedRequiredMode("H / scalar / TE", "TE11", "H", "scalar", "TE", "TE_m1_p1"),
    MixedRequiredMode("H / scalar / TE", "TE12", "H", "scalar", "TE", "TE_m2_p1"),
    MixedRequiredMode("H / scalar / TE", "TE13", "H", "scalar", "TE", "TE_m3_p1"),
    MixedRequiredMode("H / scalar / TE", "TE20", "H", "scalar", "TE", "TE_m0_p2"),
    MixedRequiredMode("H / scalar / TE", "TE21", "H", "scalar", "TE", "TE_m1_p2"),
    MixedRequiredMode("H / scalar / TE", "TE22", "H", "scalar", "TE", "TE_m2_p2"),
    MixedRequiredMode("H / scalar / TE", "TE23", "H", "scalar", "TE", "TE_m3_p2"),
    MixedRequiredMode("E / scalar / TM", "TM10", "E", "scalar", "TM", "TM_m0_p1"),
    MixedRequiredMode("E / scalar / TM", "TM11", "E", "scalar", "TM", "TM_m1_p1"),
    MixedRequiredMode("E / scalar / TM", "TM12", "E", "scalar", "TM", "TM_m2_p1"),
    MixedRequiredMode("E / scalar / TM", "TM13", "E", "scalar", "TM", "TM_m3_p1"),
    MixedRequiredMode("E / scalar / TM", "TM20", "E", "scalar", "TM", "TM_m0_p2"),
    MixedRequiredMode("E / scalar / TM", "TM21", "E", "scalar", "TM", "TM_m1_p2"),
    MixedRequiredMode("E / scalar / TM", "TM22", "E", "scalar", "TM", "TM_m2_p2"),
    MixedRequiredMode("E / scalar / TM", "TM23", "E", "scalar", "TM", "TM_m3_p2"),
    MixedRequiredMode("H / edge / TM", "TM10", "H", "edge", "TM", "TM_m0_p1"),
    MixedRequiredMode("H / edge / TM", "TM11", "H", "edge", "TM", "TM_m1_p1"),
    MixedRequiredMode("H / edge / TM", "TM12", "H", "edge", "TM", "TM_m2_p1"),
    MixedRequiredMode("H / edge / TM", "TM13", "H", "edge", "TM", "TM_m3_p1"),
    MixedRequiredMode("H / edge / TM", "TM20", "H", "edge", "TM", "TM_m0_p2"),
    MixedRequiredMode("H / edge / TM", "TM21", "H", "edge", "TM", "TM_m1_p2"),
    MixedRequiredMode("H / edge / TM", "TM22", "H", "edge", "TM", "TM_m2_p2"),
    MixedRequiredMode("H / edge / TM", "TM23", "H", "edge", "TM", "TM_m3_p2"),
)


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"CSV vazio: {path}")
    return rows


def _pick_best(rows: list[dict[str, str]]) -> dict[str, str]:
    return min(
        rows,
        key=lambda row: (
            float(row.get("error_percent", row.get("error_percent_analytic", "0")) or 0.0),
            -float(row.get("rho_abs", "0") or 0.0),
            int(row.get("positive_rank", "999999") or 999999),
        ),
    )


def _mode_label_rank(row: dict[str, str]) -> int:
    try:
        return int(row.get("positive_rank", "999999") or 999999)
    except (TypeError, ValueError):
        return 999999


def _parse_polar_mode_label(mode_label: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"(TE|TM)_m(\d+)_p(\d+)", mode_label)
    if not match:
        raise SystemExit(f"Rotulo modal polar invalido: {mode_label!r}")
    family, m_raw, p_raw = match.groups()
    return family, int(m_raw), int(p_raw)


def _polar_article_mode_label(family: str, mode_label: str) -> str:
    _, m, p = _parse_polar_mode_label(mode_label)
    if family == "TE" and m == 0 and p == 1:
        return "TE01"
    return f"{family}{p}{m}"


def _parse_scalar_polar_available(csv_path: Path) -> list[ModeRow]:
    rows = _read_csv_rows(csv_path)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["family"], row["mode_label"]), []).append(row)

    family_order = {"TE": 0, "TM": 1}
    ordered_keys = sorted(
        grouped.keys(),
        key=lambda key: (
            family_order.get(key[0], 999),
            _parse_polar_mode_label(key[1])[2],
            _parse_polar_mode_label(key[1])[1],
            min(_mode_label_rank(row) for row in grouped[key]),
        ),
    )

    out: list[ModeRow] = []
    for family, mode_label in ordered_keys:
        selected = _pick_best(grouped[(family, mode_label)])
        out.append(
            ModeRow(
                group=family,
                article_mode=_polar_article_mode_label(family, mode_label),
                repo_mode_label=selected["mode_label"],
                family=selected["family"],
                variable="kc",
                value=float(selected["kc_fem"]),
                ref_primary=float(selected["kc_ana"]),
                ref_secondary=None,
                err_primary_pct=float(selected["error_percent"]),
                err_secondary_pct=None,
                rho_label="rho_abs",
                rho_value=float(selected["rho_abs"]),
            )
        )
    return out


def _parse_scalar_case(csv_path: Path, required_modes: tuple[ScalarRequiredMode, ...]) -> list[ModeRow]:
    rows = _read_csv_rows(csv_path)
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["family"], row["mode_label"]), []).append(row)

    out: list[ModeRow] = []
    for spec in required_modes:
        selected = _pick_best(grouped[(spec.family, spec.mode_label)])
        out.append(
            ModeRow(
                group=spec.family,
                article_mode=spec.article_mode,
                repo_mode_label=selected["mode_label"],
                family=selected["family"],
                variable="kc",
                value=float(selected["kc_fem"]),
                ref_primary=float(selected["kc_ana"]),
                ref_secondary=None,
                err_primary_pct=float(selected["error_percent"]),
                err_secondary_pct=None,
                rho_label="rho_abs",
                rho_value=float(selected["rho_abs"]),
            )
        )
    return out


def _parse_mixed_case(csv_path: Path, required_modes: tuple[MixedRequiredMode, ...]) -> list[ModeRow]:
    rows = _read_csv_rows(csv_path)
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    for row in rows:
        key = (
            row["formulation"],
            row["dominant_block"],
            row["family"],
            row["mode_label"],
        )
        grouped.setdefault(key, []).append(row)

    out: list[ModeRow] = []
    for spec in required_modes:
        key = (spec.formulation, spec.dominant_block, spec.family, spec.mode_label)
        selected = _pick_best(grouped[key])
        out.append(
            ModeRow(
                group=spec.group,
                article_mode=spec.article_mode,
                repo_mode_label=selected["mode_label"],
                family=selected["family"],
                variable="kc",
                value=float(selected["kc_fem"]),
                ref_primary=float(selected["kc_ana"]),
                ref_secondary=None,
                err_primary_pct=float(selected["error_percent"]),
                err_secondary_pct=None,
                rho_label="rho_abs",
                rho_value=float(selected["rho_abs"]),
            )
        )
    return out


def parse_helm10_rect(case_root: Path) -> list[ModeRow]:
    return _parse_scalar_case(case_root / "csv" / "helm10_rect_modes.csv", SCALAR_RECT_MODES)


def parse_helm10_circle(case_root: Path) -> list[ModeRow]:
    return _parse_scalar_polar_available(case_root / "csv" / "helm10_circle_modes.csv")


def parse_helm10_coax(case_root: Path) -> list[ModeRow]:
    return _parse_scalar_case(case_root / "csv" / "helm10_coax_modes.csv", SCALAR_COAX_MODES)


def parse_mixed_rect(case_root: Path) -> list[ModeRow]:
    return _parse_mixed_case(case_root / "csv" / "mixed_rect_modes.csv", MIXED_RECT_MODES)


def parse_mixed_circle(case_root: Path) -> list[ModeRow]:
    return _parse_mixed_case(case_root / "csv" / "mixed_circle_modes.csv", MIXED_CIRCLE_MODES)


def parse_helmvec2_rect(case_root: Path) -> list[ModeRow]:
    rows = _read_csv_rows(case_root / "csv" / "helmvec2_rect_modes.csv")
    out: list[ModeRow] = []
    for row in rows:
        out.append(
            ModeRow(
                group="Figura 11 / Tabela 8",
                article_mode=f"modo {row['mode']}",
                repo_mode_label=f"cand{row['matched_candidate_rank']}",
                family="coupled",
                variable="k0L",
                value=float(row["k0l_fem_matched"]),
                ref_primary=float(row["ref_helmvec2"]),
                ref_secondary=float(row["ref_hayata"]),
                err_primary_pct=float(row["error_percent_helmvec2"]),
                err_secondary_pct=float(row["error_percent_hayata"]),
                rho_label="ez_ratio",
                rho_value=float(row["ez_ratio"]),
            )
        )
    return out


def parse_helmvec3_fig12(case_root: Path) -> list[ModeRow]:
    rows = _read_csv_rows(case_root / "csv" / "helmvec3_fig12_rect_table9.csv")
    out: list[ModeRow] = []
    for row in rows:
        out.append(
            ModeRow(
                group="Figura 12 / Tabela 9",
                article_mode=f"br/lambda0={row['br_over_lambda0']}",
                repo_mode_label=f"cand{row['selected_candidate_rank']}",
                family="coupled",
                variable="beta/k0",
                value=float(row["beta_over_k0_fem"]),
                ref_primary=float(row["beta_over_k0_analytic"]),
                ref_secondary=float(row["beta_over_k0_helmvec3"]),
                err_primary_pct=float(row["error_percent_analytic"]),
                err_secondary_pct=float(row["error_percent_helmvec3"]),
                rho_label="ez_ratio",
                rho_value=float(row["ez_ratio"]),
            )
        )
    return out


def parse_helmvec3_fig13(case_root: Path) -> list[ModeRow]:
    rows = _read_csv_rows(case_root / "csv" / "helmvec3_fig13_rect_table10.csv")
    out: list[ModeRow] = []
    for row in rows:
        out.append(
            ModeRow(
                group="Figura 13 / Tabela 10",
                article_mode=f"d/a={row['d_over_a']}, br/lambda0={row['br_over_lambda0']}",
                repo_mode_label=f"cand{row['selected_candidate_rank']}",
                family="coupled",
                variable="beta/k0",
                value=float(row["beta_over_k0_fem_matched"]),
                ref_primary=float(row["beta_over_k0_analytic"]),
                ref_secondary=float(row["beta_over_k0_helmvec3"]),
                err_primary_pct=float(row["error_percent_analytic"]),
                err_secondary_pct=float(row["error_percent_helmvec3"]),
                rho_label="ez_ratio",
                rho_value=float(row["ez_ratio"]),
            )
        )
    return out


CASES = (
    CaseSpec(
        key="helm10_rect",
        title="HELM10 retangular",
        section="2.1",
        article_ref="Tabela 1",
        exe="helm10_rect",
        out_rel=Path("helm10/rect"),
        base_args=("--ar-m", "1.0", "--nx", "10", "--ny", "20", "--nmodos", "10"),
        doubled_args=("--ar-m", "1.0", "--nx", "20", "--ny", "40", "--nmodos", "10"),
        primary_ref_label="Analitico (repositorio)",
        secondary_ref_label=None,
        parse_modes=parse_helm10_rect,
    ),
    CaseSpec(
        key="helm10_circle",
        title="HELM10 circular",
        section="2.1",
        article_ref="Tabela 2",
        exe="helm10_circle",
        out_rel=Path("helm10/circle"),
        base_args=("--nr", "8", "--nt", "15", "--nmodos", "10"),
        doubled_args=("--nr", "16", "--nt", "30", "--nmodos", "10"),
        primary_ref_label="Analitico (repositorio)",
        secondary_ref_label=None,
        parse_modes=parse_helm10_circle,
    ),
    CaseSpec(
        key="helm10_coax",
        title="HELM10 coaxial",
        section="2.1",
        article_ref="Tabela 3",
        exe="helm10_coax",
        out_rel=Path("helm10/coax"),
        base_args=("--nr", "10", "--nt", "17", "--nmodos", "10"),
        doubled_args=("--nr", "20", "--nt", "34", "--nmodos", "10"),
        primary_ref_label="Analitico (repositorio)",
        secondary_ref_label=None,
        parse_modes=parse_helm10_coax,
    ),
    CaseSpec(
        key="mixed_rect",
        title="HELMVEC1 retangular misto",
        section="2.2.2",
        article_ref="Tabela 6",
        exe="mixed_rect",
        out_rel=Path("helmvec1/rect"),
        base_args=("--nx", "10", "--ny", "20"),
        doubled_args=("--nx", "20", "--ny", "40"),
        primary_ref_label="Analitico (repositorio)",
        secondary_ref_label=None,
        parse_modes=parse_mixed_rect,
    ),
    CaseSpec(
        key="mixed_circle",
        title="HELMVEC1 circular misto",
        section="2.2.2",
        article_ref="Tabela 7",
        exe="mixed_circle",
        out_rel=Path("helmvec1/circle"),
        base_args=("--nr", "8", "--nt", "15"),
        doubled_args=("--nr", "16", "--nt", "30"),
        primary_ref_label="Analitico (repositorio)",
        secondary_ref_label=None,
        parse_modes=parse_mixed_circle,
    ),
    CaseSpec(
        key="helmvec2_rect",
        title="HELMVEC2 retangular parcialmente preenchido",
        section="2.2.3",
        article_ref="Figura 11 / Tabela 8",
        exe="helmvec2_rect",
        out_rel=Path("helmvec2/rect"),
        base_args=("--beta", "10", "--nx", "20", "--ny", "20"),
        doubled_args=("--beta", "10", "--nx", "40", "--ny", "40"),
        primary_ref_label="HELMVEC2 (artigo)",
        secondary_ref_label="Hayata (artigo)",
        parse_modes=parse_helmvec2_rect,
    ),
    CaseSpec(
        key="helmvec3_fig12_rect",
        title="HELMVEC3 retangular parcialmente preenchido",
        section="2.2.4",
        article_ref="Figura 12 / Tabela 9",
        exe="helmvec3_fig12_rect",
        out_rel=Path("helmvec3/fig12_rect"),
        base_args=("--nx", "10", "--ny", "5"),
        doubled_args=("--nx", "20", "--ny", "10"),
        primary_ref_label="Analitico",
        secondary_ref_label="HELMVEC3 (artigo)",
        parse_modes=parse_helmvec3_fig12,
    ),
    CaseSpec(
        key="helmvec3_fig13_rect",
        title="HELMVEC3 retangular parcialmente preenchido",
        section="2.2.4",
        article_ref="Figura 13 / Tabela 10",
        exe="helmvec3_fig13_rect",
        out_rel=Path("helmvec3/fig13_rect"),
        base_args=("--d-over-a-preview", "0.2", "--nx", "10", "--ny", "5"),
        doubled_args=("--d-over-a-preview", "0.2", "--nx", "20", "--ny", "10"),
        primary_ref_label="Analitico",
        secondary_ref_label="HELMVEC3 (artigo)",
        parse_modes=parse_helmvec3_fig13,
    ),
)


def _run_case(build_dir: Path, out_root: Path, case: CaseSpec, preset: MeshPreset, backend_label: str, backend: str) -> None:
    case_root = out_root / preset.key / backend_label.lower() / case.out_rel
    timing_csv = case_root / "run_timing.csv"
    if timing_csv.exists():
        print(f"[reuse] {preset.key} {backend_label} {case.key}")
        return

    case_root.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(build_dir / case.exe)]
    cmd.extend(case.base_args if preset.key == BASE.key else case.doubled_args)
    cmd.extend(["--backend", backend])
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str((out_root / preset.key / backend_label.lower()).resolve())
    print(f"[run] {preset.key:<7} backend={backend_label:<4} case={case.key:<18} cmd={' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, text=True, capture_output=True)
    log_path = out_root / preset.key / backend_label.lower() / f"{case.key}.stdout.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(proc.stdout + "\n\n[stderr]\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise SystemExit(f"Falha ao executar {' '.join(cmd)}\nVeja: {log_path}")


def _read_timing(path: Path) -> dict[str, str]:
    row = _read_csv_rows(path)[0]
    return {
        "assembly_ms": row.get("assembly_ms_total") or row.get("assembly_ms") or "",
        "solve_ms": row.get("solve_ms_total") or row.get("solve_ms") or "",
        "post_ms": row.get("post_ms_total") or row.get("post_ms") or "",
        "total_ms": row.get("total_ms") or "",
        "mesh_nodes": row.get("mesh_nodes") or "",
        "mesh_tris": row.get("mesh_tris") or "",
    }


def _count_edges_from_vtk(vtk_path: Path) -> tuple[int, int, int]:
    n_points = 0
    n_cells = 0
    edges: set[tuple[int, int]] = set()
    lines = vtk_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    for idx, line in enumerate(lines):
        parts = line.strip().split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            n_points = int(parts[1])
        if len(parts) >= 3 and parts[0] == "CELLS":
            n_cells = int(parts[1])
            cursor = idx + 1
            for _ in range(n_cells):
                cell_parts = [int(tok) for tok in lines[cursor].strip().split()]
                cursor += 1
                if not cell_parts:
                    continue
                nn = cell_parts[0]
                nodes = cell_parts[1 : 1 + nn]
                if nn >= 2:
                    for a in range(nn):
                        b = (a + 1) % nn
                        edge = tuple(sorted((nodes[a], nodes[b])))
                        edges.add(edge)
            break
    return n_points, n_cells, len(edges)


def _mesh_info(case_root: Path, timing: dict[str, str]) -> dict[str, int]:
    vtk_dir = case_root / "vtk"
    vtk_files = sorted(vtk_dir.glob("*.vtk"))
    if not vtk_files:
        raise SystemExit(f"Nenhum VTK encontrado em {vtk_dir}")
    vtk_nodes, vtk_tris, vtk_edges = _count_edges_from_vtk(vtk_files[0])
    nodes = int(timing["mesh_nodes"]) if timing["mesh_nodes"] else vtk_nodes
    tris = int(timing["mesh_tris"]) if timing["mesh_tris"] else vtk_tris
    return {
        "nodes": nodes,
        "tris": tris,
        "edges": vtk_edges,
    }


def _fmt_num(value: float | None, digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def _fmt_time(raw: str) -> str:
    if not raw:
        return ""
    return f"{float(raw):.3f}"


def _group_order(rows: list[ModeRow]) -> list[str]:
    seen: set[str] = set()
    order: list[str] = []
    for row in rows:
        if row.group not in seen:
            seen.add(row.group)
            order.append(row.group)
    return order


def _mode_index(rows: list[ModeRow]) -> dict[tuple[str, str], ModeRow]:
    return {(row.group, row.article_mode): row for row in rows}


def _write_report(out_path: Path, out_root: Path, presets: tuple[MeshPreset, ...]) -> None:
    lines: list[str] = []
    lines.append("# Relatorio Comparativo FEM x EFGM")
    lines.append("")
    lines.append(f"Gerado em: `{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}`")
    lines.append("")
    lines.append("Observacao: neste relatorio, `EFGM` corresponde ao backend `efgmi` atualmente usado no repositorio.")
    lines.append("Nos casos mistos/acoplados, isso significa `base nodal EFGM + base vetorial FEM`.")
    lines.append("As contagens de nos/triangulos/arestas abaixo sao lidas dos artefatos realmente gerados por cada executavel.")
    lines.append("Nos casos circulares com `--nr 8 --nt 15`, a malha polar atual do repositorio gera `225` triangulos.")
    lines.append("")

    for preset in presets:
        lines.append(f"## {preset.title}")
        lines.append("")
        lines.append("| Caso | Secao | Tabela/Figura | Discretizacao | Nos | Triangulos | Arestas |")
        lines.append("|---|---|---|---|---:|---:|---:|")

        case_cache: dict[str, dict[str, object]] = {}
        for case in CASES:
            fem_root = out_root / preset.key / "fem" / case.out_rel
            fem_timing = _read_timing(fem_root / "run_timing.csv")
            info = _mesh_info(fem_root, fem_timing)
            args = case.base_args if preset.key == BASE.key else case.doubled_args
            discret = " ".join(args)
            case_cache[case.key] = {
                "info": info,
                "discret": discret,
            }
            lines.append(
                f"| {case.title} | {case.section} | {case.article_ref} | `{discret}` | "
                f"{info['nodes']} | {info['tris']} | {info['edges']} |"
            )
        lines.append("")

        for case in CASES:
            lines.append(f"### {case.title} — {case.article_ref}")
            lines.append("")
            info = case_cache[case.key]["info"]
            discret = case_cache[case.key]["discret"]
            lines.append("| Campo | Valor |")
            lines.append("|---|---|")
            lines.append(f"| Executavel | `{case.exe}` |")
            lines.append(f"| Secao | `{case.section}` |")
            lines.append(f"| Discretizacao | `{discret}` |")
            lines.append(f"| Nos nodais | `{info['nodes']}` |")
            lines.append(f"| Triangulos | `{info['tris']}` |")
            lines.append(f"| Arestas | `{info['edges']}` |")
            lines.append(f"| Referencia primaria | `{case.primary_ref_label}` |")
            lines.append(f"| Referencia secundaria | `{case.secondary_ref_label or '-'}` |")
            lines.append("")

            fem_root = out_root / preset.key / "fem" / case.out_rel
            efgm_root = out_root / preset.key / "efgm" / case.out_rel
            fem_timing = _read_timing(fem_root / "run_timing.csv")
            efgm_timing = _read_timing(efgm_root / "run_timing.csv")
            lines.append("| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |")
            lines.append("|---|---:|---:|---:|---:|")
            lines.append(
                f"| FEM | {_fmt_time(fem_timing['assembly_ms'])} | {_fmt_time(fem_timing['solve_ms'])} | "
                f"{_fmt_time(fem_timing['post_ms'])} | {_fmt_time(fem_timing['total_ms'])} |"
            )
            lines.append(
                f"| EFGM | {_fmt_time(efgm_timing['assembly_ms'])} | {_fmt_time(efgm_timing['solve_ms'])} | "
                f"{_fmt_time(efgm_timing['post_ms'])} | {_fmt_time(efgm_timing['total_ms'])} |"
            )
            lines.append("")

            fem_modes = case.parse_modes(fem_root)
            efgm_modes = case.parse_modes(efgm_root)
            fem_idx = _mode_index(fem_modes)
            efgm_idx = _mode_index(efgm_modes)
            for group in _group_order(fem_modes):
                lines.append(f"#### {group}")
                lines.append("")
                lines.append(
                    "| Modo artigo | Label repo | Familia | Ref. primaria | Ref. secundaria | "
                    "FEM valor | FEM err p.(%) | FEM err s.(%) | FEM rho | "
                    "EFGM valor | EFGM err p.(%) | EFGM err s.(%) | EFGM rho |"
                )
                lines.append(
                    "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
                )
                group_rows = [row for row in fem_modes if row.group == group]
                for row in group_rows:
                    fem_row = fem_idx[(group, row.article_mode)]
                    efgm_row = efgm_idx[(group, row.article_mode)]
                    lines.append(
                        f"| {row.article_mode} | `{row.repo_mode_label}` | {row.family} | "
                        f"{_fmt_num(fem_row.ref_primary)} | {_fmt_num(fem_row.ref_secondary)} | "
                        f"{_fmt_num(fem_row.value)} | {_fmt_num(fem_row.err_primary_pct, 3)} | {_fmt_num(fem_row.err_secondary_pct, 3)} | {_fmt_num(fem_row.rho_value, 6)} | "
                        f"{_fmt_num(efgm_row.value)} | {_fmt_num(efgm_row.err_primary_pct, 3)} | {_fmt_num(efgm_row.err_secondary_pct, 3)} | {_fmt_num(efgm_row.rho_value, 6)} |"
                    )
                lines.append("")

        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Gera relatorio FEM x EFGM por modo e por tabela do artigo.")
    ap.add_argument("--build-dir", type=Path, default=Path("build"))
    ap.add_argument("--out-dir", type=Path, default=Path("out/fem_efgmi_mode_report"))
    ap.add_argument("--report", type=Path, default=Path("out/fem_efgmi_mode_report/FEM_EFGM_MODE_REPORT.md"))
    ap.add_argument(
        "--presets",
        default="base",
        help="Lista separada por virgula entre: base,doubled (default: base).",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    build_dir = (ROOT / args.build_dir).resolve()
    out_root = (ROOT / args.out_dir).resolve()
    report_path = (ROOT / args.report).resolve()
    preset_keys = tuple(part.strip() for part in args.presets.split(",") if part.strip())
    if not preset_keys:
        raise SystemExit("Nenhum preset informado em --presets.")
    try:
        presets = tuple(PRESET_BY_KEY[key] for key in preset_keys)
    except KeyError as exc:
        raise SystemExit(f"Preset invalido em --presets: {exc.args[0]}") from exc

    for preset in presets:
        for backend_label, backend in BACKENDS:
            for case in CASES:
                _run_case(build_dir, out_root, case, preset, backend_label, backend)

    _write_report(report_path, out_root, presets)
    print(f"Saved: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
