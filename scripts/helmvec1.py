#!/usr/bin/env python3
"""
Generate HELMVEC1 summary images directly from the modal CSV outputs.

The HELMVEC1 family currently exports didactic modal tables, not spatial field
maps. This script reads:

- `out/helmvec1/<case>/csv/mixed_<case>_modes.csv`

and produces summary images per case:

- normalized cutoff by mode
- rho by mode
- dominant energy ratio by mode
- edge energy vs scalar energy by mode
- error by mode, when analytic references are available
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]

GROUP_ORDER: Sequence[Tuple[str, str]] = (
    ("E", "TE"),
    ("E", "TM"),
    ("H", "TE"),
    ("H", "TM"),
)

GROUP_COLORS: Dict[Tuple[str, str], str] = {
    ("E", "TE"): "#1f77b4",
    ("E", "TM"): "#ff7f0e",
    ("H", "TE"): "#2ca02c",
    ("H", "TM"): "#d62728",
}


@dataclass
class ModeRow:
    case_name: str
    row: Dict[str, str]

    @property
    def formulation(self) -> str:
        return self.row["formulation"]

    @property
    def dominant_block(self) -> str:
        return self.row["dominant_block"]

    @property
    def component_label(self) -> str:
        return self.row["component_label"]

    @property
    def family(self) -> str:
        return self.row["family"]

    @property
    def mode_label(self) -> str:
        return self.row["mode_label"]

    @property
    def positive_rank(self) -> int:
        return int(self.row["positive_rank"])

    @property
    def eig_index(self) -> int:
        return int(self.row["eig_index"])

    @property
    def kc_fem(self) -> Optional[float]:
        return _to_float(self.row.get("kc_fem"))

    @property
    def normalized_kc(self) -> Optional[float]:
        if self.case_name == "rect":
            value = _to_float(self.row.get("kc_ar_fem"))
            if value is not None:
                return value
            kc = self.kc_fem
            a = _to_float(self.row.get("ar_m"))
            if kc is not None and a is not None:
                return kc * a
            return None

        if self.case_name == "circle":
            value = _to_float(self.row.get("kc_r_fem"))
            if value is not None:
                return value
            kc = self.kc_fem
            r = _to_float(self.row.get("r_m"))
            if kc is not None and r is not None:
                return kc * r
            return None

        value = _to_float(self.row.get("kc_r1_fem"))
        if value is not None:
            return value
        kc = self.kc_fem
        r1 = _to_float(self.row.get("r1_m"))
        if kc is not None and r1 is not None:
            return kc * r1
        return None

    @property
    def normalized_kc_label(self) -> str:
        if self.case_name == "rect":
            return "$k_c a_r$"
        if self.case_name == "circle":
            return "$k_c r$"
        return "$k_c r_1$"

    @property
    def error_percent(self) -> Optional[float]:
        return _to_float(self.row.get("error_percent"))

    @property
    def rho_abs(self) -> Optional[float]:
        return _to_float(self.row.get("rho_abs"))

    @property
    def edge_energy(self) -> float:
        return float(self.row["edge_energy"])

    @property
    def scalar_energy(self) -> float:
        return float(self.row["scalar_energy"])

    @property
    def dominant_energy_ratio(self) -> float:
        return float(self.row["dominant_energy_ratio"])

    def group_key(self) -> Tuple[str, str]:
        return (self.formulation, self.family)

    def subplot_title(self) -> str:
        return (
            f"{self.case_name.upper()} | {self.formulation} / {self.family} / "
            f"{self.component_label} / {self.dominant_block}"
        )

    def compact_mode_label(self) -> str:
        label = self.mode_label
        prefix = f"{self.family}_"
        if label.startswith(prefix):
            label = label[len(prefix):]
        return label.replace("_", ",")


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _read_rows(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _read_mode_rows(case_name: str, case_root: Path) -> List[ModeRow]:
    csv_path = case_root / "csv" / f"mixed_{case_name}_modes.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {csv_path}")
    return [ModeRow(case_name=case_name, row=row) for row in _read_rows(csv_path)]


def _sort_modes(rows: Iterable[ModeRow]) -> List[ModeRow]:
    return sorted(
        rows,
        key=lambda item: (
            float("inf") if item.kc_fem is None else item.kc_fem,
            item.positive_rank,
            item.mode_label,
        ),
    )


def _group_rows(rows: Sequence[ModeRow]) -> List[Tuple[Tuple[str, str], List[ModeRow]]]:
    grouped: List[Tuple[Tuple[str, str], List[ModeRow]]] = []
    by_key: Dict[Tuple[str, str], List[ModeRow]] = {}
    for row in rows:
        by_key.setdefault(row.group_key(), []).append(row)

    for key in GROUP_ORDER:
        if key in by_key:
            grouped.append((key, _sort_modes(by_key[key])))

    for key in sorted(by_key):
        if key not in GROUP_ORDER:
            grouped.append((key, _sort_modes(by_key[key])))

    return grouped


def _subplot_grid(num_groups: int) -> Tuple[int, int]:
    cols = 2 if num_groups > 1 else 1
    rows = int(math.ceil(num_groups / cols))
    return rows, cols


def _save_figure(fig: plt.Figure, out_path: Path, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0.03, 0.03, 0.98, 0.95))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _annotate_points(ax: plt.Axes, xpos: np.ndarray, values: Sequence[float], fmt: str) -> None:
    if not values:
        return
    ymin, ymax = ax.get_ylim()
    span = max(1e-12, ymax - ymin)
    offset = 0.02 * span
    for x, value in zip(xpos, values):
        ax.text(
            x,
            value + offset,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=8,
        )


def _is_constant_series(values: Sequence[float], tol: float = 1e-12) -> bool:
    if not values:
        return False
    first = values[0]
    return all(abs(value - first) <= tol for value in values)


def _plot_line_summary(
    case_name: str,
    rows: Sequence[ModeRow],
    out_path: Path,
    dpi: int,
    *,
    figure_title: str,
    y_label: str,
    value_getter,
    fmt: str,
    y_limits: Optional[Tuple[float, float]] = None,
    empty_warning: str,
    constant_note: Optional[str] = None,
) -> None:
    grouped = _group_rows(rows)
    plotted_groups = [(key, group) for key, group in grouped if any(value_getter(item) is not None for item in group)]
    if not plotted_groups:
        print(f"Warning: {case_name} nao possui dados de {empty_warning}.")
        return

    nrows, ncols = _subplot_grid(len(plotted_groups))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(max(12.0, 5.8 * ncols), max(4.3 * nrows, 4.5)),
        squeeze=False,
    )
    axes_flat = list(axes.flat)

    for ax, (key, group) in zip(axes_flat, plotted_groups):
        color = GROUP_COLORS.get(key, "#1f77b4")
        valid = [(item, value_getter(item)) for item in group]
        valid = [(item, value) for item, value in valid if value is not None]
        xpos = np.arange(len(valid))
        labels = [item.compact_mode_label() for item, _ in valid]
        values = [float(value) for _, value in valid]
        ax.plot(xpos, values, color=color, marker="o", lw=1.8, ms=5.5)
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel(y_label)
        ax.set_title(valid[0][0].subplot_title())
        ax.grid(True, axis="y", alpha=0.28)
        if y_limits is not None:
            ax.set_ylim(*y_limits)
        if constant_note is not None and _is_constant_series(values):
            ax.text(
                0.5,
                0.88,
                f"{constant_note}: {fmt.format(values[0])}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                bbox={"boxstyle": "round,pad=0.28", "fc": "white", "ec": "0.75", "alpha": 0.95},
            )
        else:
            _annotate_points(ax, xpos, values, fmt)

    for ax in axes_flat[len(plotted_groups):]:
        ax.axis("off")

    fig.suptitle(f"HELMVEC1 | {figure_title} | {case_name.upper()}", fontsize=14)
    _save_figure(fig, out_path, dpi)


def _plot_block_energy_summary(
    case_name: str,
    rows: Sequence[ModeRow],
    out_path: Path,
    dpi: int,
) -> None:
    grouped = _group_rows(rows)
    if not grouped:
        print(f"Warning: {case_name} nao possui dados de energia para o resumo.")
        return

    nrows, ncols = _subplot_grid(len(grouped))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(max(12.0, 5.8 * ncols), max(4.5 * nrows, 4.8)),
        squeeze=False,
    )
    axes_flat = list(axes.flat)

    for ax, (key, group) in zip(axes_flat, grouped):
        color_edge = "#1f77b4"
        color_scalar = "#d62728"
        xpos = np.arange(len(group))
        labels = [item.compact_mode_label() for item in group]
        edge_values = [item.edge_energy for item in group]
        scalar_values = [item.scalar_energy for item in group]
        ax.plot(xpos, edge_values, color=color_edge, marker="o", lw=1.8, ms=5.2, label="edge_energy")
        ax.plot(
            xpos,
            scalar_values,
            color=color_scalar,
            marker="s",
            lw=1.8,
            ms=4.8,
            label="scalar_energy",
        )
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("energia numerica")
        ax.set_title(group[0].subplot_title())
        ax.grid(True, axis="y", alpha=0.28)
        ax.legend(loc="best", fontsize=8)

    for ax in axes_flat[len(grouped):]:
        ax.axis("off")

    fig.suptitle(
        f"HELMVEC1 | Energias de Bloco por Modo | {case_name.upper()}",
        fontsize=14,
    )
    _save_figure(fig, out_path, dpi)


def _plot_case_images(case_root: Path, case_name: str, dpi: int) -> None:
    print(f"Processing HELMVEC1 case: {case_name}")
    rows = _read_mode_rows(case_name=case_name, case_root=case_root)
    if not rows:
        print(f"Warning: {case_root} nao possui linhas em mixed_{case_name}_modes.csv.")
        return

    img_root = case_root / "img"
    img_root.mkdir(parents=True, exist_ok=True)

    normalized_label = rows[0].normalized_kc_label
    _plot_line_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_kc_by_mode.png",
        dpi=dpi,
        figure_title="Cutoff Normalizado por Modo",
        y_label=normalized_label,
        value_getter=lambda item: item.normalized_kc,
        fmt="{:.4f}",
        empty_warning="cutoff normalizado para o resumo",
    )
    _plot_line_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_rho_by_mode.png",
        dpi=dpi,
        figure_title="Correlacao Modal por Modo",
        y_label="|rho|",
        value_getter=lambda item: item.rho_abs,
        fmt="{:.4f}",
        y_limits=(0.0, 1.05),
        empty_warning="correlacao modal para o resumo",
        constant_note="todos os modos",
    )
    _plot_line_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_dominant_ratio_by_mode.png",
        dpi=dpi,
        figure_title="Razao de Energia Dominante por Modo",
        y_label="dominant_energy_ratio",
        value_getter=lambda item: item.dominant_energy_ratio,
        fmt="{:.4f}",
        y_limits=(0.0, 1.05),
        empty_warning="razao de energia dominante para o resumo",
        constant_note="todos os modos",
    )
    _plot_block_energy_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_block_energy_by_mode.png",
        dpi=dpi,
    )
    _plot_line_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_error_by_mode.png",
        dpi=dpi,
        figure_title="Erro Percentual por Modo",
        y_label="|error_percent| [%]",
        value_getter=lambda item: abs(item.error_percent) if item.error_percent is not None else None,
        fmt="{:.2f}",
        empty_warning="erro percentual para o resumo",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do HELMVEC1 e gera imagens-resumo do sistema misto "
            "da Eq. (92): cutoff normalizado, correlacao modal, razao "
            "dominante, energias de bloco e erro quando houver referencia "
            "analitica."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("out/helmvec1"),
        help="Diretorio raiz das saidas do HELMVEC1. Padrao: out/helmvec1",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=["rect", "circle", "coax", "all"],
        default=None,
        help="Seleciona quais casos processar. Pode ser repetido.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="Resolucao das imagens.")
    args = parser.parse_args()

    root = _resolve(args.root)
    if not root.exists():
        raise SystemExit(f"Diretorio de saida nao encontrado: {root}")

    raw_cases = args.case or ["all"]
    cases: List[str] = []
    for item in raw_cases:
        if item == "all":
            cases.extend(["rect", "circle", "coax"])
        else:
            cases.append(item)
    ordered_cases = list(dict.fromkeys(cases))

    for case_name in ordered_cases:
        case_root = root / case_name
        if not case_root.exists():
            print(f"Warning: caso ausente, pulando: {case_root}")
            continue
        _plot_case_images(case_root=case_root, case_name=case_name, dpi=args.dpi)


if __name__ == "__main__":
    main()
