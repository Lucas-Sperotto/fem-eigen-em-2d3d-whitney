#!/usr/bin/env python3
"""
Generate HELMVEC1 images directly from the mixed modal CSV outputs.

The HELMVEC1 family now exports:

- `out/helmvec1/<case>/csv/mixed_<case>_modes.csv`
- one spatial CSV per exported mode
- one VTK per exported mode

The script generates:

- summary figures (normalized cutoff, rho, dominant ratio, block energies, error)
- spatial figures for edge-dominant modes:
  - magnitude of the transverse vector field
  - quiver of the transverse vector field
- spatial figures for scalar-dominant modes:
  - scalar map with isolines of the longitudinal component
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
import matplotlib.tri as mtri
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
class EdgeFieldData:
    cell_ids: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    fx: np.ndarray
    fy: np.ndarray
    fmag: np.ndarray
    x_name: str
    y_name: str
    mag_name: str
    symbol: str


@dataclass
class ScalarFieldData:
    node_ids: np.ndarray
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    scalar_name: str


@dataclass
class ModeRow:
    case_name: str
    case_root: Path
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
    def m(self) -> Optional[int]:
        return _to_int(self.row.get("m"))

    @property
    def n(self) -> Optional[int]:
        return _to_int(self.row.get("n"))

    @property
    def p(self) -> Optional[int]:
        return _to_int(self.row.get("p"))

    @property
    def kc_fem(self) -> Optional[float]:
        return _to_float(self.row.get("kc_fem"))

    @property
    def field_data_kind(self) -> str:
        return self.row.get("field_data_kind", "").strip()

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "").strip()

    @property
    def fields_csv_name(self) -> str:
        return self.row.get("fields_csv_file", "").strip()

    @property
    def fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.fields_csv_name

    @property
    def vtk_name(self) -> str:
        return self.row.get("vtk_file", "").strip()

    @property
    def vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.vtk_name

    @property
    def has_spatial_artifacts(self) -> bool:
        return (
            bool(self.field_data_kind)
            and bool(self.fields_csv_name)
            and bool(self.vtk_name)
            and self.fields_csv_path.exists()
            and self.vtk_path.exists()
        )

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

    def index_pair(self) -> Tuple[Optional[int], Optional[int]]:
        if self.case_name == "rect":
            return self.m, self.n
        return self.m, self.p

    def mode_tex(self) -> str:
        i0, i1 = self.index_pair()
        if i0 is None or i1 is None:
            return self.family
        return f"{self.family}$_{{{i0},{i1}}}$"

    def mode_line(self) -> str:
        return (
            f"{self.case_name.upper()} | {self.formulation} | "
            f"{self.mode_tex()} | {self.component_label}"
        )

    def normalized_cutoff_line(self) -> str:
        value = self.normalized_kc
        if value is None:
            return f"{self.normalized_kc_label} = ?"
        return f"{self.normalized_kc_label} = {value:.6f}"

    def spatial_stem(self) -> str:
        return (
            f"helmvec1_{self.case_name}_{self.formulation}_{self.component_label}_"
            f"{self.mode_label}_rank{self.positive_rank:02d}"
        )


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


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _read_rows(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _read_mode_rows(case_name: str, case_root: Path) -> List[ModeRow]:
    csv_path = case_root / "csv" / f"mixed_{case_name}_modes.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {csv_path}")
    return [ModeRow(case_name=case_name, case_root=case_root, row=row) for row in _read_rows(csv_path)]


def _read_edge_field_data(csv_path: Path) -> EdgeFieldData:
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"{csv_path}: CSV de campos vetoriais vazio")
    rows.sort(key=lambda row: int(row["cell_id"]))

    if {"Ex", "Ey", "Emag"}.issubset(fieldnames):
        x_name, y_name, mag_name, symbol = "Ex", "Ey", "Emag", "E"
    elif {"Hx", "Hy", "Hmag"}.issubset(fieldnames):
        x_name, y_name, mag_name, symbol = "Hx", "Hy", "Hmag", "H"
    else:
        raise ValueError(f"{csv_path}: cabecalho vetorial desconhecido")

    def col(name: str) -> np.ndarray:
        return np.array([float(row[name]) for row in rows], dtype=float)

    return EdgeFieldData(
        cell_ids=np.array([int(row["cell_id"]) for row in rows], dtype=int),
        xc=col("xc_m"),
        yc=col("yc_m"),
        fx=col(x_name),
        fy=col(y_name),
        fmag=col(mag_name),
        x_name=x_name,
        y_name=y_name,
        mag_name=mag_name,
        symbol=symbol,
    )


def _read_scalar_field_data(csv_path: Path) -> ScalarFieldData:
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"{csv_path}: CSV de campos escalares vazio")

    scalar_columns = [name for name in fieldnames if name not in {"node_id", "x_m", "y_m"}]
    if len(scalar_columns) != 1:
        raise ValueError(f"{csv_path}: esperado exatamente um campo escalar, obtido {scalar_columns}")
    scalar_name = scalar_columns[0]

    rows.sort(key=lambda row: int(row["node_id"]))

    def col(name: str) -> np.ndarray:
        return np.array([float(row[name]) for row in rows], dtype=float)

    return ScalarFieldData(
        node_ids=np.array([int(row["node_id"]) for row in rows], dtype=int),
        x=col("x_m"),
        y=col("y_m"),
        values=col(scalar_name),
        scalar_name=scalar_name,
    )


def _read_legacy_vtk_connectivity(vtk_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    lines = vtk_path.read_text(encoding="utf-8").splitlines()

    points_idx = next((i for i, line in enumerate(lines) if line.strip().startswith("POINTS ")), -1)
    if points_idx < 0:
        raise ValueError(f"{vtk_path}: bloco POINTS ausente")
    num_points = int(lines[points_idx].split()[1])
    points = np.array(
        [[float(v) for v in lines[points_idx + 1 + k].split()[:2]] for k in range(num_points)],
        dtype=float,
    )

    cells_idx = next((i for i, line in enumerate(lines) if line.strip().startswith("CELLS ")), -1)
    if cells_idx < 0:
        raise ValueError(f"{vtk_path}: bloco CELLS ausente")
    num_cells = int(lines[cells_idx].split()[1])
    triangles: List[List[int]] = []
    for k in range(num_cells):
        items = [int(v) for v in lines[cells_idx + 1 + k].split()]
        if items[0] != 3:
            raise ValueError(f"{vtk_path}: apenas triangulos sao suportados")
        triangles.append(items[1:4])

    return points, np.array(triangles, dtype=int)


def _subplot_grid(num_groups: int) -> Tuple[int, int]:
    cols = 2 if num_groups > 1 else 1
    rows = int(math.ceil(num_groups / cols))
    return rows, cols


def _save_figure(fig: plt.Figure, out_path: Path, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0.03, 0.03, 0.98, 0.96))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _set_axes_style(ax: plt.Axes, points: np.ndarray) -> None:
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_xlim(points[:, 0].min(), points[:, 0].max())
    ax.set_ylim(points[:, 1].min(), points[:, 1].max())


def _thin_indices(n: int, max_arrows: int) -> np.ndarray:
    if n <= max_arrows:
        return np.arange(n, dtype=int)
    step = max(1, n // max_arrows)
    return np.arange(0, n, step, dtype=int)


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


def _sort_modes(rows: Iterable[ModeRow]) -> List[ModeRow]:
    return sorted(
        rows,
        key=lambda item: (
            float("inf") if item.kc_fem is None else item.kc_fem,
            item.positive_rank,
            item.mode_label,
        ),
    )


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

    for ax, (_, group) in zip(axes_flat, grouped):
        xpos = np.arange(len(group))
        labels = [item.compact_mode_label() for item in group]
        edge_values = [item.edge_energy for item in group]
        scalar_values = [item.scalar_energy for item in group]
        ax.plot(xpos, edge_values, color="#1f77b4", marker="o", lw=1.8, ms=5.2, label="edge_energy")
        ax.plot(xpos, scalar_values, color="#d62728", marker="s", lw=1.8, ms=4.8, label="scalar_energy")
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel("energia numerica")
        ax.set_title(group[0].subplot_title())
        ax.grid(True, axis="y", alpha=0.28)
        ax.legend(loc="best", fontsize=8)

    for ax in axes_flat[len(grouped):]:
        ax.axis("off")

    fig.suptitle(f"HELMVEC1 | Energias de Bloco por Modo | {case_name.upper()}", fontsize=14)
    _save_figure(fig, out_path, dpi)


def _plot_edge_magnitude(
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(8.2, 6.3))
    tpc = ax.tripcolor(tri, facecolors=field_data.fmag, shading="flat", cmap="viridis")
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="white", alpha=0.35, zorder=2)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label(f"|{field_data.symbol}|")
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC1 | Magnitude do Campo Transversal\n"
        f"{meta.mode_line()}\n"
        f"{meta.normalized_cutoff_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_edge_quiver(
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(8.2, 6.3))
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="0.82", zorder=0)

    keep = _thin_indices(len(field_data.cell_ids), max_arrows=max_arrows)
    ax.quiver(
        field_data.xc[keep],
        field_data.yc[keep],
        field_data.fx[keep],
        field_data.fy[keep],
        color="#154c79",
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0032,
        zorder=3,
    )

    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC1 | Distribuicao do Campo Transversal\n"
        f"{meta.mode_line()}\n"
        f"{meta.normalized_cutoff_line()} | setas=({field_data.x_name}, {field_data.y_name})"
    )
    _save_figure(fig, out_path, dpi)


def _plot_scalar_mode(
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: ScalarFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(8.2, 6.3))
    tpc = ax.tripcolor(tri, field_data.values, shading="gouraud", cmap="coolwarm")
    levels = 14
    ax.tricontour(tri, field_data.values, levels=levels, colors="k", linewidths=0.45, alpha=0.55)
    if show_mesh:
        ax.triplot(tri, lw=0.25, color="0.45", alpha=0.25, zorder=3)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label(field_data.scalar_name)
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC1 | Funcao Longitudinal Escalar\n"
        f"{meta.mode_line()}\n"
        f"{meta.normalized_cutoff_line()} | {field_data.scalar_name}"
    )
    _save_figure(fig, out_path, dpi)


def _collect_modal_summary_data(
    rows: Sequence[ModeRow],
    value_getter,
) -> List[Tuple[str, int, str, float, Optional[float]]]:
    data: List[Tuple[str, int, str, float, Optional[float]]] = []
    for row in rows:
        value = value_getter(row)
        if value is None:
            continue
        data.append((row.family, row.positive_rank, row.mode_label, float(value), row.kc_fem))
    return data


def _plot_modal_summary(
    case_name: str,
    rows: Sequence[ModeRow],
    out_path: Path,
    dpi: int,
    *,
    value_getter,
    y_label: str,
    family_title: str,
    figure_title: str,
    empty_warning: str,
    y_limits: Optional[Tuple[float, float]] = None,
    fmt: str = "{:.2f}",
) -> None:
    data = _collect_modal_summary_data(rows, value_getter=value_getter)
    if not data:
        print(f"Warning: {case_name} nao possui dados de {empty_warning}.")
        return

    families = sorted({family for family, _, _, _, _ in data})
    fig, axes = plt.subplots(
        len(families),
        1,
        figsize=(max(11.0, 0.72 * max(len(data), 6)), 4.2 * len(families)),
        squeeze=False,
    )

    for ax, family in zip(axes[:, 0], families):
        subset = sorted(
            (item for item in data if item[0] == family),
            key=lambda item: (
                float("inf") if item[4] is None else item[4],
                item[1],
                item[2],
            ),
        )
        labels = [mode_label for _, _, mode_label, _, _ in subset]
        values = [value for _, _, _, value, _ in subset]
        xpos = np.arange(len(subset))
        color = "#2c7fb8" if family == "TE" else "#d95f0e"
        ax.plot(xpos, values, color=color, marker="o", lw=1.8, ms=5.5)
        ax.set_xticks(xpos)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel(y_label)
        ax.set_title(f"{case_name.upper()} | {family_title} | {family}")
        ax.grid(True, axis="y", alpha=0.28)
        if y_limits is not None:
            ax.set_ylim(*y_limits)
            y_offset = 0.02 * max(1e-12, y_limits[1] - y_limits[0])
        else:
            y_span = max(values) - min(values) if values else 0.0
            y_offset = 0.03 * max(1.0, max(values, default=0.0), y_span)

        for x, value in zip(xpos, values):
            ax.text(x, value + y_offset, fmt.format(value), ha="center", va="bottom", fontsize=8)

    fig.suptitle(f"HELMVEC1 | {figure_title} | {case_name.upper()}", fontsize=14)
    _save_figure(fig, out_path, dpi)


def _generate_spatial_images(
    case_root: Path,
    case_name: str,
    rows: Sequence[ModeRow],
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
) -> None:
    img_root = case_root / "img"
    magnitude_dir = img_root / "magnitude"
    quiver_dir = img_root / "quiver"
    scalar_dir = img_root / "scalar"
    magnitude_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    scalar_dir.mkdir(parents=True, exist_ok=True)

    for meta in rows:
        if not meta.has_spatial_artifacts:
            continue
        points, triangles = _read_legacy_vtk_connectivity(meta.vtk_path)
        if meta.field_data_kind == "edge_vector_cell":
            field_data = _read_edge_field_data(meta.fields_csv_path)
            _plot_edge_magnitude(
                meta,
                points,
                triangles,
                field_data,
                magnitude_dir / f"{meta.spatial_stem()}_magnitude.png",
                dpi,
                show_mesh,
            )
            _plot_edge_quiver(
                meta,
                points,
                triangles,
                field_data,
                quiver_dir / f"{meta.spatial_stem()}_quiver.png",
                dpi,
                show_mesh,
                max_arrows,
            )
        elif meta.field_data_kind == "scalar_nodal":
            field_data = _read_scalar_field_data(meta.fields_csv_path)
            _plot_scalar_mode(
                meta,
                points,
                triangles,
                field_data,
                scalar_dir / f"{meta.spatial_stem()}_scalar.png",
                dpi,
                show_mesh,
            )


def _plot_case_images(
    case_root: Path,
    case_name: str,
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
) -> None:
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
    _plot_modal_summary(
        case_name=case_name,
        rows=rows,
        out_path=img_root / f"helmvec1_{case_name}_error_by_mode.png",
        dpi=dpi,
        value_getter=lambda item: abs(item.error_percent) if item.error_percent is not None else None,
        y_label="|error_percent| [%]",
        family_title="erro percentual absoluto por modo",
        figure_title="Erro Percentual por Modo",
        empty_warning="erro percentual para o resumo",
    )
    _generate_spatial_images(
        case_root=case_root,
        case_name=case_name,
        rows=rows,
        dpi=dpi,
        show_mesh=show_mesh,
        max_arrows=max_arrows,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do HELMVEC1 e gera imagens-resumo e imagens espaciais "
            "didaticas do sistema misto da Eq. (92). Modos edge geram "
            "magnitude/quiver do campo transversal; modos scalar geram mapa "
            "da componente longitudinal."
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
    parser.add_argument("--show-mesh", action="store_true", help="Sobrepoe a malha triangular.")
    parser.add_argument(
        "--max-arrows",
        type=int,
        default=300,
        help="Numero maximo aproximado de setas por figura de quiver. Padrao: 300",
    )
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
        _plot_case_images(
            case_root=case_root,
            case_name=case_name,
            dpi=args.dpi,
            show_mesh=args.show_mesh,
            max_arrows=args.max_arrows,
        )


if __name__ == "__main__":
    main()
