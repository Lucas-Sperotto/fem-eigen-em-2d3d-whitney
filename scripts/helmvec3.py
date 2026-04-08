#!/usr/bin/env python3
"""
Generate HELMVEC3 images directly from the CSV outputs.

The HELMVEC3 family exports:

- `out/helmvec3/rect/csv/helmvec3_rect_table9.csv`
- `out/helmvec3/rect/csv/helmvec3_rect_preview.csv`
- `out/helmvec3/rect/csv/helmvec3_rect_table10.csv`
- one Et CSV + VTK per exported point
- one Ez CSV + VTK per exported point

The script produces:

- Figure 12 / Table 9 summary plots
- preview branch plot
- Table 10 branch and error plots
- magnitude + quiver of Et for each exported point
- scalar map of Ez for each exported point
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


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


def _label_number(value: float) -> str:
    return f"{value:g}".replace(".", "_").replace("-", "m")


@dataclass
class TimingRow:
    row: Dict[str, str]

    @property
    def backend(self) -> str:
        return self.row.get("backend", "")

    @property
    def nx(self) -> Optional[int]:
        value = _to_float(self.row.get("nx"))
        return int(value) if value is not None else None

    @property
    def ny(self) -> Optional[int]:
        value = _to_float(self.row.get("ny"))
        return int(value) if value is not None else None

    @property
    def d13_preview_over_a(self) -> Optional[float]:
        return _to_float(self.row.get("d13_preview_over_a"))

    @property
    def a(self) -> Optional[float]:
        return _to_float(self.row.get("a"))

    @property
    def b(self) -> Optional[float]:
        return _to_float(self.row.get("b"))

    @property
    def d12(self) -> Optional[float]:
        return _to_float(self.row.get("d12"))


@dataclass
class Table9Row:
    case_root: Path
    row: Dict[str, str]

    @property
    def br_over_lambda0(self) -> float:
        return float(self.row["br_over_lambda0"])

    @property
    def beta_over_k0_fem(self) -> float:
        return float(self.row["beta_over_k0_fem"])

    @property
    def beta_over_k0_analytic(self) -> float:
        return float(self.row["beta_over_k0_analytic"])

    @property
    def beta_over_k0_helmvec3(self) -> float:
        return float(self.row["beta_over_k0_helmvec3"])

    @property
    def selected_candidate_rank(self) -> int:
        return int(self.row.get("selected_candidate_rank", "0"))

    @property
    def selected_eig_index(self) -> int:
        return int(self.row.get("selected_eig_index", "-1"))

    @property
    def ez_ratio(self) -> Optional[float]:
        return _to_float(self.row.get("ez_ratio"))

    @property
    def error_percent_analytic(self) -> float:
        return float(self.row["error_percent_analytic"])

    @property
    def error_percent_helmvec3(self) -> float:
        return float(self.row["error_percent_helmvec3"])

    @property
    def match_status(self) -> str:
        return self.row.get("match_status", "")

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "")

    @property
    def et_fields_csv_file(self) -> str:
        return self.row.get("et_fields_csv_file", "")

    @property
    def ez_fields_csv_file(self) -> str:
        return self.row.get("ez_fields_csv_file", "")

    @property
    def et_vtk_file(self) -> str:
        return self.row.get("et_vtk_file", "")

    @property
    def ez_vtk_file(self) -> str:
        return self.row.get("ez_vtk_file", "")

    @property
    def et_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.et_fields_csv_file

    @property
    def ez_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.ez_fields_csv_file

    @property
    def et_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.et_vtk_file

    @property
    def ez_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.ez_vtk_file

    @property
    def has_spatial_artifacts(self) -> bool:
        return (
            self.match_status == "matched"
            and self.et_fields_csv_file
            and self.ez_fields_csv_file
            and self.et_vtk_file
            and self.ez_vtk_file
            and self.et_fields_csv_path.exists()
            and self.ez_fields_csv_path.exists()
            and self.et_vtk_path.exists()
            and self.ez_vtk_path.exists()
        )

    def spatial_stem(self) -> str:
        if self.et_fields_csv_file.endswith("_Et_fields.csv"):
            return self.et_fields_csv_file[: -len("_Et_fields.csv")]
        return f"helmvec3_rect_table9_br{self.br_over_lambda0:g}"

    def spatial_group(self, timing: Optional[TimingRow]) -> str:
        if timing is None or timing.a is None or timing.d12 is None:
            return "figure12"
        return f"figure12_da_{_label_number(timing.d12 / timing.a)}"

    def title_line(self) -> str:
        return rf"RECT | Figura 12 | $b_r/\lambda_0 = {self.br_over_lambda0:g}$"

    def value_line(self) -> str:
        if self.ez_ratio is None:
            return rf"$\beta/k_0 = {self.beta_over_k0_fem:.6f}$"
        return rf"$\beta/k_0 = {self.beta_over_k0_fem:.6f}$ | Ez-ratio={self.ez_ratio:.4f}"

    def interface_spec(self, timing: Optional[TimingRow]) -> Tuple[str, Optional[float]]:
        if timing is None:
            return ("horizontal", None)
        return ("horizontal", timing.d12)


@dataclass
class PreviewRow:
    case_root: Path
    row: Dict[str, str]

    @property
    def d_over_a_preview(self) -> float:
        return float(self.row["d_over_a_preview"])

    @property
    def br_over_lambda0(self) -> float:
        return float(self.row["br_over_lambda0"])

    @property
    def beta_over_k0_fem_branch(self) -> float:
        return float(self.row["beta_over_k0_fem_branch"])

    @property
    def selected_candidate_rank(self) -> int:
        return int(self.row.get("selected_candidate_rank", "0"))

    @property
    def selected_eig_index(self) -> int:
        return int(self.row.get("selected_eig_index", "-1"))

    @property
    def ez_ratio(self) -> Optional[float]:
        return _to_float(self.row.get("ez_ratio"))

    @property
    def branch_status(self) -> str:
        return self.row.get("branch_status", "")

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "")

    @property
    def et_fields_csv_file(self) -> str:
        return self.row.get("et_fields_csv_file", "")

    @property
    def ez_fields_csv_file(self) -> str:
        return self.row.get("ez_fields_csv_file", "")

    @property
    def et_vtk_file(self) -> str:
        return self.row.get("et_vtk_file", "")

    @property
    def ez_vtk_file(self) -> str:
        return self.row.get("ez_vtk_file", "")

    @property
    def et_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.et_fields_csv_file

    @property
    def ez_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.ez_fields_csv_file

    @property
    def et_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.et_vtk_file

    @property
    def ez_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.ez_vtk_file

    @property
    def has_spatial_artifacts(self) -> bool:
        return (
            self.branch_status == "tracked_branch"
            and self.et_fields_csv_file
            and self.ez_fields_csv_file
            and self.et_vtk_file
            and self.ez_vtk_file
            and self.et_fields_csv_path.exists()
            and self.ez_fields_csv_path.exists()
            and self.et_vtk_path.exists()
            and self.ez_vtk_path.exists()
        )

    def spatial_stem(self) -> str:
        if self.et_fields_csv_file.endswith("_Et_fields.csv"):
            return self.et_fields_csv_file[: -len("_Et_fields.csv")]
        return f"helmvec3_rect_preview_da{self.d_over_a_preview:g}_br{self.br_over_lambda0:g}"

    def spatial_group(self, timing: Optional[TimingRow]) -> str:
        return f"preview_da_{_label_number(self.d_over_a_preview)}"

    def title_line(self) -> str:
        return rf"RECT | Preview | $d/a = {self.d_over_a_preview:g}$ | $b_r/\lambda_0 = {self.br_over_lambda0:g}$"

    def value_line(self) -> str:
        if self.ez_ratio is None:
            return rf"$\beta/k_0 = {self.beta_over_k0_fem_branch:.6f}$"
        return rf"$\beta/k_0 = {self.beta_over_k0_fem_branch:.6f}$ | Ez-ratio={self.ez_ratio:.4f}"

    def interface_spec(self, timing: Optional[TimingRow]) -> Tuple[str, Optional[float]]:
        if timing is None or timing.a is None:
            return ("vertical", None)
        return ("vertical", self.d_over_a_preview * timing.a)


@dataclass
class Table10Row:
    case_root: Path
    row: Dict[str, str]

    @property
    def d_over_a(self) -> float:
        return float(self.row["d_over_a"])

    @property
    def br_over_lambda0(self) -> float:
        return float(self.row["br_over_lambda0"])

    @property
    def beta_over_k0_fem_matched(self) -> float:
        return float(self.row["beta_over_k0_fem_matched"])

    @property
    def beta_over_k0_analytic(self) -> float:
        return float(self.row["beta_over_k0_analytic"])

    @property
    def beta_over_k0_helmvec3(self) -> float:
        return float(self.row["beta_over_k0_helmvec3"])

    @property
    def selected_candidate_rank(self) -> int:
        return int(self.row.get("selected_candidate_rank", "0"))

    @property
    def selected_eig_index(self) -> int:
        return int(self.row.get("selected_eig_index", "-1"))

    @property
    def ez_ratio(self) -> Optional[float]:
        return _to_float(self.row.get("ez_ratio"))

    @property
    def error_percent_analytic(self) -> float:
        return float(self.row["error_percent_analytic"])

    @property
    def error_percent_helmvec3(self) -> float:
        return float(self.row["error_percent_helmvec3"])

    @property
    def match_status(self) -> str:
        return self.row.get("match_status", "")

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "")

    @property
    def et_fields_csv_file(self) -> str:
        return self.row.get("et_fields_csv_file", "")

    @property
    def ez_fields_csv_file(self) -> str:
        return self.row.get("ez_fields_csv_file", "")

    @property
    def et_vtk_file(self) -> str:
        return self.row.get("et_vtk_file", "")

    @property
    def ez_vtk_file(self) -> str:
        return self.row.get("ez_vtk_file", "")

    @property
    def et_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.et_fields_csv_file

    @property
    def ez_fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.ez_fields_csv_file

    @property
    def et_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.et_vtk_file

    @property
    def ez_vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.ez_vtk_file

    @property
    def has_spatial_artifacts(self) -> bool:
        return (
            self.match_status == "matched"
            and self.et_fields_csv_file
            and self.ez_fields_csv_file
            and self.et_vtk_file
            and self.ez_vtk_file
            and self.et_fields_csv_path.exists()
            and self.ez_fields_csv_path.exists()
            and self.et_vtk_path.exists()
            and self.ez_vtk_path.exists()
        )

    def spatial_stem(self) -> str:
        if self.et_fields_csv_file.endswith("_Et_fields.csv"):
            return self.et_fields_csv_file[: -len("_Et_fields.csv")]
        return f"helmvec3_rect_table10_da{self.d_over_a:g}_br{self.br_over_lambda0:g}"

    def spatial_group(self, timing: Optional[TimingRow]) -> str:
        return f"table10_da_{_label_number(self.d_over_a)}"

    def title_line(self) -> str:
        return rf"RECT | Tabela 10 | $d/a = {self.d_over_a:g}$ | $b_r/\lambda_0 = {self.br_over_lambda0:g}$"

    def value_line(self) -> str:
        if self.ez_ratio is None:
            return rf"$\beta/k_0 = {self.beta_over_k0_fem_matched:.6f}$"
        return rf"$\beta/k_0 = {self.beta_over_k0_fem_matched:.6f}$ | Ez-ratio={self.ez_ratio:.4f}"

    def interface_spec(self, timing: Optional[TimingRow]) -> Tuple[str, Optional[float]]:
        if timing is None or timing.a is None:
            return ("vertical", None)
        return ("vertical", self.d_over_a * timing.a)


@dataclass
class EdgeFieldData:
    cell_ids: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    ex: np.ndarray
    ey: np.ndarray
    emag: np.ndarray


@dataclass
class ScalarFieldData:
    node_ids: np.ndarray
    x: np.ndarray
    y: np.ndarray
    ez: np.ndarray


def _save_figure(fig: plt.Figure, out_path: Path, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.03, 0.03, 0.98, 0.95))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _load_timing(case_root: Path) -> Optional[TimingRow]:
    timing_path = case_root / "run_timing.csv"
    if not timing_path.exists():
        return None
    rows = _read_rows(timing_path)
    if not rows:
        return None
    return TimingRow(rows[0])


def _title_suffix(timing: Optional[TimingRow]) -> str:
    if timing is None:
        return "RECT"
    parts = ["RECT"]
    if timing.nx is not None and timing.ny is not None:
        parts.append(f"mesh={timing.nx}x{timing.ny}")
    if timing.backend:
        parts.append(f"backend={timing.backend}")
    return " | ".join(parts)


def _plot_table9(rows: List[Table9Row], img_root: Path, dpi: int, title_suffix: str) -> None:
    x = np.array([item.br_over_lambda0 for item in rows], dtype=float)
    fem = [item.beta_over_k0_fem for item in rows]
    ana = [item.beta_over_k0_analytic for item in rows]
    ref = [item.beta_over_k0_helmvec3 for item in rows]

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    ax.plot(x, fem, color="#1f77b4", marker="o", lw=1.9, ms=5.5, label="FEM")
    ax.plot(x, ana, color="#ff7f0e", marker="s", lw=1.6, ms=5.0, ls="--", label="Analitico (ref)")
    ax.plot(x, ref, color="#2ca02c", marker="d", lw=1.6, ms=4.8, ls="-.", label="HELMVEC3 (ref)")
    ax.set_xlabel(r"$b_r / \lambda_0$")
    ax.set_ylabel(r"$\beta / k_0$")
    ax.set_title(f"HELMVEC3 | Figura 12 / Tabela 9 | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _save_figure(fig, img_root / "helmvec3_rect_table9_beta_over_k0.png", dpi)


def _plot_table9_error(rows: List[Table9Row], img_root: Path, dpi: int, title_suffix: str) -> None:
    x = np.array([item.br_over_lambda0 for item in rows], dtype=float)
    err_ana = [abs(item.error_percent_analytic) for item in rows]
    err_ref = [abs(item.error_percent_helmvec3) for item in rows]

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    ax.plot(x, err_ana, color="#ff7f0e", marker="s", lw=1.8, ms=5.2, label="|err| vs Analitico (%)")
    ax.plot(x, err_ref, color="#2ca02c", marker="d", lw=1.8, ms=5.0, label="|err| vs HELMVEC3 (%)")
    ax.set_xlabel(r"$b_r / \lambda_0$")
    ax.set_ylabel("Erro relativo absoluto [%]")
    ax.set_title(f"HELMVEC3 | Erro da Tabela 9 | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _save_figure(fig, img_root / "helmvec3_rect_table9_error_by_point.png", dpi)


def _plot_preview(
    rows: List[PreviewRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
    preview_da: Optional[float],
) -> None:
    x = np.array([item.br_over_lambda0 for item in rows], dtype=float)
    y = [item.beta_over_k0_fem_branch for item in rows]
    da_label = f"d/a={preview_da:g}" if preview_da is not None else "preview"

    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    ax.plot(x, y, color="#1f77b4", marker="o", lw=1.9, ms=5.5)
    ax.set_xlabel(r"$b_r / \lambda_0$")
    ax.set_ylabel(r"$\beta / k_0$")
    ax.set_title(f"HELMVEC3 | Preview de Ramo | {title_suffix} | {da_label}")
    ax.grid(True, alpha=0.28)
    _save_figure(fig, img_root / "helmvec3_rect_preview_branch.png", dpi)


def _plot_table10_branches(rows: List[Table10Row], img_root: Path, dpi: int, title_suffix: str) -> None:
    grouped: Dict[float, List[Table10Row]] = {}
    for row in rows:
        grouped.setdefault(row.d_over_a, []).append(row)

    fig, ax = plt.subplots(figsize=(10.2, 6.0))
    cmap = plt.get_cmap("tab10")
    for idx, d_over_a in enumerate(sorted(grouped.keys())):
        pts = sorted(grouped[d_over_a], key=lambda item: item.br_over_lambda0)
        x = [item.br_over_lambda0 for item in pts]
        y = [item.beta_over_k0_fem_matched for item in pts]
        ax.plot(
            x,
            y,
            marker="o",
            lw=1.6,
            ms=4.4,
            color=cmap(idx % 10),
            label=f"d/a={d_over_a:g}",
        )
    ax.set_xlabel(r"$b_r / \lambda_0$")
    ax.set_ylabel(r"$\beta / k_0$")
    ax.set_title(f"HELMVEC3 | Figura 13 / Tabela 10 | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend(ncol=2, fontsize=8)
    _save_figure(fig, img_root / "helmvec3_rect_table10_fem_branches.png", dpi)


def _plot_table10_error(rows: List[Table10Row], img_root: Path, dpi: int, title_suffix: str) -> None:
    grouped: Dict[float, List[Table10Row]] = {}
    for row in rows:
        grouped.setdefault(row.d_over_a, []).append(row)

    fig, ax = plt.subplots(figsize=(10.2, 6.0))
    cmap = plt.get_cmap("tab10")
    for idx, d_over_a in enumerate(sorted(grouped.keys())):
        pts = sorted(grouped[d_over_a], key=lambda item: item.br_over_lambda0)
        x = [item.br_over_lambda0 for item in pts]
        y = [abs(item.error_percent_analytic) for item in pts]
        ax.plot(
            x,
            y,
            marker="o",
            lw=1.6,
            ms=4.4,
            color=cmap(idx % 10),
            label=f"d/a={d_over_a:g}",
        )
    ax.set_xlabel(r"$b_r / \lambda_0$")
    ax.set_ylabel("|err| vs Analitico [%]")
    ax.set_title(f"HELMVEC3 | Erro da Tabela 10 | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend(ncol=2, fontsize=8)
    _save_figure(fig, img_root / "helmvec3_rect_table10_error_by_branch.png", dpi)


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


def _read_edge_field_data(csv_path: Path) -> EdgeFieldData:
    rows = _read_rows(csv_path)
    rows.sort(key=lambda row: int(row["cell_id"]))
    return EdgeFieldData(
        cell_ids=np.array([int(row["cell_id"]) for row in rows], dtype=int),
        xc=np.array([float(row["xc_m"]) for row in rows], dtype=float),
        yc=np.array([float(row["yc_m"]) for row in rows], dtype=float),
        ex=np.array([float(row["Ex"]) for row in rows], dtype=float),
        ey=np.array([float(row["Ey"]) for row in rows], dtype=float),
        emag=np.array([float(row["Emag"]) for row in rows], dtype=float),
    )


def _read_scalar_field_data(csv_path: Path) -> ScalarFieldData:
    rows = _read_rows(csv_path)
    rows.sort(key=lambda row: int(row["node_id"]))
    return ScalarFieldData(
        node_ids=np.array([int(row["node_id"]) for row in rows], dtype=int),
        x=np.array([float(row["x_m"]) for row in rows], dtype=float),
        y=np.array([float(row["y_m"]) for row in rows], dtype=float),
        ez=np.array([float(row["Ez"]) for row in rows], dtype=float),
    )


def _set_axes_style(ax: plt.Axes, points: np.ndarray) -> None:
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_xlim(points[:, 0].min(), points[:, 0].max())
    ax.set_ylim(points[:, 1].min(), points[:, 1].max())


def _boundary_edges(triangles: np.ndarray) -> List[Tuple[int, int]]:
    edge_counts: Dict[Tuple[int, int], int] = {}
    for tri in triangles:
        i0, i1, i2 = int(tri[0]), int(tri[1]), int(tri[2])
        for a, b in ((i0, i1), (i1, i2), (i2, i0)):
            edge = (a, b) if a < b else (b, a)
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    return [edge for edge, count in edge_counts.items() if count == 1]


def _plot_domain_outline(
    ax: plt.Axes,
    points: np.ndarray,
    triangles: np.ndarray,
    *,
    color: str = "k",
    linewidth: float = 1.0,
    alpha: float = 0.9,
    zorder: int = 4,
) -> None:
    for i, j in _boundary_edges(triangles):
        ax.plot(
            [points[i, 0], points[j, 0]],
            [points[i, 1], points[j, 1]],
            color=color,
            lw=linewidth,
            alpha=alpha,
            zorder=zorder,
        )


def _plot_material_interface(
    ax: plt.Axes,
    points: np.ndarray,
    interface_kind: str,
    interface_pos: Optional[float],
    *,
    color: str = "#8b0000",
    linewidth: float = 1.4,
    alpha: float = 0.9,
    linestyle: str = "--",
    zorder: int = 5,
) -> None:
    if interface_pos is None:
        return
    xmin = float(points[:, 0].min())
    xmax = float(points[:, 0].max())
    ymin = float(points[:, 1].min())
    ymax = float(points[:, 1].max())

    if interface_kind == "horizontal":
        ax.plot(
            [xmin, xmax],
            [interface_pos, interface_pos],
            color=color,
            lw=linewidth,
            alpha=alpha,
            ls=linestyle,
            zorder=zorder,
        )
    elif interface_kind == "vertical":
        ax.plot(
            [interface_pos, interface_pos],
            [ymin, ymax],
            color=color,
            lw=linewidth,
            alpha=alpha,
            ls=linestyle,
            zorder=zorder,
        )


def _thin_indices(n: int, max_arrows: int) -> np.ndarray:
    if n <= max_arrows:
        return np.arange(n, dtype=int)
    step = max(1, n // max_arrows)
    return np.arange(0, n, step, dtype=int)


def _plot_et_magnitude(
    title_line: str,
    value_line: str,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    interface_kind: str,
    interface_pos: Optional[float],
    out_path: Path,
    dpi: int,
    show_mesh: bool,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    tpc = ax.tripcolor(tri, facecolors=field_data.emag, shading="flat", cmap="viridis")
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="white", alpha=0.35, zorder=2)
    _plot_domain_outline(ax, points, triangles, color="white", linewidth=1.1, alpha=0.95, zorder=3)
    _plot_material_interface(ax, points, interface_kind, interface_pos, color="#fff3b0", linewidth=1.5, alpha=0.95, linestyle="--", zorder=4)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label("|E|")
    _set_axes_style(ax, points)
    ax.set_title("HELMVEC3 | Magnitude do Campo Transversal\n" f"{title_line}\n" f"{value_line}")
    _save_figure(fig, out_path, dpi)


def _plot_et_quiver(
    title_line: str,
    value_line: str,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    interface_kind: str,
    interface_pos: Optional[float],
    out_path: Path,
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="0.82", zorder=0)
    _plot_domain_outline(ax, points, triangles, color="0.15", linewidth=1.1, alpha=0.95, zorder=2)
    _plot_material_interface(ax, points, interface_kind, interface_pos, color="#b22222", linewidth=1.5, alpha=0.9, linestyle="--", zorder=3)

    keep = _thin_indices(len(field_data.cell_ids), max_arrows=max_arrows)
    ax.quiver(
        field_data.xc[keep],
        field_data.yc[keep],
        field_data.ex[keep],
        field_data.ey[keep],
        color="#154c79",
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0032,
        zorder=4,
    )
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC3 | Distribuicao do Campo Transversal\n"
        f"{title_line}\n"
        f"{value_line} | setas=(Ex, Ey)"
    )
    _save_figure(fig, out_path, dpi)


def _plot_ez_scalar(
    title_line: str,
    value_line: str,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: ScalarFieldData,
    interface_kind: str,
    interface_pos: Optional[float],
    out_path: Path,
    dpi: int,
    show_mesh: bool,
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    tpc = ax.tripcolor(tri, field_data.ez, shading="gouraud", cmap="coolwarm")
    ax.tricontour(tri, field_data.ez, levels=14, colors="k", linewidths=0.45, alpha=0.55)
    if show_mesh:
        ax.triplot(tri, lw=0.25, color="0.45", alpha=0.25, zorder=3)
    _plot_domain_outline(ax, points, triangles, color="0.1", linewidth=1.0, alpha=0.95, zorder=4)
    _plot_material_interface(ax, points, interface_kind, interface_pos, color="#4b0000", linewidth=1.4, alpha=0.9, linestyle="--", zorder=5)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label("Ez")
    _set_axes_style(ax, points)
    ax.set_title("HELMVEC3 | Funcao Longitudinal Escalar\n" f"{title_line}\n" f"{value_line} | Ez")
    _save_figure(fig, out_path, dpi)


def _generate_spatial_images_for_rows(
    rows: List[object],
    timing: Optional[TimingRow],
    img_root: Path,
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
) -> None:
    magnitude_dir = img_root / "magnitude"
    quiver_dir = img_root / "quiver"
    scalar_dir = img_root / "scalar"
    magnitude_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    scalar_dir.mkdir(parents=True, exist_ok=True)

    for meta in rows:
        if not meta.has_spatial_artifacts:
            continue
        points_et, triangles_et = _read_legacy_vtk_connectivity(meta.et_vtk_path)
        points_ez, triangles_ez = _read_legacy_vtk_connectivity(meta.ez_vtk_path)
        edge_data = _read_edge_field_data(meta.et_fields_csv_path)
        scalar_data = _read_scalar_field_data(meta.ez_fields_csv_path)
        interface_kind, interface_pos = meta.interface_spec(timing)
        group_dir = meta.spatial_group(timing)

        _plot_et_magnitude(
            meta.title_line(),
            meta.value_line(),
            points_et,
            triangles_et,
            edge_data,
            interface_kind,
            interface_pos,
            magnitude_dir / group_dir / f"{meta.spatial_stem()}_Et_magnitude.png",
            dpi,
            show_mesh,
        )
        _plot_et_quiver(
            meta.title_line(),
            meta.value_line(),
            points_et,
            triangles_et,
            edge_data,
            interface_kind,
            interface_pos,
            quiver_dir / group_dir / f"{meta.spatial_stem()}_Et_quiver.png",
            dpi,
            show_mesh,
            max_arrows,
        )
        _plot_ez_scalar(
            meta.title_line(),
            meta.value_line(),
            points_ez,
            triangles_ez,
            scalar_data,
            interface_kind,
            interface_pos,
            scalar_dir / group_dir / f"{meta.spatial_stem()}_Ez_scalar.png",
            dpi,
            show_mesh,
        )


def _plot_case(case_root: Path, dpi: int, show_mesh: bool, max_arrows: int) -> None:
    print("Processing HELMVEC3 case: rect")
    table9_path = case_root / "csv" / "helmvec3_rect_table9.csv"
    preview_path = case_root / "csv" / "helmvec3_rect_preview.csv"
    table10_path = case_root / "csv" / "helmvec3_rect_table10.csv"
    if not table9_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {table9_path}")
    if not preview_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {preview_path}")
    if not table10_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {table10_path}")

    table9_rows = [Table9Row(case_root, row) for row in _read_rows(table9_path)]
    preview_rows = [PreviewRow(case_root, row) for row in _read_rows(preview_path)]
    table10_rows = [Table10Row(case_root, row) for row in _read_rows(table10_path)]
    table9_rows.sort(key=lambda item: item.br_over_lambda0)
    preview_rows.sort(key=lambda item: item.br_over_lambda0)
    table10_rows.sort(key=lambda item: (item.d_over_a, item.br_over_lambda0))
    timing = _load_timing(case_root)
    title_suffix = _title_suffix(timing)
    preview_da = timing.d13_preview_over_a if timing is not None else None

    img_root = case_root / "img"
    img_root.mkdir(parents=True, exist_ok=True)

    _plot_table9(table9_rows, img_root, dpi, title_suffix)
    _plot_table9_error(table9_rows, img_root, dpi, title_suffix)
    _plot_preview(preview_rows, img_root, dpi, title_suffix, preview_da)
    _plot_table10_branches(table10_rows, img_root, dpi, title_suffix)
    _plot_table10_error(table10_rows, img_root, dpi, title_suffix)

    _generate_spatial_images_for_rows(table9_rows, timing, img_root, dpi, show_mesh, max_arrows)
    _generate_spatial_images_for_rows(preview_rows, timing, img_root, dpi, show_mesh, max_arrows)
    _generate_spatial_images_for_rows(table10_rows, timing, img_root, dpi, show_mesh, max_arrows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do HELMVEC3 e gera imagens da Eq. (136): "
            "Tabela 9, preview de ramo, Tabela 10 e visualizacoes espaciais Et/Ez."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("out/helmvec3"),
        help="Diretorio raiz das saidas do HELMVEC3. Padrao: out/helmvec3",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=["rect", "all"],
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
            cases.append("rect")
        else:
            cases.append(item)
    ordered_cases = list(dict.fromkeys(cases))

    for case_name in ordered_cases:
        if case_name != "rect":
            print(f"Warning: caso nao suportado ainda em HELMVEC3: {case_name}")
            continue
        case_root = root / case_name
        if not case_root.exists():
            print(f"Warning: caso ausente, pulando: {case_root}")
            continue
        _plot_case(case_root, args.dpi, args.show_mesh, args.max_arrows)


if __name__ == "__main__":
    main()
