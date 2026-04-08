#!/usr/bin/env python3
"""
Generate HELMVEC2 images directly from the modal CSV outputs.

The HELMVEC2 family exports:

- `out/helmvec2/rect/csv/helmvec2_rect_modes.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_candidates.csv`
- one Et CSV + VTK per matched mode
- one Ez CSV + VTK per matched mode

The script produces:

- Figure 11 / Table 8 matched comparison
- error by mode against HELMVEC2 and Hayata references
- candidate spectrum before/after matching
- magnitude + quiver of Et for each matched mode
- scalar map of Ez for each matched mode
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


@dataclass
class ModeRow:
    case_root: Path
    row: Dict[str, str]

    @property
    def mode(self) -> int:
        return int(self.row["mode"])

    @property
    def matched_candidate_rank(self) -> int:
        return int(self.row["matched_candidate_rank"])

    @property
    def matched_eig_index(self) -> int:
        return int(self.row.get("matched_eig_index", "0"))

    @property
    def k0l_fem_matched(self) -> Optional[float]:
        return _to_float(self.row.get("k0l_fem_matched"))

    @property
    def ref_helmvec2(self) -> Optional[float]:
        return _to_float(self.row.get("ref_helmvec2"))

    @property
    def ref_hayata(self) -> Optional[float]:
        return _to_float(self.row.get("ref_hayata"))

    @property
    def error_percent_helmvec2(self) -> Optional[float]:
        return _to_float(self.row.get("error_percent_helmvec2"))

    @property
    def error_percent_hayata(self) -> Optional[float]:
        return _to_float(self.row.get("error_percent_hayata"))

    @property
    def ez_ratio(self) -> Optional[float]:
        return _to_float(self.row.get("ez_ratio"))

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
            and self.et_fields_csv_path.exists()
            and self.ez_fields_csv_path.exists()
            and self.et_vtk_path.exists()
            and self.ez_vtk_path.exists()
        )

    def spatial_stem(self) -> str:
        return f"helmvec2_rect_mode{self.mode:02d}_cand{self.matched_candidate_rank:02d}"

    def title_line(self) -> str:
        if self.ez_ratio is None:
            return f"RECT | modo {self.mode:02d}"
        return f"RECT | modo {self.mode:02d} | Ez-ratio={self.ez_ratio:.4f}"

    def k0_line(self) -> str:
        if self.k0l_fem_matched is None:
            return "$k_0 L = ?$"
        return f"$k_0 L = {self.k0l_fem_matched:.6f}$"


@dataclass
class CandidateRow:
    row: Dict[str, str]

    @property
    def candidate_rank(self) -> int:
        return int(self.row["candidate_rank"])

    @property
    def eig_index(self) -> int:
        return int(self.row["eig_index"])

    @property
    def k0l(self) -> Optional[float]:
        return _to_float(self.row.get("k0l"))

    @property
    def passes_physical_filter(self) -> bool:
        return int(self.row["passes_physical_filter"]) != 0


@dataclass
class TimingRow:
    row: Dict[str, str]

    @property
    def betaL(self) -> Optional[float]:
        return _to_float(self.row.get("betaL"))

    @property
    def nx(self) -> Optional[int]:
        value = _to_float(self.row.get("nx"))
        return int(value) if value is not None else None

    @property
    def ny(self) -> Optional[int]:
        value = _to_float(self.row.get("ny"))
        return int(value) if value is not None else None

    @property
    def backend(self) -> str:
        return self.row.get("backend", "")

    @property
    def L(self) -> Optional[float]:
        return _to_float(self.row.get("L"))

    @property
    def eps_top(self) -> Optional[float]:
        return _to_float(self.row.get("eps_top"))

    @property
    def eps_bottom(self) -> Optional[float]:
        return _to_float(self.row.get("eps_bottom"))


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


def _annotate_points(ax: plt.Axes, xpos: np.ndarray, values: List[float], fmt: str) -> None:
    if not values:
        return
    ymin, ymax = ax.get_ylim()
    span = max(1e-12, ymax - ymin)
    offset = 0.02 * span
    for x, value in zip(xpos, values):
        ax.text(x, value + offset, fmt.format(value), ha="center", va="bottom", fontsize=8)


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
    if timing.betaL is not None:
        parts.append(f"betaL={timing.betaL:g}")
    if timing.nx is not None and timing.ny is not None:
        parts.append(f"mesh={timing.nx}x{timing.ny}")
    if timing.backend:
        parts.append(f"backend={timing.backend}")
    return " | ".join(parts)


def _plot_matched_modes(
    modes: List[ModeRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
) -> None:
    x = np.array([item.mode for item in modes], dtype=float)
    fem = [item.k0l_fem_matched for item in modes]
    ref_helmvec2 = [item.ref_helmvec2 for item in modes]
    ref_hayata = [item.ref_hayata for item in modes]

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.plot(x, fem, color="#1f77b4", marker="o", lw=1.9, ms=5.5, label="FEM matched")
    ax.plot(x, ref_helmvec2, color="#ff7f0e", marker="s", lw=1.6, ms=5.0, ls="--", label="HELMVEC2 (ref)")
    ax.plot(x, ref_hayata, color="#2ca02c", marker="d", lw=1.6, ms=4.8, ls="-.", label="Hayata (ref)")
    ax.set_xlabel("Modo")
    ax.set_ylabel("k0 L")
    ax.set_title(f"HELMVEC2 | Figura 11 / Tabela 8 | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _annotate_points(ax, x, [value for value in fem if value is not None], "{:.4f}")
    _save_figure(fig, img_root / "helmvec2_rect_table8_k0l_by_mode.png", dpi)


def _plot_mode_errors(
    modes: List[ModeRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
) -> None:
    x = np.array([item.mode for item in modes], dtype=float)
    err_helmvec2 = [abs(item.error_percent_helmvec2 or 0.0) for item in modes]
    err_hayata = [abs(item.error_percent_hayata or 0.0) for item in modes]

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.plot(x, err_helmvec2, color="#ff7f0e", marker="s", lw=1.8, ms=5.2, label="|err| vs HELMVEC2 (%)")
    ax.plot(x, err_hayata, color="#2ca02c", marker="d", lw=1.8, ms=5.0, label="|err| vs Hayata (%)")
    ax.set_xlabel("Modo")
    ax.set_ylabel("Erro relativo absoluto [%]")
    ax.set_title(f"HELMVEC2 | Erro por Modo | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _annotate_points(ax, x, err_helmvec2, "{:.3f}")
    _save_figure(fig, img_root / "helmvec2_rect_error_by_mode.png", dpi)


def _plot_candidate_spectrum(
    candidates: List[CandidateRow],
    modes: List[ModeRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
) -> None:
    matched_by_rank = {
        item.matched_candidate_rank: item.mode
        for item in modes
        if item.matched_candidate_rank > 0 and item.match_status == "matched"
    }
    x_all = np.array([item.candidate_rank for item in candidates], dtype=float)
    y_all = np.array([item.k0l for item in candidates], dtype=float)
    phys_mask = np.array([item.passes_physical_filter for item in candidates], dtype=bool)

    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    if np.any(~phys_mask):
        ax.scatter(
            x_all[~phys_mask],
            y_all[~phys_mask],
            color="#bbbbbb",
            s=18,
            alpha=0.65,
            label="reprovado no filtro fisico",
        )
    ax.scatter(
        x_all[phys_mask],
        y_all[phys_mask],
        color="#1f77b4",
        s=20,
        alpha=0.55,
        label="candidato fisico",
    )

    highlighted_x: List[float] = []
    highlighted_y: List[float] = []
    for candidate in candidates:
        mode_rank = matched_by_rank.get(candidate.candidate_rank)
        if mode_rank is None or candidate.k0l is None:
            continue
        highlighted_x.append(float(candidate.candidate_rank))
        highlighted_y.append(float(candidate.k0l))
        ax.text(
            float(candidate.candidate_rank),
            float(candidate.k0l),
            f"M{mode_rank}",
            fontsize=8,
            ha="left",
            va="bottom",
        )

    if highlighted_x:
        ax.scatter(
            highlighted_x,
            highlighted_y,
            color="#d62728",
            s=34,
            zorder=3,
            label="candidato casado na Tabela 8",
        )

    ax.set_xlabel("Rank do candidato")
    ax.set_ylabel("k0 L")
    ax.set_title(f"HELMVEC2 | Espectro de Candidatos | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _save_figure(fig, img_root / "helmvec2_rect_candidates_k0l_by_rank.png", dpi)


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
    timing: Optional[TimingRow],
    *,
    color: str = "#8b0000",
    linewidth: float = 1.4,
    alpha: float = 0.9,
    linestyle: str = "--",
    zorder: int = 5,
) -> None:
    if timing is None:
        return
    if timing.eps_top is None or timing.eps_bottom is None:
        return
    if abs(timing.eps_top - timing.eps_bottom) < 1e-12:
        return

    xmin = float(points[:, 0].min())
    xmax = float(points[:, 0].max())
    ymin = float(points[:, 1].min())
    ymax = float(points[:, 1].max())
    ymid = 0.5 * (ymin + ymax)

    ax.plot(
        [xmin, xmax],
        [ymid, ymid],
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
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
    timing: Optional[TimingRow],
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    tpc = ax.tripcolor(tri, facecolors=field_data.emag, shading="flat", cmap="viridis")
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="white", alpha=0.35, zorder=2)
    _plot_domain_outline(ax, points, triangles, color="white", linewidth=1.1, alpha=0.95, zorder=3)
    _plot_material_interface(ax, points, timing, color="#fff3b0", linewidth=1.5, alpha=0.95, linestyle="--", zorder=4)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label("|E|")
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC2 | Magnitude do Campo Transversal\n"
        f"{meta.title_line()}\n"
        f"{meta.k0_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_et_quiver(
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: EdgeFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
    timing: Optional[TimingRow],
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    if show_mesh:
        ax.triplot(tri, lw=0.28, color="0.82", zorder=0)
    _plot_domain_outline(ax, points, triangles, color="0.15", linewidth=1.1, alpha=0.95, zorder=2)
    _plot_material_interface(ax, points, timing, color="#b22222", linewidth=1.5, alpha=0.9, linestyle="--", zorder=3)

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
        zorder=3,
    )
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC2 | Distribuicao do Campo Transversal\n"
        f"{meta.title_line()}\n"
        f"{meta.k0_line()} | setas=(Ex, Ey)"
    )
    _save_figure(fig, out_path, dpi)


def _plot_ez_scalar(
    meta: ModeRow,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: ScalarFieldData,
    out_path: Path,
    dpi: int,
    show_mesh: bool,
    timing: Optional[TimingRow],
) -> None:
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    fig, ax = plt.subplots(figsize=(9.4, 7.2))
    tpc = ax.tripcolor(tri, field_data.ez, shading="gouraud", cmap="coolwarm")
    ax.tricontour(tri, field_data.ez, levels=14, colors="k", linewidths=0.45, alpha=0.55)
    if show_mesh:
        ax.triplot(tri, lw=0.25, color="0.45", alpha=0.25, zorder=3)
    _plot_domain_outline(ax, points, triangles, color="0.1", linewidth=1.0, alpha=0.95, zorder=4)
    _plot_material_interface(ax, points, timing, color="#4b0000", linewidth=1.4, alpha=0.9, linestyle="--", zorder=5)
    cbar = fig.colorbar(tpc, ax=ax, shrink=0.94)
    cbar.set_label("Ez")
    _set_axes_style(ax, points)
    ax.set_title(
        "HELMVEC2 | Funcao Longitudinal Escalar\n"
        f"{meta.title_line()}\n"
        f"{meta.k0_line()} | Ez"
    )
    _save_figure(fig, out_path, dpi)


def _generate_spatial_images(
    case_root: Path,
    modes: List[ModeRow],
    dpi: int,
    show_mesh: bool,
    max_arrows: int,
    timing: Optional[TimingRow],
) -> None:
    img_root = case_root / "img"
    magnitude_dir = img_root / "magnitude"
    quiver_dir = img_root / "quiver"
    scalar_dir = img_root / "scalar"
    magnitude_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    scalar_dir.mkdir(parents=True, exist_ok=True)

    for meta in modes:
        if not meta.has_spatial_artifacts:
            continue
        edge_points, edge_triangles = _read_legacy_vtk_connectivity(meta.et_vtk_path)
        scalar_points, scalar_triangles = _read_legacy_vtk_connectivity(meta.ez_vtk_path)
        edge_data = _read_edge_field_data(meta.et_fields_csv_path)
        scalar_data = _read_scalar_field_data(meta.ez_fields_csv_path)

        _plot_et_magnitude(
            meta,
            edge_points,
            edge_triangles,
            edge_data,
            magnitude_dir / f"{meta.spatial_stem()}_Et_magnitude.png",
            dpi,
            show_mesh,
            timing,
        )
        _plot_et_quiver(
            meta,
            edge_points,
            edge_triangles,
            edge_data,
            quiver_dir / f"{meta.spatial_stem()}_Et_quiver.png",
            dpi,
            show_mesh,
            max_arrows,
            timing,
        )
        _plot_ez_scalar(
            meta,
            scalar_points,
            scalar_triangles,
            scalar_data,
            scalar_dir / f"{meta.spatial_stem()}_Ez_scalar.png",
            dpi,
            show_mesh,
            timing,
        )


def _plot_case(case_root: Path, dpi: int, show_mesh: bool, max_arrows: int) -> None:
    print("Processing HELMVEC2 case: rect")
    modes_path = case_root / "csv" / "helmvec2_rect_modes.csv"
    candidates_path = case_root / "csv" / "helmvec2_rect_candidates.csv"
    if not modes_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {modes_path}")
    if not candidates_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {candidates_path}")

    modes = [ModeRow(case_root, row) for row in _read_rows(modes_path)]
    candidates = [CandidateRow(row) for row in _read_rows(candidates_path)]
    modes.sort(key=lambda item: item.mode)
    candidates.sort(key=lambda item: item.candidate_rank)
    timing = _load_timing(case_root)
    title_suffix = _title_suffix(timing)

    img_root = case_root / "img"
    img_root.mkdir(parents=True, exist_ok=True)

    _plot_matched_modes(modes, img_root, dpi, title_suffix)
    _plot_mode_errors(modes, img_root, dpi, title_suffix)
    _plot_candidate_spectrum(candidates, modes, img_root, dpi, title_suffix)
    _generate_spatial_images(case_root, modes, dpi, show_mesh, max_arrows, timing)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do HELMVEC2 e gera imagens-resumo da Eq. (119), "
            "alem das visualizacoes espaciais Et/Ez dos modos casados da Tabela 8."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("out/helmvec2"),
        help="Diretorio raiz das saidas do HELMVEC2. Padrao: out/helmvec2",
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
            print(f"Warning: caso nao suportado ainda em HELMVEC2: {case_name}")
            continue
        case_root = root / case_name
        if not case_root.exists():
            print(f"Warning: caso ausente, pulando: {case_root}")
            continue
        _plot_case(case_root, args.dpi, args.show_mesh, args.max_arrows)


if __name__ == "__main__":
    main()
