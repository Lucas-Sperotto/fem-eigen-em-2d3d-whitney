#!/usr/bin/env python3
"""
Generate FEM3D images directly from the 3D modal CSV outputs.

The FEM3D family exports, per solver/case:

- `out/fem3d{0,1}/<case>/csv/<solver>_<case>_modes.csv`
- one `*_E_fields.csv` per matched mode
- one tetrahedral `*.vtk` per matched mode

This script produces:

- `k0` comparison by mode
- error by mode
- projected magnitude views of the reconstructed 3D field
- projected quiver views of the reconstructed 3D field
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


ROOT = Path(__file__).resolve().parents[1]

SOLVERS: Sequence[str] = ("fem3d0", "fem3d1")
CASES: Sequence[str] = ("air", "half", "cyl", "sphere")
PROJECTIONS: Sequence[Tuple[str, int, int, str, str, str]] = (
    ("xy", 0, 1, "x [m]", "y [m]", "(Ex, Ey)"),
    ("xz", 0, 2, "x [m]", "z [m]", "(Ex, Ez)"),
    ("yz", 1, 2, "y [m]", "z [m]", "(Ey, Ez)"),
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


def _read_rows(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


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
    def nz(self) -> Optional[int]:
        value = _to_float(self.row.get("nz"))
        return int(value) if value is not None else None


@dataclass
class ModeRow:
    solver: str
    case_name: str
    case_root: Path
    row: Dict[str, str]

    @property
    def reference_index(self) -> int:
        return int(self.row["reference_index"])

    @property
    def mode_label(self) -> str:
        return self.row["mode_label"]

    @property
    def k0_analytic(self) -> Optional[float]:
        return _to_float(self.row.get("k0_analytic"))

    @property
    def k0_fem(self) -> Optional[float]:
        return _to_float(self.row.get("k0_fem"))

    @property
    def ref_paper(self) -> Optional[float]:
        return _to_float(self.row.get("ref_paper"))

    @property
    def error_percent_analytic(self) -> Optional[float]:
        return _to_float(self.row.get("error_percent_analytic"))

    @property
    def error_percent_ref_paper(self) -> Optional[float]:
        return _to_float(self.row.get("error_percent_ref_paper"))

    @property
    def matched_eig_index(self) -> int:
        return int(self.row.get("matched_eig_index", "-1"))

    @property
    def match_status(self) -> str:
        return self.row.get("match_status", "")

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "")

    @property
    def fields_csv_file(self) -> str:
        return self.row.get("fields_csv_file", "")

    @property
    def vtk_file(self) -> str:
        return self.row.get("vtk_file", "")

    @property
    def fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.fields_csv_file

    @property
    def vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.vtk_file

    @property
    def has_spatial_artifacts(self) -> bool:
        return (
            self.match_status == "matched"
            and self.fields_csv_file
            and self.vtk_file
            and self.fields_csv_path.exists()
            and self.vtk_path.exists()
        )

    def spatial_stem(self) -> str:
        if self.fields_csv_file.endswith("_E_fields.csv"):
            return self.fields_csv_file[: -len("_E_fields.csv")]
        return f"{self.solver}_{self.case_name}_mode{self.reference_index:02d}_{self.mode_label}"

    def title_line(self) -> str:
        return f"{self.case_name.upper()} | {self.mode_label}"

    def value_line(self) -> str:
        if self.k0_fem is None:
            return "$k_0 = ?$"
        return f"$k_0 = {self.k0_fem:.6f}$"


@dataclass
class FieldData:
    tet_ids: np.ndarray
    xc: np.ndarray
    yc: np.ndarray
    zc: np.ndarray
    ex: np.ndarray
    ey: np.ndarray
    ez: np.ndarray
    emag: np.ndarray


def _load_timing(case_root: Path) -> Optional[TimingRow]:
    timing_path = case_root / "run_timing.csv"
    if not timing_path.exists():
        return None
    rows = _read_rows(timing_path)
    if not rows:
        return None
    return TimingRow(rows[0])


def _title_suffix(solver: str, case_name: str, timing: Optional[TimingRow]) -> str:
    parts = [solver.upper(), case_name.upper()]
    if timing is not None and timing.nx is not None and timing.ny is not None and timing.nz is not None:
        parts.append(f"mesh={timing.nx}x{timing.ny}x{timing.nz}")
    if timing is not None and timing.backend:
        parts.append(f"backend={timing.backend}")
    return " | ".join(parts)


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
    fig.subplots_adjust(left=0.05, right=0.96, bottom=0.08, top=0.84, wspace=0.28)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _read_mode_rows(solver: str, case_name: str, case_root: Path) -> List[ModeRow]:
    csv_path = case_root / "csv" / f"{solver}_{case_name}_modes.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Arquivo ausente: {csv_path}")
    rows = [ModeRow(solver=solver, case_name=case_name, case_root=case_root, row=row) for row in _read_rows(csv_path)]
    rows.sort(key=lambda item: item.reference_index)
    return rows


def _read_field_data(csv_path: Path) -> FieldData:
    rows = _read_rows(csv_path)
    rows.sort(key=lambda row: int(row["tet_id"]))
    return FieldData(
        tet_ids=np.array([int(row["tet_id"]) for row in rows], dtype=int),
        xc=np.array([float(row["xc_m"]) for row in rows], dtype=float),
        yc=np.array([float(row["yc_m"]) for row in rows], dtype=float),
        zc=np.array([float(row["zc_m"]) for row in rows], dtype=float),
        ex=np.array([float(row["Ex"]) for row in rows], dtype=float),
        ey=np.array([float(row["Ey"]) for row in rows], dtype=float),
        ez=np.array([float(row["Ez"]) for row in rows], dtype=float),
        emag=np.array([float(row["Emag"]) for row in rows], dtype=float),
    )


def _read_legacy_vtk_tets(vtk_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    lines = vtk_path.read_text(encoding="utf-8").splitlines()

    points_idx = next((i for i, line in enumerate(lines) if line.strip().startswith("POINTS ")), -1)
    if points_idx < 0:
        raise ValueError(f"{vtk_path}: bloco POINTS ausente")
    num_points = int(lines[points_idx].split()[1])
    points = np.array(
        [[float(v) for v in lines[points_idx + 1 + k].split()[:3]] for k in range(num_points)],
        dtype=float,
    )

    cells_idx = next((i for i, line in enumerate(lines) if line.strip().startswith("CELLS ")), -1)
    if cells_idx < 0:
        raise ValueError(f"{vtk_path}: bloco CELLS ausente")
    num_cells = int(lines[cells_idx].split()[1])
    tets: List[List[int]] = []
    for k in range(num_cells):
        items = [int(v) for v in lines[cells_idx + 1 + k].split()]
        if items[0] != 4:
            raise ValueError(f"{vtk_path}: apenas tetraedros sao suportados")
        tets.append(items[1:5])

    return points, np.array(tets, dtype=int)


def _boundary_nodes(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    face_counts: Dict[Tuple[int, int, int], int] = {}
    for tet in tets:
        i0, i1, i2, i3 = [int(v) for v in tet]
        faces = (
            tuple(sorted((i0, i1, i2))),
            tuple(sorted((i0, i1, i3))),
            tuple(sorted((i0, i2, i3))),
            tuple(sorted((i1, i2, i3))),
        )
        for face in faces:
            face_counts[face] = face_counts.get(face, 0) + 1

    node_ids = set()
    for face, count in face_counts.items():
        if count != 1:
            continue
        node_ids.update(face)

    ordered = np.array(sorted(node_ids), dtype=int)
    return points[ordered]


def _boundary_faces(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    face_counts: Dict[Tuple[int, int, int], int] = {}
    face_owner: Dict[Tuple[int, int, int], Tuple[int, int, int]] = {}
    for tet in tets:
        i0, i1, i2, i3 = [int(v) for v in tet]
        faces = (
            (i0, i1, i2),
            (i0, i1, i3),
            (i0, i2, i3),
            (i1, i2, i3),
        )
        for face in faces:
            key = tuple(sorted(face))
            face_counts[key] = face_counts.get(key, 0) + 1
            face_owner[key] = face

    boundary: List[np.ndarray] = []
    for key, count in face_counts.items():
        if count != 1:
            continue
        face = face_owner[key]
        boundary.append(points[np.array(face, dtype=int)])
    if not boundary:
        return np.empty((0, 3, 3), dtype=float)
    return np.array(boundary, dtype=float)


def _unique_2d(points: np.ndarray) -> np.ndarray:
    if len(points) == 0:
        return points
    rounded = np.round(points, decimals=12)
    _, idx = np.unique(rounded, axis=0, return_index=True)
    return points[np.sort(idx)]


def _cross(o: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    return float((a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]))


def _convex_hull(points: np.ndarray) -> np.ndarray:
    pts = _unique_2d(points)
    if len(pts) <= 2:
        return pts

    ordered = np.array(sorted((float(p[0]), float(p[1])) for p in pts), dtype=float)
    lower: List[np.ndarray] = []
    for p in ordered:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0.0:
            lower.pop()
        lower.append(p)

    upper: List[np.ndarray] = []
    for p in reversed(ordered):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0.0:
            upper.pop()
        upper.append(p)

    hull = lower[:-1] + upper[:-1]
    return np.array(hull, dtype=float)


def _project_points(points: np.ndarray, axis0: int, axis1: int) -> np.ndarray:
    return points[:, [axis0, axis1]]


def _plot_outline(ax: plt.Axes, hull: np.ndarray) -> None:
    if len(hull) == 0:
        return
    if len(hull) == 1:
        ax.scatter(hull[:, 0], hull[:, 1], color="0.1", s=12, zorder=3)
        return
    closed = np.vstack([hull, hull[0]])
    ax.plot(closed[:, 0], closed[:, 1], color="0.1", lw=1.2, alpha=0.95, zorder=3)


def _plot_material_interface(ax: plt.Axes, case_name: str, plane: str, boundary_points: np.ndarray) -> None:
    if case_name != "half":
        return
    zmin = float(boundary_points[:, 2].min())
    zmax = float(boundary_points[:, 2].max())
    zmid = 0.5 * (zmin + zmax)

    if plane == "xz":
        xmin = float(boundary_points[:, 0].min())
        xmax = float(boundary_points[:, 0].max())
        ax.plot([xmin, xmax], [zmid, zmid], color="#8b0000", lw=1.3, ls="--", alpha=0.9, zorder=4)
    elif plane == "yz":
        ymin = float(boundary_points[:, 1].min())
        ymax = float(boundary_points[:, 1].max())
        ax.plot([ymin, ymax], [zmid, zmid], color="#8b0000", lw=1.3, ls="--", alpha=0.9, zorder=4)


def _set_axes_style(ax: plt.Axes, coords: np.ndarray, xlabel: str, ylabel: str) -> None:
    ax.set_aspect("equal")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    xmin, xmax = float(coords[:, 0].min()), float(coords[:, 0].max())
    ymin, ymax = float(coords[:, 1].min()), float(coords[:, 1].max())
    dx = max(1e-9, xmax - xmin)
    dy = max(1e-9, ymax - ymin)
    ax.set_xlim(xmin - 0.04 * dx, xmax + 0.04 * dx)
    ax.set_ylim(ymin - 0.04 * dy, ymax + 0.04 * dy)


def _thin_indices(n: int, max_arrows: int) -> np.ndarray:
    if n <= max_arrows:
        return np.arange(n, dtype=int)
    step = max(1, n // max_arrows)
    return np.arange(0, n, step, dtype=int)


def _set_axes_3d_equal(ax: plt.Axes, points: np.ndarray) -> None:
    xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
    ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
    zmin, zmax = float(points[:, 2].min()), float(points[:, 2].max())

    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    zmid = 0.5 * (zmin + zmax)
    radius = 0.55 * max(xmax - xmin, ymax - ymin, zmax - zmin, 1e-9)

    ax.set_xlim(xmid - radius, xmid + radius)
    ax.set_ylim(ymid - radius, ymid + radius)
    ax.set_zlim(zmid - radius, zmid + radius)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")


def _plot_material_interface_3d(ax: plt.Axes, case_name: str, boundary_points: np.ndarray) -> None:
    if case_name != "half":
        return

    xmin, xmax = float(boundary_points[:, 0].min()), float(boundary_points[:, 0].max())
    ymin, ymax = float(boundary_points[:, 1].min()), float(boundary_points[:, 1].max())
    zmin, zmax = float(boundary_points[:, 2].min()), float(boundary_points[:, 2].max())
    zmid = 0.5 * (zmin + zmax)

    xx, yy = np.meshgrid(np.array([xmin, xmax]), np.array([ymin, ymax]))
    zz = np.full_like(xx, zmid, dtype=float)
    ax.plot_surface(xx, yy, zz, color="#8b0000", alpha=0.12, linewidth=0.0, shade=False, zorder=0)
    ax.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], [zmid] * 5, color="#8b0000", lw=1.2, alpha=0.9)


def _plot_geometry_surface_3d(ax: plt.Axes, boundary_faces: np.ndarray) -> None:
    if boundary_faces.size == 0:
        return
    surface = Poly3DCollection(
        boundary_faces,
        facecolor=(0.72, 0.72, 0.72, 0.10),
        edgecolor=(0.35, 0.35, 0.35, 0.45),
        linewidths=0.35,
    )
    ax.add_collection3d(surface)


def _plot_magnitude_3d(
    meta: ModeRow,
    field_data: FieldData,
    boundary_points: np.ndarray,
    boundary_faces: np.ndarray,
    out_path: Path,
    dpi: int,
) -> None:
    coords3 = np.column_stack([field_data.xc, field_data.yc, field_data.zc])
    fig = plt.figure(figsize=(8.8, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    _plot_geometry_surface_3d(ax, boundary_faces)
    ax.scatter(
        boundary_points[:, 0],
        boundary_points[:, 1],
        boundary_points[:, 2],
        color="0.55",
        s=3,
        alpha=0.16,
        depthshade=False,
    )
    scatter = ax.scatter(
        coords3[:, 0],
        coords3[:, 1],
        coords3[:, 2],
        c=field_data.emag,
        cmap="viridis",
        s=16,
        alpha=0.95,
        depthshade=False,
    )
    _plot_material_interface_3d(ax, meta.case_name, boundary_points)
    _set_axes_3d_equal(ax, boundary_points)
    ax.set_title(
        "FEM3D | Distribuicao 3D da Magnitude\n"
        f"{meta.solver.upper()} | {meta.title_line()}\n"
        f"{meta.value_line()}"
    )
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.82, pad=0.08)
    cbar.set_label("|E|")
    _save_figure(fig, out_path, dpi)


def _plot_quiver_3d(
    meta: ModeRow,
    field_data: FieldData,
    boundary_points: np.ndarray,
    boundary_faces: np.ndarray,
    out_path: Path,
    dpi: int,
    max_arrows: int,
) -> None:
    coords3 = np.column_stack([field_data.xc, field_data.yc, field_data.zc])
    keep = _thin_indices(len(field_data.tet_ids), max_arrows=max_arrows)

    fig = plt.figure(figsize=(8.8, 7.0))
    ax = fig.add_subplot(111, projection="3d")
    _plot_geometry_surface_3d(ax, boundary_faces)
    ax.scatter(
        boundary_points[:, 0],
        boundary_points[:, 1],
        boundary_points[:, 2],
        color="0.58",
        s=3,
        alpha=0.14,
        depthshade=False,
    )
    ax.quiver(
        coords3[keep, 0],
        coords3[keep, 1],
        coords3[keep, 2],
        field_data.ex[keep],
        field_data.ey[keep],
        field_data.ez[keep],
        length=0.12 * max(
            float(boundary_points[:, 0].max() - boundary_points[:, 0].min()),
            float(boundary_points[:, 1].max() - boundary_points[:, 1].min()),
            float(boundary_points[:, 2].max() - boundary_points[:, 2].min()),
            1e-9,
        ),
        normalize=False,
        color="#154c79",
        linewidth=0.8,
    )
    _plot_material_interface_3d(ax, meta.case_name, boundary_points)
    _set_axes_3d_equal(ax, boundary_points)
    ax.set_title(
        "FEM3D | Distribuicao 3D do Campo Eletrico\n"
        f"{meta.solver.upper()} | {meta.title_line()}\n"
        f"{meta.value_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_k0_by_mode(
    solver: str,
    case_name: str,
    modes: List[ModeRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
) -> None:
    xpos = np.arange(1, len(modes) + 1, dtype=float)
    labels = [item.mode_label for item in modes]
    fem = [item.k0_fem for item in modes]
    ana = [item.k0_analytic for item in modes]
    ref = [item.ref_paper for item in modes]

    fig, ax = plt.subplots(figsize=(10.6, 5.6))
    ax.plot(xpos, fem, color="#1f77b4", marker="o", lw=1.9, ms=5.4, label="FEM")
    ax.plot(xpos, ana, color="#ff7f0e", marker="s", lw=1.6, ms=5.0, ls="--", label="Analitico")
    ax.plot(xpos, ref, color="#2ca02c", marker="d", lw=1.6, ms=4.8, ls="-.", label="Ref. paper")
    ax.set_xticks(xpos, labels, rotation=35, ha="right")
    ax.set_ylabel("k0")
    ax.set_title(f"{solver.upper()} | k0 por Modo | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _annotate_points(ax, xpos, [value if value is not None else np.nan for value in fem], "{:.3f}")
    _save_figure(fig, img_root / f"{solver}_{case_name}_k0_by_mode.png", dpi)


def _plot_error_by_mode(
    solver: str,
    case_name: str,
    modes: List[ModeRow],
    img_root: Path,
    dpi: int,
    title_suffix: str,
) -> None:
    xpos = np.arange(1, len(modes) + 1, dtype=float)
    labels = [item.mode_label for item in modes]
    err_ana = [abs(item.error_percent_analytic or 0.0) for item in modes]
    err_ref = [abs(item.error_percent_ref_paper or 0.0) for item in modes]

    fig, ax = plt.subplots(figsize=(10.6, 5.6))
    ax.plot(xpos, err_ana, color="#ff7f0e", marker="s", lw=1.8, ms=5.2, label="|err| vs Analitico [%]")
    ax.plot(xpos, err_ref, color="#2ca02c", marker="d", lw=1.8, ms=5.0, label="|err| vs Ref. paper [%]")
    ax.set_xticks(xpos, labels, rotation=35, ha="right")
    ax.set_ylabel("Erro relativo absoluto [%]")
    ax.set_title(f"{solver.upper()} | Erro por Modo | {title_suffix}")
    ax.grid(True, alpha=0.28)
    ax.legend()
    _annotate_points(ax, xpos, err_ana, "{:.3f}")
    _save_figure(fig, img_root / f"{solver}_{case_name}_error_by_mode.png", dpi)


def _plot_magnitude_projections(
    meta: ModeRow,
    field_data: FieldData,
    boundary_points: np.ndarray,
    out_path: Path,
    dpi: int,
) -> None:
    coords3 = np.column_stack([field_data.xc, field_data.yc, field_data.zc])
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.1))
    scatter = None
    for ax, (plane, axis0, axis1, xlabel, ylabel, _) in zip(axes, PROJECTIONS):
        coords2 = _project_points(coords3, axis0, axis1)
        boundary2 = _project_points(boundary_points, axis0, axis1)
        hull = _convex_hull(boundary2)
        scatter = ax.scatter(coords2[:, 0], coords2[:, 1], c=field_data.emag, cmap="viridis", s=18, alpha=0.92)
        _plot_outline(ax, hull)
        _plot_material_interface(ax, meta.case_name, plane, boundary_points)
        _set_axes_style(ax, coords2, xlabel, ylabel)
        ax.set_title(plane.upper())
    cbar = fig.colorbar(scatter, ax=axes, shrink=0.92)
    cbar.set_label("|E|")
    fig.suptitle(
        "FEM3D | Magnitude do Campo Eletrico\n"
        f"{meta.solver.upper()} | {meta.title_line()}\n"
        f"{meta.value_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_quiver_projections(
    meta: ModeRow,
    field_data: FieldData,
    boundary_points: np.ndarray,
    out_path: Path,
    dpi: int,
    max_arrows: int,
) -> None:
    coords3 = np.column_stack([field_data.xc, field_data.yc, field_data.zc])
    comps = np.column_stack([field_data.ex, field_data.ey, field_data.ez])
    keep = _thin_indices(len(field_data.tet_ids), max_arrows=max_arrows)

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.1))
    for ax, (plane, axis0, axis1, xlabel, ylabel, vec_label) in zip(axes, PROJECTIONS):
        coords2 = _project_points(coords3, axis0, axis1)
        boundary2 = _project_points(boundary_points, axis0, axis1)
        hull = _convex_hull(boundary2)
        ax.quiver(
            coords2[keep, 0],
            coords2[keep, 1],
            comps[keep, axis0],
            comps[keep, axis1],
            color="#154c79",
            angles="xy",
            scale_units="xy",
            scale=None,
            width=0.0030,
            zorder=2,
        )
        _plot_outline(ax, hull)
        _plot_material_interface(ax, meta.case_name, plane, boundary_points)
        _set_axes_style(ax, coords2, xlabel, ylabel)
        ax.set_title(f"{plane.upper()} | {vec_label}")
    fig.suptitle(
        "FEM3D | Projecoes do Campo Eletrico\n"
        f"{meta.solver.upper()} | {meta.title_line()}\n"
        f"{meta.value_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _generate_spatial_images(
    solver: str,
    case_name: str,
    case_root: Path,
    modes: List[ModeRow],
    dpi: int,
    max_arrows: int,
) -> None:
    img_root = case_root / "img"
    magnitude_dir = img_root / "magnitude"
    quiver_dir = img_root / "quiver"
    scatter3d_dir = img_root / "3d_scatter"
    quiver3d_dir = img_root / "3d_quiver"
    magnitude_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    scatter3d_dir.mkdir(parents=True, exist_ok=True)
    quiver3d_dir.mkdir(parents=True, exist_ok=True)

    for meta in modes:
        if not meta.has_spatial_artifacts:
            continue
        points, tets = _read_legacy_vtk_tets(meta.vtk_path)
        boundary_points = _boundary_nodes(points, tets)
        boundary_faces = _boundary_faces(points, tets)
        field_data = _read_field_data(meta.fields_csv_path)

        _plot_magnitude_projections(
            meta,
            field_data,
            boundary_points,
            magnitude_dir / f"{meta.spatial_stem()}_magnitude_proj.png",
            dpi,
        )
        _plot_quiver_projections(
            meta,
            field_data,
            boundary_points,
            quiver_dir / f"{meta.spatial_stem()}_quiver_proj.png",
            dpi,
            max_arrows,
        )
        _plot_magnitude_3d(
            meta,
            field_data,
            boundary_points,
            boundary_faces,
            scatter3d_dir / f"{meta.spatial_stem()}_scatter3d.png",
            dpi,
        )
        _plot_quiver_3d(
            meta,
            field_data,
            boundary_points,
            boundary_faces,
            quiver3d_dir / f"{meta.spatial_stem()}_quiver3d.png",
            dpi,
            max_arrows,
        )


def _plot_case(
    solver: str,
    case_name: str,
    case_root: Path,
    dpi: int,
    max_arrows: int,
) -> None:
    print(f"Processing FEM3D case: {solver}/{case_name}")
    modes = _read_mode_rows(solver, case_name, case_root)
    timing = _load_timing(case_root)
    title_suffix = _title_suffix(solver, case_name, timing)

    img_root = case_root / "img"
    img_root.mkdir(parents=True, exist_ok=True)
    _plot_k0_by_mode(solver, case_name, modes, img_root, dpi, title_suffix)
    _plot_error_by_mode(solver, case_name, modes, img_root, dpi, title_suffix)
    _generate_spatial_images(solver, case_name, case_root, modes, dpi, max_arrows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do FEM3D0/FEM3D1 e gera figuras-resumo e "
            "projecoes ortogonais dos campos vetoriais 3D exportados por modo."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("out"),
        help="Diretorio raiz das saidas do projeto. Padrao: out",
    )
    parser.add_argument(
        "--solver",
        action="append",
        choices=["fem3d0", "fem3d1", "all"],
        default=None,
        help="Seleciona quais solvers 3D processar. Pode ser repetido.",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=["air", "half", "cyl", "sphere", "all"],
        default=None,
        help="Seleciona quais casos 3D processar. Pode ser repetido.",
    )
    parser.add_argument("--dpi", type=int, default=180, help="Resolucao das imagens.")
    parser.add_argument(
        "--max-arrows",
        type=int,
        default=350,
        help="Numero maximo aproximado de setas por figura de quiver. Padrao: 350",
    )
    args = parser.parse_args()

    root = _resolve(args.root)
    if not root.exists():
        raise SystemExit(f"Diretorio de saida nao encontrado: {root}")

    solver_args = args.solver or ["all"]
    selected_solvers: List[str] = []
    for item in solver_args:
        if item == "all":
            selected_solvers.extend(SOLVERS)
        else:
            selected_solvers.append(item)
    selected_solvers = list(dict.fromkeys(selected_solvers))

    case_args = args.case or ["all"]
    selected_cases: List[str] = []
    for item in case_args:
        if item == "all":
            selected_cases.extend(CASES)
        else:
            selected_cases.append(item)
    selected_cases = list(dict.fromkeys(selected_cases))

    for solver in selected_solvers:
        for case_name in selected_cases:
            case_root = root / solver / case_name
            if not case_root.exists():
                print(f"Warning: caso ausente, pulando: {case_root}")
                continue
            try:
                _plot_case(solver, case_name, case_root, args.dpi, args.max_arrows)
            except FileNotFoundError as exc:
                print(f"Warning: {exc}")


if __name__ == "__main__":
    main()
