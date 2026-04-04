#!/usr/bin/env python3
"""
Generate HELM10 images directly from the CSV outputs produced by the solvers.

The script scans `out/helm10/<case>/csv`, reads:

- `helm10_<case>_modes.csv`
- `helm10_<case>_fields_<mode>.csv`

and produces, for each exported mode:

- an isopotential plot for the scalar potential `psi`
- a quiver plot for the transverse field `(Ex, Ey)`

It also creates summary images per case using:

- `error_percent`
- `rho_abs`

When a matching VTK file exists, the script reuses its triangle connectivity.
This is especially helpful for the coaxial case, where a naive Delaunay
triangulation of nodal CSV coordinates would incorrectly fill the inner hole.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


@dataclass
class FieldData:
    node_ids: np.ndarray
    x: np.ndarray
    y: np.ndarray
    psi: np.ndarray
    dpsi_dx: np.ndarray
    dpsi_dy: np.ndarray
    ex: np.ndarray
    ey: np.ndarray


@dataclass(frozen=True)
class QuiverSource:
    key: str
    title: str
    file_suffix: str
    u_label: str
    v_label: str


QUIVER_SOURCES: Dict[str, QuiverSource] = {
    "field": QuiverSource(
        key="field",
        title="Distribuicao do Campo Eletrico",
        file_suffix="field",
        u_label="Ex",
        v_label="Ey",
    ),
    "gradient": QuiverSource(
        key="gradient",
        title="Distribuicao do Gradiente do Potencial",
        file_suffix="gradient",
        u_label="dpsi_dx",
        v_label="dpsi_dy",
    ),
}


@dataclass
class ModeMeta:
    case_name: str
    case_root: Path
    row: Dict[str, str]

    @property
    def family(self) -> str:
        return self.row["family"]

    @property
    def longitudinal_label(self) -> str:
        return self.row["longitudinal_label"]

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
    def fields_csv_name(self) -> str:
        return self.row["fields_csv_file"]

    @property
    def fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.fields_csv_name

    @property
    def index2_name(self) -> str:
        return "n" if self.case_name == "rect" else "p"

    @property
    def index2_value(self) -> int:
        return int(self.row[self.index2_name])

    @property
    def m(self) -> int:
        return int(self.row["m"])

    @property
    def error_percent(self) -> float:
        return abs(_to_float(self.row.get("error_percent")) or 0.0)

    @property
    def rho_abs(self) -> Optional[float]:
        return _to_float(self.row.get("rho_abs"))

    @property
    def kc_fem(self) -> Optional[float]:
        return _to_float(self.row.get("kc_fem"))

    @property
    def kc_ana(self) -> Optional[float]:
        return _to_float(self.row.get("kc_ana"))

    @property
    def field_status(self) -> str:
        return self.row.get("field_status", "")

    @property
    def geometry_tuple(self) -> Tuple[float, ...]:
        if self.case_name == "rect":
            return (
                _to_float(self.row["ar_m"]) or 0.0,
                _to_float(self.row["b_m"]) or 0.0,
            )
        if self.case_name == "circle":
            return (_to_float(self.row["r_m"]) or 0.0,)
        return (
            _to_float(self.row["r1_m"]) or 0.0,
            _to_float(self.row["r2_m"]) or 0.0,
        )

    def vtk_path(self) -> Path:
        rank = f"{self.positive_rank:02d}"
        return (
            self.case_root
            / "vtk"
            / (
                f"{self.family.lower()}_{self.case_name}_m{self.m}_"
                f"{self.index2_name}{self.index2_value}_rank{rank}_sv.vtk"
            )
        )

    def short_name(self) -> str:
        return f"{self.mode_label}_rank{self.positive_rank:02d}"

    def base_stem(self) -> str:
        return f"helm10_{self.case_name}_{self.short_name()}"

    def plot_mode_line(self) -> str:
        return f"{self.case_name.upper()} | {self.family}$_{{{self.m},{self.index2_value}}}$"

    def normalized_cutoff_line(self) -> str:
        if self.case_name == "rect":
            value = _to_float(self.row.get("kc_ar_fem"))
            if value is None and self.kc_fem is not None:
                a, _ = self.geometry_tuple
                value = self.kc_fem * a
            return f"$k_c a_r = {value:.6f}$" if value is not None else "$k_c a_r = ?$"

        if self.case_name == "circle":
            value = _to_float(self.row.get("kc_r_fem"))
            if value is None and self.kc_fem is not None:
                (radius,) = self.geometry_tuple
                value = self.kc_fem * radius
            return f"$k_c r = {value:.6f}$" if value is not None else "$k_c r = ?$"

        r1, _ = self.geometry_tuple
        value = self.kc_fem * r1 if self.kc_fem is not None else None
        return f"$k_c r_1 = {value:.6f}$" if value is not None else "$k_c r_1 = ?$"

    def potential_mode_line(self) -> str:
        return f"{self.plot_mode_line()} | {self.longitudinal_label}"


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


def _read_field_data(csv_path: Path) -> FieldData:
    rows = _read_rows(csv_path)
    if not rows:
        raise ValueError(f"{csv_path}: CSV de campos vazio")

    # We sort by node_id to align the arrays with the original mesh numbering.
    rows.sort(key=lambda row: int(row["node_id"]))

    def col(name: str) -> np.ndarray:
        return np.array([float(row[name]) for row in rows], dtype=float)

    return FieldData(
        node_ids=np.array([int(row["node_id"]) for row in rows], dtype=int),
        x=col("x_m"),
        y=col("y_m"),
        psi=col("psi"),
        dpsi_dx=col("dpsi_dx"),
        dpsi_dy=col("dpsi_dy"),
        ex=col("Ex"),
        ey=col("Ey"),
    )


def _read_legacy_vtk_connectivity(vtk_path: Path) -> Tuple[int, np.ndarray]:
    lines = vtk_path.read_text(encoding="utf-8").splitlines()

    points_idx = _find_line_index(lines, "POINTS ")
    if points_idx < 0:
        raise ValueError(f"{vtk_path}: bloco POINTS ausente")
    num_points = int(lines[points_idx].split()[1])

    cells_idx = _find_line_index(lines, "CELLS ")
    if cells_idx < 0:
        raise ValueError(f"{vtk_path}: bloco CELLS ausente")
    num_cells = int(lines[cells_idx].split()[1])

    triangles: List[List[int]] = []
    for k in range(num_cells):
        items = [int(v) for v in lines[cells_idx + 1 + k].split()]
        if items[0] != 3:
            raise ValueError(f"{vtk_path}: apenas triangulos sao suportados")
        triangles.append(items[1:4])

    return num_points, np.array(triangles, dtype=int)


def _find_line_index(lines: Sequence[str], prefix: str) -> int:
    for idx, line in enumerate(lines):
        if line.strip().startswith(prefix):
            return idx
    return -1


def _point_inside_geometry(case_name: str, meta: ModeMeta, x: float, y: float) -> bool:
    scale = max(1.0, *meta.geometry_tuple)
    tol = 1.0e-8 * scale

    if case_name == "rect":
        a, b = meta.geometry_tuple
        return -tol <= x <= a + tol and -tol <= y <= b + tol

    if case_name == "circle":
        (radius,) = meta.geometry_tuple
        return math.hypot(x, y) <= radius + tol

    r1, r2 = meta.geometry_tuple
    radius = math.hypot(x, y)
    return r1 - tol <= radius <= r2 + tol


def _masked_triangulation_from_csv(meta: ModeMeta, field_data: FieldData) -> mtri.Triangulation:
    tri = mtri.Triangulation(field_data.x, field_data.y)
    mask = []
    for triangle in tri.triangles:
        pts = [(field_data.x[idx], field_data.y[idx]) for idx in triangle]

        # We test vertices, edge midpoints and centroid. This avoids filling the
        # circular/coaxial exterior when the fallback Delaunay path is used.
        samples = list(pts)
        for i0, i1 in ((0, 1), (1, 2), (2, 0)):
            samples.append(
                (
                    0.5 * (pts[i0][0] + pts[i1][0]),
                    0.5 * (pts[i0][1] + pts[i1][1]),
                )
            )
        samples.append(
            (
                (pts[0][0] + pts[1][0] + pts[2][0]) / 3.0,
                (pts[0][1] + pts[1][1] + pts[2][1]) / 3.0,
            )
        )
        keep = all(
            _point_inside_geometry(meta.case_name, meta, x, y) for x, y in samples
        )
        mask.append(not keep)

    tri.set_mask(np.array(mask, dtype=bool))
    return tri


def _build_triangulation(meta: ModeMeta, field_data: FieldData, prefer_vtk: bool) -> mtri.Triangulation:
    if prefer_vtk:
        vtk_path = meta.vtk_path()
        if vtk_path.exists():
            num_points, triangles = _read_legacy_vtk_connectivity(vtk_path)
            if num_points == len(field_data.node_ids):
                return mtri.Triangulation(field_data.x, field_data.y, triangles)
            print(
                f"Warning: {vtk_path.name} possui {num_points} pontos, "
                f"mas {meta.fields_csv_name} possui {len(field_data.node_ids)} linhas; "
                "caindo para triangulacao a partir do CSV."
            )

    return _masked_triangulation_from_csv(meta, field_data)


def _draw_geometry_outline(ax: plt.Axes, meta: ModeMeta) -> None:
    if meta.case_name == "rect":
        a, b = meta.geometry_tuple
        xs = [0.0, a, a, 0.0, 0.0]
        ys = [0.0, 0.0, b, b, 0.0]
        ax.plot(xs, ys, color="black", lw=1.2, zorder=4)
        return

    theta = np.linspace(0.0, 2.0 * math.pi, 512)
    if meta.case_name == "circle":
        (radius,) = meta.geometry_tuple
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="black", lw=1.2, zorder=4)
        return

    r1, r2 = meta.geometry_tuple
    ax.plot(r1 * np.cos(theta), r1 * np.sin(theta), color="black", lw=1.2, zorder=4)
    ax.plot(r2 * np.cos(theta), r2 * np.sin(theta), color="black", lw=1.2, zorder=4)


def _set_axes_style(ax: plt.Axes, meta: ModeMeta) -> None:
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(True, alpha=0.18)

    if meta.case_name == "rect":
        a, b = meta.geometry_tuple
        pad = 0.04 * max(a, b, 1.0)
        ax.set_xlim(-pad, a + pad)
        ax.set_ylim(-pad, b + pad)
        return

    if meta.case_name == "circle":
        (radius,) = meta.geometry_tuple
        pad = 0.04 * max(radius, 1.0)
        ax.set_xlim(-radius - pad, radius + pad)
        ax.set_ylim(-radius - pad, radius + pad)
        return

    _, r2 = meta.geometry_tuple
    pad = 0.04 * max(r2, 1.0)
    ax.set_xlim(-r2 - pad, r2 + pad)
    ax.set_ylim(-r2 - pad, r2 + pad)


def _plot_bounds(meta: ModeMeta) -> Tuple[float, float, float, float]:
    if meta.case_name == "rect":
        a, b = meta.geometry_tuple
        return 0.0, a, 0.0, b

    if meta.case_name == "circle":
        (radius,) = meta.geometry_tuple
        return -radius, radius, -radius, radius

    _, r2 = meta.geometry_tuple
    return -r2, r2, -r2, r2


def _build_polar_quiver_points(
    inner_radius: float,
    outer_radius: float,
    num_circles: int = 13,
    points_per_circle: int = 51,
) -> Tuple[np.ndarray, np.ndarray]:
    radii = np.linspace(inner_radius, outer_radius, max(2, num_circles))
    theta = np.linspace(0.0, 2.0 * math.pi, max(3, points_per_circle), endpoint=False)
    rr, tt = np.meshgrid(radii, theta, indexing="ij")
    gx = rr * np.cos(tt)
    gy = rr * np.sin(tt)
    return gx, gy


def _build_quiver_grid(
    meta: ModeMeta,
    tri: mtri.Triangulation,
    field_data: FieldData,
    nx: int,
    ny: int,
    source: QuiverSource,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if meta.case_name == "circle":
        (radius,) = meta.geometry_tuple
        # For the circular guide, a polar sampling follows the geometry better
        # than a Cartesian grid. We use 13 concentric circles with 51 angular
        # points each. The innermost radius is R/13 to avoid collapsing all 51
        # samples into a single point at the center.
        gx, gy = _build_polar_quiver_points(
            inner_radius=radius / 13.0,
            outer_radius=radius,
            num_circles=13,
            points_per_circle=51,
        )
    elif meta.case_name == "coax":
        # For the annular coaxial section, a polar sampling also follows the
        # geometry better than a Cartesian grid. The current didactic choice is
        # 13 concentric circles with 51 angular points each.
        r1, r2 = meta.geometry_tuple
        gx, gy = _build_polar_quiver_points(
            inner_radius=r1,
            outer_radius=r2,
            num_circles=13,
            points_per_circle=51,
        )
    else:
        xmin, xmax, ymin, ymax = _plot_bounds(meta)
        gx, gy = np.meshgrid(
            np.linspace(xmin, xmax, max(2, nx)),
            np.linspace(ymin, ymax, max(2, ny)),
        )

    # We interpolate nodal values onto a regular didactic grid so the quiver
    # pattern stays visually stable across different meshes.
    if source.key == "field":
        u_nodal = field_data.ex
        v_nodal = field_data.ey
    else:
        u_nodal = field_data.dpsi_dx
        v_nodal = field_data.dpsi_dy

    interp_u = mtri.LinearTriInterpolator(tri, u_nodal)
    interp_v = mtri.LinearTriInterpolator(tri, v_nodal)
    u_grid = interp_u(gx, gy)
    v_grid = interp_v(gx, gy)

    inside_geom = np.vectorize(
        lambda x, y: _point_inside_geometry(meta.case_name, meta, float(x), float(y))
    )(gx, gy)

    u_mask = np.ma.getmaskarray(u_grid)
    v_mask = np.ma.getmaskarray(v_grid)
    valid = inside_geom & ~u_mask & ~v_mask

    return (
        gx[valid].astype(float),
        gy[valid].astype(float),
        np.asarray(u_grid[valid], dtype=float),
        np.asarray(v_grid[valid], dtype=float),
    )


def _save_figure(fig: plt.Figure, out_path: Path, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {out_path}")


def _normalize_display_vectors(u: np.ndarray, v: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    magnitude = np.sqrt(u * u + v * v)
    vmax = float(np.max(magnitude)) if magnitude.size else 0.0
    if vmax <= 0.0:
        return u, v
    return u / vmax, v / vmax


def _plot_isopotential(
    meta: ModeMeta,
    field_data: FieldData,
    tri: mtri.Triangulation,
    out_path: Path,
    dpi: int,
    levels: int,
    show_mesh: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.2))

    psi_min = float(np.min(field_data.psi))
    psi_max = float(np.max(field_data.psi))
    if math.isclose(psi_min, psi_max, rel_tol=1.0e-12, abs_tol=1.0e-12):
        ax.tripcolor(tri, field_data.psi, shading="gouraud", cmap="viridis")
    else:
        contourf = ax.tricontourf(tri, field_data.psi, levels=levels, cmap="viridis")
        ax.tricontour(tri, field_data.psi, levels=levels, colors="white", linewidths=0.55, alpha=0.7)
        fig.colorbar(contourf, ax=ax, label="psi")

    if show_mesh:
        ax.triplot(tri, color="0.80", lw=0.25, alpha=0.8, zorder=3)

    _draw_geometry_outline(ax, meta)
    _set_axes_style(ax, meta)
    ax.set_title(
        f"HELM10 | Funcao Potencial Escalar\n"
        f"{meta.potential_mode_line()}\n"
        f"{meta.normalized_cutoff_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_quiver(
    meta: ModeMeta,
    field_data: FieldData,
    tri: mtri.Triangulation,
    out_path: Path,
    dpi: int,
    quiver_nx: int,
    quiver_ny: int,
    show_mesh: bool,
    source: QuiverSource,
) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.2))

    if show_mesh:
        ax.triplot(tri, color="0.82", lw=0.20, alpha=0.65, zorder=1)

    qx, qy, qu, qv = _build_quiver_grid(
        meta,
        tri,
        field_data,
        quiver_nx,
        quiver_ny,
        source,
    )
    qun, qvn = _normalize_display_vectors(qu, qv)
    ax.quiver(
        qx,
        qy,
        qun,
        qvn,
        color="#0f4c5c",
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.0032,
        zorder=3,
    )

    _draw_geometry_outline(ax, meta)
    _set_axes_style(ax, meta)
    ax.set_title(
        f"HELM10 | {source.title}\n"
        f"{meta.plot_mode_line()}\n"
        f"{meta.normalized_cutoff_line()} | setas=({source.u_label}, {source.v_label})"
    )
    _save_figure(fig, out_path, dpi)


def _collect_modal_summary_data(
    rows: List[Dict[str, str]],
    value_key: str,
    use_abs: bool,
) -> List[Tuple[str, int, str, float, Optional[float]]]:
    data = []
    for row in rows:
        try:
            family = row["family"]
            mode_label = row["mode_label"]
            positive_rank = int(row["positive_rank"])
        except (KeyError, ValueError):
            continue
        value = _to_float(row.get(value_key))
        if value is None:
            continue
        kc_fem = _to_float(row.get("kc_fem"))
        data.append((family, positive_rank, mode_label, abs(value) if use_abs else value, kc_fem))
    return data


def _plot_modal_summary(
    case_root: Path,
    case_name: str,
    rows: List[Dict[str, str]],
    out_path: Path,
    dpi: int,
    *,
    value_key: str,
    use_abs: bool,
    y_label: str,
    family_title: str,
    figure_title: str,
    empty_warning: str,
    y_limits: Optional[Tuple[float, float]] = None,
) -> None:
    data = _collect_modal_summary_data(rows, value_key=value_key, use_abs=use_abs)

    if not data:
        print(f"Warning: {case_root} nao possui dados de {empty_warning}.")
        return

    families = sorted({family for family, _, _, _, _ in data})
    fig, axes = plt.subplots(
        len(families),
        1,
        figsize=(max(10.0, 0.6 * max(len(data), 6)), 3.6 * len(families)),
        squeeze=False,
    )

    for ax, family in zip(axes[:, 0], families):
        # Ordenamos por kc_fem para que a leitura siga a sequencia do espectro
        # numerico calculado. Em degenerescencias, usamos o rank positivo como
        # desempate estavel para manter a ordem reproduzivel.
        subset = sorted(
            (item for item in data if item[0] == family),
            key=lambda item: (
                float("inf") if item[4] is None else item[4],
                item[1],
                item[2],
            ),
        )
        labels = [mode_label for _, _, mode_label, _, _ in subset]
        values = [error_abs for _, _, _, error_abs, _ in subset]
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
        y_span = max(values) - min(values) if values else 0.0
        if y_limits is not None:
            y_offset = 0.02 * max(1e-12, y_limits[1] - y_limits[0])
        else:
            y_offset = 0.03 * max(1.0, max(values, default=0.0), y_span)
        for x, value in zip(xpos, values):
            ax.text(
                x,
                value + y_offset,
                f"{value:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    fig.suptitle(f"HELM10 | {figure_title} | {case_name.upper()}", fontsize=14)
    _save_figure(fig, out_path, dpi)


def _plot_error_summary(case_root: Path, case_name: str, rows: List[Dict[str, str]], out_path: Path, dpi: int) -> None:
    _plot_modal_summary(
        case_root=case_root,
        case_name=case_name,
        rows=rows,
        out_path=out_path,
        dpi=dpi,
        value_key="error_percent",
        use_abs=True,
        y_label="|error_percent| [%]",
        family_title="erro percentual absoluto por modo",
        figure_title="resumo de erro por modo",
        empty_warning="erro para o resumo",
    )


def _plot_rho_summary(case_root: Path, case_name: str, rows: List[Dict[str, str]], out_path: Path, dpi: int) -> None:
    _plot_modal_summary(
        case_root=case_root,
        case_name=case_name,
        rows=rows,
        out_path=out_path,
        dpi=dpi,
        value_key="rho_abs",
        use_abs=True,
        y_label="|rho|",
        family_title="correlacao modal por modo",
        figure_title="resumo de correlacao modal por modo",
        empty_warning="correlacao modal para o resumo",
        y_limits=(0.0, 1.05),
    )


def _iter_case_modes(case_root: Path, case_name: str) -> Iterable[ModeMeta]:
    modes_csv = case_root / "csv" / f"helm10_{case_name}_modes.csv"
    if not modes_csv.exists():
        raise FileNotFoundError(f"Arquivo ausente: {modes_csv}")
    for row in _read_rows(modes_csv):
        yield ModeMeta(case_name=case_name, case_root=case_root, row=row)


def _generate_case_images(
    case_root: Path,
    case_name: str,
    dpi: int,
    levels: int,
    quiver_nx: int,
    quiver_ny: int,
    show_mesh: bool,
    prefer_vtk: bool,
    quiver_source_key: str,
) -> None:
    print(f"Processing HELM10 case: {case_name}")
    if quiver_source_key == "all":
        quiver_sources = [QUIVER_SOURCES["field"], QUIVER_SOURCES["gradient"]]
    else:
        quiver_sources = [QUIVER_SOURCES[quiver_source_key]]
    img_root = case_root / "img"
    iso_dir = img_root / "isopotential"
    quiver_dir = img_root / "quiver"
    iso_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    img_root.mkdir(parents=True, exist_ok=True)

    mode_rows = _read_rows(case_root / "csv" / f"helm10_{case_name}_modes.csv")
    for meta in _iter_case_modes(case_root, case_name):
        field_data = _read_field_data(meta.fields_csv_path)
        tri = _build_triangulation(meta, field_data, prefer_vtk=prefer_vtk)

        iso_path = iso_dir / f"{meta.base_stem()}_psi.png"
        _plot_isopotential(
            meta=meta,
            field_data=field_data,
            tri=tri,
            out_path=iso_path,
            dpi=dpi,
            levels=levels,
            show_mesh=show_mesh,
        )
        for quiver_source in quiver_sources:
            quiver_path = quiver_dir / f"{meta.base_stem()}_{quiver_source.file_suffix}.png"
            _plot_quiver(
                meta=meta,
                field_data=field_data,
                tri=tri,
                out_path=quiver_path,
                dpi=dpi,
                quiver_nx=quiver_nx,
                quiver_ny=quiver_ny,
                show_mesh=show_mesh,
                source=quiver_source,
            )

    _plot_error_summary(
        case_root=case_root,
        case_name=case_name,
        rows=mode_rows,
        out_path=img_root / f"helm10_{case_name}_error_by_mode.png",
        dpi=dpi,
    )
    _plot_rho_summary(
        case_root=case_root,
        case_name=case_name,
        rows=mode_rows,
        out_path=img_root / f"helm10_{case_name}_rho_by_mode.png",
        dpi=dpi,
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate HELM10 images from the CSV outputs stored under out/helm10. "
            "The script saves isopotential plots, quiver plots and error summaries."
        )
    )
    parser.add_argument(
        "--root",
        default="out/helm10",
        help="Diretorio raiz com as saidas do HELM10. Padrao: out/helm10",
    )
    parser.add_argument(
        "--case",
        action="append",
        choices=["rect", "circle", "coax", "all"],
        help="Caso a processar. Pode ser repetido. Padrao: all",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="Resolucao das imagens salvas. Padrao: 220",
    )
    parser.add_argument(
        "--levels",
        type=int,
        default=18,
        help="Numero de niveis de isolinhas para o potencial escalar. Padrao: 18",
    )
    parser.add_argument(
        "--quiver-nx",
        type=int,
        default=20,
        help="Numero de pontos na direcao x da grade regular usada no quiver. Padrao: 20",
    )
    parser.add_argument(
        "--quiver-ny",
        type=int,
        default=10,
        help="Numero de pontos na direcao y da grade regular usada no quiver. Padrao: 10",
    )
    parser.add_argument(
        "--show-mesh",
        action="store_true",
        help="Sobrepoe a malha triangular nas imagens.",
    )
    parser.add_argument(
        "--quiver-source",
        choices=["all"] + sorted(QUIVER_SOURCES.keys()),
        default="all",
        help=(
            "Fonte vetorial do quiver. "
            "Use `field` para Ex/Ey, `gradient` para dpsi_dx/dpsi_dy, "
            "ou `all` para salvar ambos. Padrao: all"
        ),
    )
    parser.add_argument(
        "--no-vtk",
        action="store_true",
        help="Forca a triangulacao apenas a partir dos CSVs, sem reaproveitar o VTK.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    root = _resolve(Path(args.root))
    selected = args.case or ["all"]
    case_names = ["rect", "circle", "coax"] if "all" in selected else selected

    for case_name in case_names:
        case_root = root / case_name
        if not case_root.exists():
            raise FileNotFoundError(f"Caso ausente: {case_root}")
        _generate_case_images(
            case_root=case_root,
            case_name=case_name,
            dpi=args.dpi,
            levels=args.levels,
            quiver_nx=args.quiver_nx,
            quiver_ny=args.quiver_ny,
            show_mesh=args.show_mesh,
            prefer_vtk=not args.no_vtk,
            quiver_source_key=args.quiver_source,
        )


if __name__ == "__main__":
    main()
