#!/usr/bin/env python3
"""
Generate HELMVEC images directly from the CSV outputs produced by the solvers.

The script scans `out/helmvec/<case>/csv`, reads:

- `edge_<case>_modes.csv`
- `edge_<case>_fields_<mode>.csv`

and produces, for each exported mode:

- a magnitude map of the transverse field on the triangular cells
- a quiver plot of the transverse field on cell centroids

It also creates summary images per case using:

- `error_percent`
- `rho_abs`
"""

from __future__ import annotations

import argparse
import csv
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


ROOT = Path(__file__).resolve().parents[1]


@dataclass
class FieldData:
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
class ModeMeta:
    case_name: str
    case_root: Path
    row: Dict[str, str]

    @property
    def family(self) -> str:
        return self.row["family"]

    @property
    def transverse_label(self) -> str:
        return self.row["transverse_label"]

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
    def m(self) -> int:
        return int(self.row["m"])

    @property
    def index2_name(self) -> str:
        return "n" if self.case_name == "rect" else "p"

    @property
    def index2_value(self) -> int:
        return int(self.row[self.index2_name])

    @property
    def fields_csv_name(self) -> str:
        return self.row["fields_csv_file"]

    @property
    def fields_csv_path(self) -> Path:
        return self.case_root / "csv" / self.fields_csv_name

    @property
    def vtk_name(self) -> str:
        return self.row["vtk_file"]

    @property
    def vtk_path(self) -> Path:
        return self.case_root / "vtk" / self.vtk_name

    @property
    def kc_fem(self) -> Optional[float]:
        return _to_float(self.row.get("kc_fem"))

    def short_name(self) -> str:
        return f"{self.mode_label}_rank{self.positive_rank:02d}"

    def base_stem(self) -> str:
        return f"helmvec_{self.case_name}_{self.short_name()}"

    def plot_mode_line(self) -> str:
        return (
            f"{self.case_name.upper()} | "
            f"{self.family}$_{{{self.m},{self.index2_value}}}$ | "
            f"{self.transverse_label}"
        )

    def normalized_cutoff_line(self) -> str:
        if self.case_name == "rect":
            value = _to_float(self.row.get("kc_ar_fem"))
            if value is None and self.kc_fem is not None:
                value = self.kc_fem * (_to_float(self.row.get("ar_m")) or 0.0)
            return f"$k_c a_r = {value:.6f}$" if value is not None else "$k_c a_r = ?$"
        if self.case_name == "circle":
            value = _to_float(self.row.get("kc_r_fem"))
            if value is None and self.kc_fem is not None:
                value = self.kc_fem * (_to_float(self.row.get("r_m")) or 0.0)
            return f"$k_c r = {value:.6f}$" if value is not None else "$k_c r = ?$"
        value = _to_float(self.row.get("kc_r1_fem"))
        if value is None and self.kc_fem is not None:
            value = self.kc_fem * (_to_float(self.row.get("r1_m")) or 0.0)
        return f"$k_c r_1 = {value:.6f}$" if value is not None else "$k_c r_1 = ?$"


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
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    if not rows:
        raise ValueError(f"{csv_path}: CSV de campos vazio")
    rows.sort(key=lambda row: int(row["cell_id"]))

    if {"Ex", "Ey", "Emag"}.issubset(fieldnames):
        x_name, y_name, mag_name, symbol = "Ex", "Ey", "Emag", "E"
    elif {"Hx", "Hy", "Hmag"}.issubset(fieldnames):
        x_name, y_name, mag_name, symbol = "Hx", "Hy", "Hmag", "H"
    else:
        x_name, y_name, mag_name, symbol = "Fx", "Fy", "Fmag", "F"

    def col(name: str) -> np.ndarray:
        return np.array([float(row[name]) for row in rows], dtype=float)

    return FieldData(
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


def _save_figure(fig: plt.Figure, out_path: Path, dpi: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.97))
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


def _plot_magnitude(
    meta: ModeMeta,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: FieldData,
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
        "HELMVEC | Magnitude do Campo Transversal\n"
        f"{meta.plot_mode_line()}\n"
        f"{meta.normalized_cutoff_line()}"
    )
    _save_figure(fig, out_path, dpi)


def _plot_quiver(
    meta: ModeMeta,
    points: np.ndarray,
    triangles: np.ndarray,
    field_data: FieldData,
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
        "HELMVEC | Distribuicao do Campo Transversal\n"
        f"{meta.plot_mode_line()}\n"
        f"{meta.normalized_cutoff_line()} | setas=({field_data.x_name}, {field_data.y_name})"
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
            ax.text(
                x,
                value + y_offset,
                f"{value:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    fig.suptitle(f"HELMVEC | {figure_title} | {case_name.upper()}", fontsize=14)
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
    modes_csv = case_root / "csv" / f"edge_{case_name}_modes.csv"
    if not modes_csv.exists():
        raise FileNotFoundError(f"Arquivo ausente: {modes_csv}")
    for row in _read_rows(modes_csv):
        yield ModeMeta(case_name=case_name, case_root=case_root, row=row)


def _generate_case_images(case_root: Path, case_name: str, dpi: int, show_mesh: bool, max_arrows: int) -> None:
    print(f"Processing HELMVEC case: {case_name}")
    img_root = case_root / "img"
    magnitude_dir = img_root / "magnitude"
    quiver_dir = img_root / "quiver"
    magnitude_dir.mkdir(parents=True, exist_ok=True)
    quiver_dir.mkdir(parents=True, exist_ok=True)
    img_root.mkdir(parents=True, exist_ok=True)

    mode_rows = _read_rows(case_root / "csv" / f"edge_{case_name}_modes.csv")
    for meta in _iter_case_modes(case_root, case_name):
        field_data = _read_field_data(meta.fields_csv_path)
        points, triangles = _read_legacy_vtk_connectivity(meta.vtk_path)
        magnitude_path = magnitude_dir / f"{meta.base_stem()}_magnitude.png"
        quiver_path = quiver_dir / f"{meta.base_stem()}_quiver.png"
        _plot_magnitude(meta, points, triangles, field_data, magnitude_path, dpi, show_mesh)
        _plot_quiver(meta, points, triangles, field_data, quiver_path, dpi, show_mesh, max_arrows)

    _plot_error_summary(
        case_root=case_root,
        case_name=case_name,
        rows=mode_rows,
        out_path=img_root / f"helmvec_{case_name}_error_by_mode.png",
        dpi=dpi,
    )
    _plot_rho_summary(
        case_root=case_root,
        case_name=case_name,
        rows=mode_rows,
        out_path=img_root / f"helmvec_{case_name}_rho_by_mode.png",
        dpi=dpi,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Le os CSVs do HELMVEC e gera mapas de magnitude, quivers e "
            "resumos de erro/rho por modo."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("out/helmvec"),
        help="Diretorio raiz das saidas do HELMVEC. Padrao: out/helmvec",
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
        _generate_case_images(
            case_root=case_root,
            case_name=case_name,
            dpi=args.dpi,
            show_mesh=args.show_mesh,
            max_arrows=args.max_arrows,
        )


if __name__ == "__main__":
    main()
