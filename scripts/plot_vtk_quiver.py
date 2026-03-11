#!/usr/bin/env python3
"""
Plot 2D vector fields from legacy ASCII VTK files generated in this repo.

Main modes:
1) Single file plot:
   python3 scripts/plot_vtk_quiver.py out/2d/2.1_scalar/circle/te_circle_sv.vtk --out out/img/te_circle.png

2) Batch generation:
   python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out/2d --out-dir out/img_all
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class VtkTriData:
    points: np.ndarray
    tris: np.ndarray
    point_vectors: Optional[np.ndarray] = None
    cell_vectors: Optional[np.ndarray] = None
    vector_name: str = "V"


@dataclass
class ModeRow:
    geometry: str
    formulation: str
    pol: str
    field: str
    mode_rank: int
    m: int
    n_or_p: int
    index2_name: str
    kc_ana: float
    kc_fem: float
    err_pct: float
    rho: float
    a: Optional[float] = None
    b: Optional[float] = None
    r1: Optional[float] = None
    r2: Optional[float] = None
    R: Optional[float] = None

    def as_csv_row(self) -> Dict[str, object]:
        row = {
            "geometry": self.geometry,
            "formulation": self.formulation,
            "polarization": self.pol,
            "field": self.field,
            "mode_rank": self.mode_rank,
            "m": self.m,
            self.index2_name: self.n_or_p,
            "index2_name": self.index2_name,
            "kc_analytic": self.kc_ana,
            "kc_fem": self.kc_fem,
            "rel_error_pct": self.err_pct,
            "rho_abs": self.rho,
            "a": self.a,
            "b": self.b,
            "R": self.R,
            "r1": self.r1,
            "r2": self.r2,
            "kc_times_a": self.kc_fem * self.a if self.a else None,
            "kc_times_R": self.kc_fem * self.R if self.R else None,
            "kc_times_r1": self.kc_fem * self.r1 if self.r1 else None,
            "kc_times_r2": self.kc_fem * self.r2 if self.r2 else None,
        }
        return row


def _find_line_index(lines: List[str], prefix: str) -> int:
    for i, ln in enumerate(lines):
        if ln.strip().startswith(prefix):
            return i
    return -1


def _read_legacy_vtk_ascii(path: Path) -> VtkTriData:
    lines = path.read_text(encoding="utf-8").splitlines()

    i_points = _find_line_index(lines, "POINTS ")
    if i_points < 0:
        raise ValueError(f"{path}: missing POINTS block")
    n_points = int(lines[i_points].split()[1])
    pts = np.array(
        [[float(v) for v in lines[i_points + 1 + k].split()[:3]] for k in range(n_points)],
        dtype=float,
    )
    points = pts[:, :2]

    i_cells = _find_line_index(lines, "CELLS ")
    if i_cells < 0:
        raise ValueError(f"{path}: missing CELLS block")
    n_cells = int(lines[i_cells].split()[1])
    tris = []
    for k in range(n_cells):
        vals = [int(v) for v in lines[i_cells + 1 + k].split()]
        if vals[0] != 3:
            raise ValueError(f"{path}: only triangular cells are supported")
        tris.append(vals[1:4])
    tris_arr = np.array(tris, dtype=int)

    point_vectors = None
    cell_vectors = None
    vector_name = "V"

    i_pd = _find_line_index(lines, "POINT_DATA ")
    if i_pd >= 0:
        for j in range(i_pd + 1, len(lines)):
            ln = lines[j].strip()
            if ln.startswith("VECTORS "):
                parts = ln.split()
                vector_name = parts[1]
                vec = np.array(
                    [[float(v) for v in lines[j + 1 + k].split()[:3]] for k in range(n_points)],
                    dtype=float,
                )
                point_vectors = vec[:, :2]
                break

    i_cd = _find_line_index(lines, "CELL_DATA ")
    if i_cd >= 0:
        for j in range(i_cd + 1, len(lines)):
            ln = lines[j].strip()
            if ln.startswith("VECTORS "):
                parts = ln.split()
                vector_name = parts[1]
                vec = np.array(
                    [[float(v) for v in lines[j + 1 + k].split()[:3]] for k in range(n_cells)],
                    dtype=float,
                )
                cell_vectors = vec[:, :2]
                break

    return VtkTriData(
        points=points,
        tris=tris_arr,
        point_vectors=point_vectors,
        cell_vectors=cell_vectors,
        vector_name=vector_name,
    )


def _cell_centers(points: np.ndarray, tris: np.ndarray) -> np.ndarray:
    return points[tris].mean(axis=1)


def _thin_indices(n: int, step: int) -> np.ndarray:
    return np.arange(0, n, max(1, step), dtype=int)


def _metadata_from_filename(vtk_name: str) -> Tuple[str, str, str]:
    lower = vtk_name.lower()
    geometry = "rect" if "rect" in lower else ("circle" if "circle" in lower else ("coax" if "coax" in lower else "unknown"))
    formulation = "edge" if lower.startswith("edge_") else "scalar"
    if "_et" in lower or "_te_" in lower or lower.startswith("te"):
        pol = "TE"
    elif "_ht" in lower or "_tm_" in lower or lower.startswith("tm"):
        pol = "TM"
    else:
        pol = "?"
    return geometry, formulation, pol


def _field_label(formulation: str, pol: str, vtk_name: str) -> str:
    lower = vtk_name.lower()
    if formulation == "edge":
        if "_et" in lower:
            return "Et"
        if "_ht" in lower:
            return "Ht"
    if pol == "TE":
        return "Et"
    if pol == "TM":
        return "Ht"
    return "Ft"


def _build_title(vtk_name: str, row: Optional[ModeRow]) -> str:
    geometry, formulation, pol = _metadata_from_filename(vtk_name)
    field = _field_label(formulation, pol, vtk_name)
    lhs = f"{formulation.upper()} | {geometry.upper()} | {pol} | {field}"

    tags = _mode_tags_from_filename(vtk_name)
    if row is None:
        if tags is not None:
            m, idx2_name, idx2_val, rank = tags
            return f"{lhs}\nmode=({m},{idx2_val}) [{idx2_name}] | rank={rank}"
        return lhs
    idx2 = row.index2_name
    suffix = f"mode=({row.m},{row.n_or_p}) [{idx2}] | kc_fem={row.kc_fem:.6f}"
    if row.geometry == "rect" and row.a:
        suffix += f" | kc*a={row.kc_fem * row.a:.6f}"
    elif row.geometry == "circle" and row.R:
        suffix += f" | kc*R={row.kc_fem * row.R:.6f}"
    elif row.geometry == "coax":
        if row.r1:
            suffix += f" | kc*r1={row.kc_fem * row.r1:.6f}"
        if row.r2:
            suffix += f" | kc*r2={row.kc_fem * row.r2:.6f}"
    return f"{lhs}\n{suffix}"


def _mode_tags_from_filename(vtk_name: str) -> Optional[Tuple[int, str, int, int]]:
    """
    Parse filenames like:
    - te_rect_m1_n0_rank01_sv.vtk
    - edge_circle_tm_m2_p1_rank03_Ht.vtk
    Returns: (m, index2_name, index2_value, rank)
    """
    m = re.search(r"_m(\d+)_(n|p)(\d+)_rank(\d+)", vtk_name.lower())
    if not m:
        return None
    return int(m.group(1)), m.group(2), int(m.group(3)), int(m.group(4))


def plot_quiver(
    vtk_path: Path,
    out_path: Optional[Path],
    stride: int,
    scale: float,
    show_mesh: bool,
    title: Optional[str],
    mode_row: Optional[ModeRow] = None,
    dpi: int = 210,
) -> None:
    data = _read_legacy_vtk_ascii(vtk_path)

    tri = mtri.Triangulation(data.points[:, 0], data.points[:, 1], data.tris)
    fig, ax = plt.subplots(figsize=(8.2, 6.3))

    if show_mesh:
        ax.triplot(tri, lw=0.28, color="0.83", zorder=0)

    if data.point_vectors is not None:
        p = data.points
        v = data.point_vectors
    elif data.cell_vectors is not None:
        p = _cell_centers(data.points, data.tris)
        v = data.cell_vectors
    else:
        raise ValueError(f"{vtk_path}: no VECTORS found in POINT_DATA or CELL_DATA")

    mag = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2)
    idx = _thin_indices(len(p), stride)

    q = ax.quiver(
        p[idx, 0],
        p[idx, 1],
        v[idx, 0],
        v[idx, 1],
        mag[idx],
        angles="xy",
        scale_units="xy",
        scale=scale,
        cmap="viridis",
        width=0.003,
        pivot="mid",
        zorder=2,
    )
    cbar = fig.colorbar(q, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("|F|")

    ttl = title or _build_title(vtk_path.name, mode_row)
    ax.set_title(ttl, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    xmin, xmax = data.points[:, 0].min(), data.points[:, 0].max()
    ymin, ymax = data.points[:, 1].min(), data.points[:, 1].max()
    pad_x = 0.03 * max(1e-12, xmax - xmin)
    pad_y = 0.03 * max(1e-12, ymax - ymin)
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)

    fig.tight_layout()
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi)
        print(f"Saved: {out_path}")
    else:
        plt.show()
    plt.close(fig)


def _run_cmd(cmd: List[str], cwd: Path, env: Optional[Dict[str, str]] = None) -> str:
    p = subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=True, env=env)
    return p.stdout


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _parse_mode_table(
    text: str,
    geometry: str,
    formulation: str,
    index2_name: str,
    dims: Dict[str, float],
) -> List[ModeRow]:
    rows: List[ModeRow] = []
    current_pol: Optional[str] = None

    if formulation == "scalar":
        field_for_pol = {"TE": "Et", "TM": "Ht"}
    else:
        field_for_pol = {"TE": "Et", "TM": "Ht"}

    line_re = re.compile(
        r"^\s*(\d+)\s+\((\d+),(\d+)\)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$"
    )

    for ln in text.splitlines():
        l = ln.strip()
        if l.startswith("[TE "):
            current_pol = "TE"
            continue
        if l.startswith("[TM "):
            current_pol = "TM"
            continue
        m = line_re.match(ln)
        if not m or current_pol is None:
            continue
        rows.append(
            ModeRow(
                geometry=geometry,
                formulation=formulation,
                pol=current_pol,
                field=field_for_pol[current_pol],
                mode_rank=int(m.group(1)),
                m=int(m.group(2)),
                n_or_p=int(m.group(3)),
                index2_name=index2_name,
                kc_ana=float(m.group(4)),
                kc_fem=float(m.group(5)),
                err_pct=float(m.group(6)),
                rho=float(m.group(7)),
                a=dims.get("a"),
                b=dims.get("b"),
                R=dims.get("R"),
                r1=dims.get("r1"),
                r2=dims.get("r2"),
            )
        )
    return rows


def _parse_dims(text: str) -> Dict[str, float]:
    dims: Dict[str, float] = {}
    m_rect = re.search(r"\ba=([0-9.+-eE]+)\s+b=([0-9.+-eE]+)\b", text)
    if m_rect:
        dims["a"] = float(m_rect.group(1))
        dims["b"] = float(m_rect.group(2))
    m_circle = re.search(r"\bR=([0-9.+-eE]+)\b", text)
    if m_circle:
        dims["R"] = float(m_circle.group(1))
    m_coax = re.search(r"\br1=([0-9.+-eE]+)\s+r2=([0-9.+-eE]+)\b", text)
    if m_coax:
        dims["r1"] = float(m_coax.group(1))
        dims["r2"] = float(m_coax.group(2))
    return dims


def _collect_all_mode_rows(build_dir: Path, mode_export: int, out_root: Path) -> List[ModeRow]:
    runs = [
        (["./helm10_rect", "14", "14", f"{mode_export}"], "rect", "scalar", "n"),
        (["./helm10_circle", "10", "48", f"{mode_export}"], "circle", "scalar", "p"),
        (["./helm10_coax", "10", "48", f"{mode_export}"], "coax", "scalar", "p"),
        (["./edge_rect", "14", "14", f"{mode_export}"], "rect", "edge", "n"),
        (["./edge_circle", "10", "48", f"{mode_export}"], "circle", "edge", "p"),
        (["./edge_coax", "10", "48", f"{mode_export}"], "coax", "edge", "p"),
    ]
    all_rows: List[ModeRow] = []
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(out_root)
    for cmd, geometry, formulation, idx2 in runs:
        print(f"Running: {' '.join(cmd)}")
        out = _run_cmd(cmd, cwd=build_dir, env=env)
        dims = _parse_dims(out)
        all_rows.extend(_parse_mode_table(out, geometry, formulation, idx2, dims))
    return all_rows


def _write_csv(rows: List[ModeRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "geometry",
        "formulation",
        "polarization",
        "field",
        "mode_rank",
        "m",
        "n",
        "p",
        "index2_name",
        "kc_analytic",
        "kc_fem",
        "rel_error_pct",
        "rho_abs",
        "a",
        "b",
        "R",
        "r1",
        "r2",
        "kc_times_a",
        "kc_times_R",
        "kc_times_r1",
        "kc_times_r2",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            d = r.as_csv_row()
            d.setdefault("n", "")
            d.setdefault("p", "")
            if r.index2_name == "n":
                d["n"] = r.n_or_p
                d["p"] = ""
            else:
                d["p"] = r.n_or_p
                d["n"] = ""
            w.writerow(d)
    print(f"Saved CSV: {path}")


def _first_mode_lookup(rows: List[ModeRow]) -> Dict[Tuple[str, str, str], ModeRow]:
    out: Dict[Tuple[str, str, str], ModeRow] = {}
    for r in rows:
        key = (r.geometry, r.formulation, r.pol)
        old = out.get(key)
        if old is None or r.mode_rank < old.mode_rank:
            out[key] = r
    return out


def _rank_mode_lookup(rows: List[ModeRow]) -> Dict[Tuple[str, str, str, int], ModeRow]:
    out: Dict[Tuple[str, str, str, int], ModeRow] = {}
    for r in rows:
        key = (r.geometry, r.formulation, r.pol, r.mode_rank)
        out[key] = r
    return out


def _plot_all_images(
    vtk_root: Path,
    out_dir: Path,
    rows: List[ModeRow],
    stride: int,
    scale: float,
    show_mesh: bool,
    dpi: int,
    max_rank: int,
) -> None:
    first = _first_mode_lookup(rows)
    by_rank = _rank_mode_lookup(rows)

    if not vtk_root.exists():
        raise SystemExit(f"VTK root directory not found: {vtk_root}")

    vtk_files = sorted(vtk_root.rglob("*.vtk"))
    if not vtk_files:
        raise SystemExit(f"No VTK files found under: {vtk_root}")

    for vtk_path in vtk_files:
        name = vtk_path.name
        geometry, formulation, pol = _metadata_from_filename(name)

        mode_row: Optional[ModeRow]
        tags = _mode_tags_from_filename(name)
        if tags is not None:
            _, _, _, rank = tags
            if max_rank > 0 and rank > max_rank:
                continue
            mode_row = by_rank.get((geometry, formulation, pol, rank))
        else:
            mode_row = first.get((geometry, formulation, pol))

        rel = vtk_path.relative_to(vtk_root)
        img_path = out_dir / rel.with_suffix(".png")
        plot_quiver(
            vtk_path=vtk_path,
            out_path=img_path,
            stride=stride,
            scale=scale,
            show_mesh=show_mesh,
            title=None,
            mode_row=mode_row,
            dpi=dpi,
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot quiver from legacy ASCII VTK.")
    p.add_argument("vtk", nargs="?", type=Path, default=None, help="Input .vtk file (single-file mode).")
    p.add_argument("--out", type=Path, default=None, help="Output image path for single-file mode.")
    p.add_argument("--stride", type=int, default=2, help="Arrow subsampling step.")
    p.add_argument("--scale", type=float, default=22.0, help="Quiver scale (smaller => longer arrows).")
    p.add_argument("--no-mesh", action="store_true", help="Hide triangular mesh lines.")
    p.add_argument("--title", type=str, default=None, help="Custom title (single-file mode).")
    p.add_argument("--all-img", "--all_img", dest="all_img", action="store_true", help="Run all solvers, generate all images and CSV.")
    p.add_argument("--build-dir", type=Path, default=Path("build"), help="Build/output directory.")
    p.add_argument("--vtk-root", type=Path, default=Path("out/2d"), help="Root folder containing VTK files for batch plotting.")
    p.add_argument("--out-dir", type=Path, default=Path("out/img_all"), help="Folder for batch images.")
    p.add_argument("--csv", type=Path, default=Path("out/img_all/mode_summary.csv"), help="CSV output path in batch mode.")
    p.add_argument("--mode-export", type=int, default=8, help="Mode export count forwarded to 2D executables in --all-img mode.")
    p.add_argument("--max-rank", type=int, default=8, help="Maximum rank to render from filenames with _rankXX_ (<=0 means all).")
    p.add_argument("--dpi", type=int, default=210, help="Output image DPI.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    build_dir = _resolve(args.build_dir)
    vtk_root = _resolve(args.vtk_root)
    out_dir = _resolve(args.out_dir)
    csv_path = _resolve(args.csv)

    if args.all_img:
        rows = _collect_all_mode_rows(build_dir, args.mode_export, vtk_root.parent)
        _write_csv(rows, csv_path)
        _plot_all_images(
            vtk_root=vtk_root,
            out_dir=out_dir,
            rows=rows,
            stride=args.stride,
            scale=args.scale,
            show_mesh=not args.no_mesh,
            dpi=args.dpi,
            max_rank=args.max_rank,
        )
        return

    if args.vtk is None:
        raise SystemExit("Provide a .vtk file or use --all-img (alias: --all_img).")

    plot_quiver(
        vtk_path=_resolve(args.vtk),
        out_path=_resolve(args.out) if args.out else None,
        stride=args.stride,
        scale=args.scale,
        show_mesh=not args.no_mesh,
        title=args.title,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
