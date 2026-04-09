#!/usr/bin/env python3
"""
Diagnose mesh/material staircasing for Section 2.2.4 Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_interface_staircasing.csv
- out/validation/DIAG_224_INTERFACE_STAIRCASING.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_CSV = Path("out/validation/diag_224_interface_staircasing.csv")
DEFAULT_OUT_MD = Path("out/validation/DIAG_224_INTERFACE_STAIRCASING.md")
DEFAULT_CONVERGENCE_CSV = Path("out/validation/diag_224_table10_convergence.csv")

A = 1.0
B = 0.45
SECTION = "2.2.4"
CASE = "Figure13_Table10"


@dataclass(frozen=True)
class PrioritySpec:
    priority_point: str
    d_over_a: float
    br_over_lambda0: float

    @property
    def mode(self) -> str:
        return f"d/a={self.d_over_a},br/lambda0={self.br_over_lambda0}"

    @property
    def x_split(self) -> float:
        return self.d_over_a * A


@dataclass(frozen=True)
class MeshSpec:
    nx: int
    ny: int

    @property
    def mesh_label(self) -> str:
        return f"nx{self.nx}_ny{self.ny}"

    @property
    def hx(self) -> float:
        return A / self.nx

    @property
    def hy(self) -> float:
        return B / self.ny

    @property
    def mesh_nodes(self) -> int:
        return (self.nx + 1) * (self.ny + 1)

    @property
    def mesh_tris(self) -> int:
        return 2 * self.nx * self.ny


@dataclass(frozen=True)
class Point2D:
    x: float
    y: float


@dataclass(frozen=True)
class Triangle:
    cell_i: int
    cell_j: int
    local_id: int
    vertices: tuple[Point2D, Point2D, Point2D]

    @property
    def area(self) -> float:
        return abs(_polygon_area(self.vertices))

    @property
    def centroid_x(self) -> float:
        return sum(point.x for point in self.vertices) / 3.0


PRIORITY_SPECS = (
    PrioritySpec("P1", 0.167, 0.5),
    PrioritySpec("P2", 0.286, 0.5),
    PrioritySpec("P3", 0.5, 0.4),
)

MESH_SPECS = (
    MeshSpec(10, 5),
    MeshSpec(16, 8),
    MeshSpec(20, 10),
)


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _to_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _polygon_area(points: tuple[Point2D, ...] | list[Point2D]) -> float:
    if len(points) < 3:
        return 0.0
    area2 = 0.0
    for idx, point in enumerate(points):
        nxt = points[(idx + 1) % len(points)]
        area2 += point.x * nxt.y - nxt.x * point.y
    return 0.5 * area2


def _clip_polygon_left_of_vertical(points: tuple[Point2D, ...], x_split: float) -> list[Point2D]:
    if not points:
        return []

    def inside(point: Point2D) -> bool:
        return point.x <= x_split + 1e-15

    def intersection(start: Point2D, end: Point2D) -> Point2D:
        dx = end.x - start.x
        if abs(dx) < 1e-15:
            return Point2D(x_split, start.y)
        t = (x_split - start.x) / dx
        return Point2D(x_split, start.y + t * (end.y - start.y))

    output = list(points)
    clipped: list[Point2D] = []
    prev = output[-1]
    prev_inside = inside(prev)
    for current in output:
        current_inside = inside(current)
        if current_inside:
            if not prev_inside:
                clipped.append(intersection(prev, current))
            clipped.append(current)
        elif prev_inside:
            clipped.append(intersection(prev, current))
        prev = current
        prev_inside = current_inside
    return clipped


def _make_rect_triangles(mesh: MeshSpec) -> list[Triangle]:
    triangles: list[Triangle] = []
    for j in range(mesh.ny):
        y0 = j * mesh.hy
        y1 = (j + 1) * mesh.hy
        for i in range(mesh.nx):
            x0 = i * mesh.hx
            x1 = (i + 1) * mesh.hx
            n00 = Point2D(x0, y0)
            n10 = Point2D(x1, y0)
            n01 = Point2D(x0, y1)
            n11 = Point2D(x1, y1)
            triangles.append(Triangle(i, j, 0, (n00, n10, n11)))
            triangles.append(Triangle(i, j, 1, (n00, n11, n01)))
    return triangles


def _load_convergence_rows(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if path is None or not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row.get("priority_point", "").strip(), row.get("mesh_label", "").strip())
        if key == ("", ""):
            continue
        index[key] = row
    return index


def _cut_column_metrics(spec: PrioritySpec, mesh: MeshSpec) -> dict[str, str]:
    x_split = spec.x_split
    hx = mesh.hx
    position = x_split / hx
    nearest = round(position)
    aligned = abs(position - nearest) < 1e-12

    if aligned:
        return {
            "interface_aligned_to_mesh_edge": "1",
            "cut_column_index": "",
            "cut_column_x0": f"{x_split:.16g}",
            "cut_column_x1": f"{x_split:.16g}",
            "cut_column_local_phase": "",
            "cut_column_exact_fill_fraction": "",
            "cut_column_centroid_fill_fraction": "",
            "cut_column_fill_fraction_error": "",
            "cut_column_centroid_left_count": "0",
            "cut_column_centroid_xs": "",
            "cut_column_mode": "aligned",
        }

    col = int(math.floor(position))
    x0 = col * hx
    x1 = (col + 1) * hx
    phase = (x_split - x0) / hx
    centroid_xs = (x0 + hx / 3.0, x0 + 2.0 * hx / 3.0)
    left_count = sum(1 for value in centroid_xs if value < x_split)
    fill_fraction = 0.5 * left_count

    if left_count == 0:
        mode = "empty"
    elif left_count == 1:
        mode = "half"
    else:
        mode = "full"

    return {
        "interface_aligned_to_mesh_edge": "0",
        "cut_column_index": str(col),
        "cut_column_x0": f"{x0:.16g}",
        "cut_column_x1": f"{x1:.16g}",
        "cut_column_local_phase": f"{phase:.16g}",
        "cut_column_exact_fill_fraction": f"{phase:.16g}",
        "cut_column_centroid_fill_fraction": f"{fill_fraction:.16g}",
        "cut_column_fill_fraction_error": f"{abs(fill_fraction - phase):.16g}",
        "cut_column_centroid_left_count": str(left_count),
        "cut_column_centroid_xs": "|".join(f"{value:.12g}" for value in centroid_xs),
        "cut_column_mode": mode,
    }


def build_rows(convergence_rows: dict[tuple[str, str], dict[str, str]]) -> list[dict[str, str]]:
    domain_area = A * B
    triangles_by_mesh = {mesh.mesh_label: _make_rect_triangles(mesh) for mesh in MESH_SPECS}
    output_rows: list[dict[str, str]] = []

    for spec in PRIORITY_SPECS:
        target_fill_area = spec.x_split * B
        target_fill_fraction = target_fill_area / domain_area
        for mesh in MESH_SPECS:
            triangles = triangles_by_mesh[mesh.mesh_label]
            effective_fill_area = 0.0
            symmetric_diff_area = 0.0
            cut_triangle_count = 0
            cut_cells: set[tuple[int, int]] = set()
            represented_left_count = 0

            for triangle in triangles:
                tri_area = triangle.area
                represented = tri_area if triangle.centroid_x < spec.x_split else 0.0
                if represented > 0.0:
                    represented_left_count += 1
                clipped = _clip_polygon_left_of_vertical(triangle.vertices, spec.x_split)
                exact_area = abs(_polygon_area(clipped))
                effective_fill_area += represented
                symmetric_diff_area += abs(represented - exact_area)
                if 1e-15 < exact_area < tri_area - 1e-15:
                    cut_triangle_count += 1
                    cut_cells.add((triangle.cell_i, triangle.cell_j))

            effective_fill_fraction = effective_fill_area / domain_area
            equivalent_interface_x = effective_fill_area / B
            interface_x_abs_error = abs(equivalent_interface_x - spec.x_split)
            fill_fraction_abs_error = abs(effective_fill_fraction - target_fill_fraction)

            cut_metrics = _cut_column_metrics(spec, mesh)
            conv_row = convergence_rows.get((spec.priority_point, mesh.mesh_label), {})

            output_rows.append(
                {
                    "section": SECTION,
                    "case": CASE,
                    "priority_point": spec.priority_point,
                    "mode": spec.mode,
                    "d_over_a_target": f"{spec.d_over_a:.16g}",
                    "br_over_lambda0": f"{spec.br_over_lambda0:.16g}",
                    "mesh_label": mesh.mesh_label,
                    "nx": str(mesh.nx),
                    "ny": str(mesh.ny),
                    "mesh_nodes": str(mesh.mesh_nodes),
                    "mesh_tris": str(mesh.mesh_tris),
                    "hx": f"{mesh.hx:.16g}",
                    "hy": f"{mesh.hy:.16g}",
                    "x_split_target": f"{spec.x_split:.16g}",
                    "target_fill_area": f"{target_fill_area:.16g}",
                    "target_fill_fraction": f"{target_fill_fraction:.16g}",
                    "effective_fill_area_centroid": f"{effective_fill_area:.16g}",
                    "effective_fill_fraction_centroid": f"{effective_fill_fraction:.16g}",
                    "equivalent_interface_x_centroid": f"{equivalent_interface_x:.16g}",
                    "interface_x_abs_error": f"{interface_x_abs_error:.16g}",
                    "fill_fraction_abs_error": f"{fill_fraction_abs_error:.16g}",
                    "symmetric_diff_area": f"{symmetric_diff_area:.16g}",
                    "symmetric_diff_fraction_domain": f"{(symmetric_diff_area / domain_area):.16g}",
                    "symmetric_diff_fraction_target": (
                        f"{(symmetric_diff_area / target_fill_area):.16g}"
                        if target_fill_area > 0.0
                        else ""
                    ),
                    "represented_left_triangles": str(represented_left_count),
                    "cut_triangle_count": str(cut_triangle_count),
                    "cut_cell_count": str(len(cut_cells)),
                    **cut_metrics,
                    "convergence_csv": conv_row.get("table10_csv", ""),
                    "fem": conv_row.get("fem", ""),
                    "err_primary_pct": conv_row.get("err_primary_pct", ""),
                    "err_secondary_pct": conv_row.get("err_secondary_pct", ""),
                    "selected_rank": conv_row.get("selected_rank", ""),
                    "selected_eig_index": conv_row.get("selected_eig_index", ""),
                    "ez_ratio": conv_row.get("ez_ratio", ""),
                    "match_status": conv_row.get("match_status", ""),
                    "field_status": conv_row.get("field_status", ""),
                }
            )

    return output_rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "section",
        "case",
        "priority_point",
        "mode",
        "d_over_a_target",
        "br_over_lambda0",
        "mesh_label",
        "nx",
        "ny",
        "mesh_nodes",
        "mesh_tris",
        "hx",
        "hy",
        "x_split_target",
        "target_fill_area",
        "target_fill_fraction",
        "effective_fill_area_centroid",
        "effective_fill_fraction_centroid",
        "equivalent_interface_x_centroid",
        "interface_x_abs_error",
        "fill_fraction_abs_error",
        "symmetric_diff_area",
        "symmetric_diff_fraction_domain",
        "symmetric_diff_fraction_target",
        "represented_left_triangles",
        "cut_triangle_count",
        "cut_cell_count",
        "interface_aligned_to_mesh_edge",
        "cut_column_index",
        "cut_column_x0",
        "cut_column_x1",
        "cut_column_local_phase",
        "cut_column_exact_fill_fraction",
        "cut_column_centroid_fill_fraction",
        "cut_column_fill_fraction_error",
        "cut_column_centroid_left_count",
        "cut_column_centroid_xs",
        "cut_column_mode",
        "convergence_csv",
        "fem",
        "err_primary_pct",
        "err_secondary_pct",
        "selected_rank",
        "selected_eig_index",
        "ez_ratio",
        "match_status",
        "field_status",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_md(
    path: Path,
    rows: list[dict[str, str]],
    *,
    csv_path: Path,
    convergence_path: Path | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_by_point: dict[str, list[dict[str, str]]] = {spec.priority_point: [] for spec in PRIORITY_SPECS}
    for row in rows:
        rows_by_point[row["priority_point"]].append(row)
    for point_rows in rows_by_point.values():
        point_rows.sort(key=lambda row: (int(row["nx"]), int(row["ny"])))

    lines: list[str] = []
    lines.append("# Figure 13 / Table 10 Interface Staircasing Diagnostic")
    lines.append("")
    lines.append(f"- Generated at: `{dt.datetime.now().isoformat(timespec='seconds')}`")
    lines.append(f"- Section: `{SECTION}`")
    lines.append(f"- Case: `{CASE}`")
    lines.append(f"- Geometry rule under test: `eps_step_x` by triangle centroid.")
    if convergence_path is not None and convergence_path.exists():
        lines.append(f"- Joined convergence CSV: `{_rel_from_root(convergence_path)}`")
    lines.append("")
    lines.append("This diagnostic measures how the centroid-based material assignment represents the target vertical interface on each mesh.")
    lines.append("")

    for spec in PRIORITY_SPECS:
        point_rows = rows_by_point[spec.priority_point]
        lines.append(f"## {spec.priority_point} `{spec.mode}`")
        lines.append("")
        lines.append("| Mesh | Target d/a | Effective d/a | |delta d/a| | Symmetric diff/domain | Cut tris | Cut cells | Aligned | Cut-cell mode | Err. ana % | Err. ref % |")
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: |")
        for row in point_rows:
            lines.append(
                "| "
                f"{row['mesh_label']} | "
                f"{row['d_over_a_target']} | "
                f"{row['effective_fill_fraction_centroid']} | "
                f"{row['fill_fraction_abs_error']} | "
                f"{row['symmetric_diff_fraction_domain']} | "
                f"{row['cut_triangle_count']} | "
                f"{row['cut_cell_count']} | "
                f"{'yes' if row['interface_aligned_to_mesh_edge'] == '1' else 'no'} | "
                f"{row['cut_column_mode']} | "
                f"{row['err_primary_pct']} | "
                f"{row['err_secondary_pct']} |"
            )
        if point_rows:
            first = point_rows[0]
            last = point_rows[-1]
            lines.append("")
            lines.append(
                f"- Geometric |delta d/a| trend: `{first['fill_fraction_abs_error']}` -> `{last['fill_fraction_abs_error']}`."
            )
            lines.append(
                f"- Symmetric difference/domain trend: `{first['symmetric_diff_fraction_domain']}` -> `{last['symmetric_diff_fraction_domain']}`."
            )
            if first.get("err_primary_pct") and last.get("err_primary_pct"):
                lines.append(
                    f"- Primary error trend from convergence diagnostic: `{first['err_primary_pct']}` -> `{last['err_primary_pct']}`."
                )
            if first.get("err_secondary_pct") and last.get("err_secondary_pct"):
                lines.append(
                    f"- Secondary error trend from convergence diagnostic: `{first['err_secondary_pct']}` -> `{last['err_secondary_pct']}`."
                )
        lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append(f"- CSV: `[{csv_path.name}]({csv_path.name})`")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Diagnose centroid-based interface staircasing for Figure 13 / Table 10."
    )
    ap.add_argument(
        "--out-csv",
        type=Path,
        default=DEFAULT_OUT_CSV,
        help="Path to the consolidated staircasing CSV.",
    )
    ap.add_argument(
        "--out-md",
        type=Path,
        default=DEFAULT_OUT_MD,
        help="Path to the short Markdown summary.",
    )
    ap.add_argument(
        "--convergence-csv",
        type=Path,
        default=DEFAULT_CONVERGENCE_CSV,
        help="Optional diag_224_table10_convergence.csv to join numeric errors.",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)
    convergence_csv = _resolve(args.convergence_csv) if args.convergence_csv else None

    convergence_rows = _load_convergence_rows(convergence_csv)
    rows = build_rows(convergence_rows)
    rows.sort(key=lambda row: (row["priority_point"], int(row["nx"]), int(row["ny"])))

    write_csv(out_csv, rows)
    write_md(out_md, rows, csv_path=out_csv, convergence_path=convergence_csv)

    print(f"[diag] wrote {_rel_from_root(out_csv)}")
    print(f"[diag] wrote {_rel_from_root(out_md)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
