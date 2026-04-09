#!/usr/bin/env python3
"""
Run a focused convergence diagnostic for Section 2.2.4 Figure 13 / Table 10.

Default outputs:
- out/validation/diag_224_table10_convergence.csv
- out/validation/DIAG_224_TABLE10_CONVERGENCE.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD_DIR = ROOT / "build"
DEFAULT_OUT_ROOT = ROOT / "out"
SECTION = "2.2.4"
CASE = "Figure13_Table10"
EXECUTABLE = "helmvec3_fig13_rect"
DEFAULT_MESHES = "10x5,16x8,20x10"
DEFAULT_PREVIEW_D_OVER_A = "0.20"


@dataclass(frozen=True)
class PriorityPoint:
    key: str
    d_over_a: float
    br_over_lambda0: float


PRIORITY_POINTS = (
    PriorityPoint("P1", 0.167, 0.5),
    PriorityPoint("P2", 0.286, 0.5),
    PriorityPoint("P3", 0.5, 0.4),
)


@dataclass(frozen=True)
class MeshSpec:
    nx: int
    ny: int

    @property
    def label(self) -> str:
        return f"nx{self.nx}_ny{self.ny}"

    @property
    def cells(self) -> int:
        return 2 * self.nx * self.ny


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _to_float(raw: str, *, field: str, path: Path) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError) as exc:
        raise SystemExit(
            f"Invalid numeric value for {field!r} in {_rel_from_root(path)}: {raw!r}"
        ) from exc


def _to_int(raw: str, *, field: str, path: Path) -> int:
    try:
        return int(raw)
    except (TypeError, ValueError) as exc:
        raise SystemExit(
            f"Invalid integer value for {field!r} in {_rel_from_root(path)}: {raw!r}"
        ) from exc


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _format_mode(point: PriorityPoint) -> str:
    return f"d/a={point.d_over_a},br/lambda0={point.br_over_lambda0}"


def _parse_meshes(raw: str) -> list[MeshSpec]:
    meshes: list[MeshSpec] = []
    for token in raw.split(","):
        text = token.strip().lower()
        if not text:
            continue
        match = re.fullmatch(r"(\d+)x(\d+)", text)
        if not match:
            raise SystemExit(
                f"Invalid mesh token {token!r}. Expected comma-separated entries like 10x5,20x10,30x15."
            )
        nx = int(match.group(1))
        ny = int(match.group(2))
        if nx <= 0 or ny <= 0:
            raise SystemExit(f"Mesh dimensions must be positive: {token!r}")
        meshes.append(MeshSpec(nx=nx, ny=ny))
    if not meshes:
        raise SystemExit("At least one mesh must be provided via --meshes.")
    return meshes


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    return rows


def _read_run_timing(path: Path) -> dict[str, str]:
    rows = _read_csv_rows(path)
    if len(rows) != 1:
        raise SystemExit(
            f"Expected exactly one timing row in {_rel_from_root(path)}, found {len(rows)}."
        )
    return rows[0]


def _extract_priority_rows(
    rows: list[dict[str, str]],
    *,
    path: Path,
) -> dict[tuple[float, float], dict[str, str]]:
    grouped: dict[tuple[float, float], list[dict[str, str]]] = {}
    for row in rows:
        key = _point_key(
            _to_float(row["d_over_a"], field="d_over_a", path=path),
            _to_float(row["br_over_lambda0"], field="br_over_lambda0", path=path),
        )
        grouped.setdefault(key, []).append(row)

    selected: dict[tuple[float, float], dict[str, str]] = {}
    for point in PRIORITY_POINTS:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        candidates = grouped.get(key, [])
        if not candidates:
            raise SystemExit(
                f"Missing priority point {_format_mode(point)} in {_rel_from_root(path)}."
            )
        if len(candidates) != 1:
            raise SystemExit(
                f"Expected one row for {_format_mode(point)} in {_rel_from_root(path)}, found {len(candidates)}."
            )
        selected[key] = candidates[0]
    return selected


def _parse_candidate_log(path: Path) -> dict[tuple[float, float], list[float]]:
    if not path.exists():
        raise FileNotFoundError(path)

    point_candidates: dict[tuple[float, float], list[float]] = {}
    current_d: float | None = None
    block_re = re.compile(r"\[debug\]\s+Figure13 block d/a=([0-9.]+)")
    cand_re = re.compile(r"\[debug\]\s+s=([0-9.]+)\s+cands:\s*(.*)$")

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            block_match = block_re.search(line)
            if block_match:
                current_d = float(block_match.group(1))
                continue

            cand_match = cand_re.search(line)
            if cand_match and current_d is not None:
                br_over_lambda0 = float(cand_match.group(1))
                candidate_text = cand_match.group(2).strip()
                values = [float(token) for token in candidate_text.split()] if candidate_text else []
                point_candidates[_point_key(current_d, br_over_lambda0)] = values

    for point in PRIORITY_POINTS:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        if key not in point_candidates:
            raise SystemExit(
                f"Missing debug candidate list for {_format_mode(point)} in {_rel_from_root(path)}. "
                "Re-run the diagnostic with --debug-candidates enabled."
            )
    return point_candidates


def _closest(values: list[float], target: float) -> tuple[float | None, float | None]:
    if not values:
        return (None, None)
    best = min(values, key=lambda value: abs(value - target))
    return (best, abs(best - target))


def _run_case(
    mesh: MeshSpec,
    *,
    build_dir: Path,
    out_root: Path,
    backend: str,
    preview_d_over_a: str,
) -> dict[str, Path]:
    executable = build_dir / EXECUTABLE
    if not executable.exists():
        raise SystemExit(
            f"Missing executable: {_rel_from_root(executable)}. Build the project first."
        )

    case_root = out_root / "diag_224_table10_convergence" / mesh.label
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(case_root)
    cmd = [
        str(executable),
        "--d-over-a-preview",
        preview_d_over_a,
        "--nx",
        str(mesh.nx),
        "--ny",
        str(mesh.ny),
        "--backend",
        backend,
        "--debug-candidates",
    ]
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        stdout_tail = proc.stdout[-2000:].strip()
        stderr_tail = proc.stderr[-2000:].strip()
        raise SystemExit(
            "Diagnostic run failed for "
            f"{mesh.label} with exit code {proc.returncode}.\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout tail:\n{stdout_tail}\n"
            f"stderr tail:\n{stderr_tail}"
        )

    fig13_root = case_root / "helmvec3" / "fig13_rect"
    paths = {
        "root": fig13_root,
        "table10_csv": fig13_root / "csv" / "helmvec3_fig13_rect_table10.csv",
        "run_log": fig13_root / "run.log",
        "run_timing_csv": fig13_root / "run_timing.csv",
    }
    for label, path in paths.items():
        if label == "root":
            continue
        if not path.exists():
            raise SystemExit(
                f"Expected diagnostic artifact not found for {mesh.label}: {_rel_from_root(path)}"
            )
    return paths


def build_rows(
    *,
    meshes: list[MeshSpec],
    build_dir: Path,
    out_root: Path,
    backend: str,
    preview_d_over_a: str,
) -> list[dict[str, str]]:
    output_rows: list[dict[str, str]] = []

    for mesh in meshes:
        paths = _run_case(
            mesh,
            build_dir=build_dir,
            out_root=out_root,
            backend=backend,
            preview_d_over_a=preview_d_over_a,
        )
        table10_rows = _read_csv_rows(paths["table10_csv"])
        timing_row = _read_run_timing(paths["run_timing_csv"])
        priority_rows = _extract_priority_rows(table10_rows, path=paths["table10_csv"])
        candidate_lists = _parse_candidate_log(paths["run_log"])

        for point in PRIORITY_POINTS:
            key = _point_key(point.d_over_a, point.br_over_lambda0)
            row = priority_rows[key]
            candidates = candidate_lists[key]

            analytic = _to_float(
                row["beta_over_k0_analytic"],
                field="beta_over_k0_analytic",
                path=paths["table10_csv"],
            )
            helmvec3 = _to_float(
                row["beta_over_k0_helmvec3"],
                field="beta_over_k0_helmvec3",
                path=paths["table10_csv"],
            )
            fem = _to_float(
                row["beta_over_k0_fem_matched"],
                field="beta_over_k0_fem_matched",
                path=paths["table10_csv"],
            )
            closest_analytic, closest_analytic_abs_err = _closest(candidates, analytic)
            closest_helmvec3, closest_helmvec3_abs_err = _closest(candidates, helmvec3)
            selected_delta = min((abs(candidate - fem) for candidate in candidates), default=math.nan)

            output_rows.append(
                {
                    "backend": backend,
                    "section": SECTION,
                    "case": CASE,
                    "priority_point": point.key,
                    "mode": _format_mode(point),
                    "mesh_label": mesh.label,
                    "nx": str(mesh.nx),
                    "ny": str(mesh.ny),
                    "mesh_cells_expected": str(mesh.cells),
                    "mesh_nodes": timing_row.get("mesh_nodes", ""),
                    "mesh_tris": timing_row.get("mesh_tris", ""),
                    "preview_d_over_a": preview_d_over_a,
                    "fem": row["beta_over_k0_fem_matched"],
                    "ref_primary": row["beta_over_k0_analytic"],
                    "ref_secondary": row["beta_over_k0_helmvec3"],
                    "err_primary_pct": row["error_percent_analytic"],
                    "err_secondary_pct": row["error_percent_helmvec3"],
                    "selected_rank": row["selected_candidate_rank"],
                    "selected_eig_index": row["selected_eig_index"],
                    "ez_ratio": row["ez_ratio"],
                    "match_status": row["match_status"],
                    "field_status": row["field_status"],
                    "candidate_count": str(len(candidates)),
                    "candidate_min": f"{min(candidates):.16g}" if candidates else "",
                    "candidate_max": f"{max(candidates):.16g}" if candidates else "",
                    "closest_candidate_to_analytic": (
                        f"{closest_analytic:.16g}" if closest_analytic is not None else ""
                    ),
                    "closest_candidate_err_abs_analytic": (
                        f"{closest_analytic_abs_err:.16g}"
                        if closest_analytic_abs_err is not None
                        else ""
                    ),
                    "closest_candidate_to_helmvec3": (
                        f"{closest_helmvec3:.16g}" if closest_helmvec3 is not None else ""
                    ),
                    "closest_candidate_err_abs_helmvec3": (
                        f"{closest_helmvec3_abs_err:.16g}"
                        if closest_helmvec3_abs_err is not None
                        else ""
                    ),
                    "selected_candidate_abs_delta_in_list": (
                        f"{selected_delta:.16g}" if math.isfinite(selected_delta) else ""
                    ),
                    "candidate_beta_ratios": "|".join(f"{value:.12g}" for value in candidates),
                    "et_fields_csv_file": row["et_fields_csv_file"],
                    "ez_fields_csv_file": row["ez_fields_csv_file"],
                    "et_vtk_file": row["et_vtk_file"],
                    "ez_vtk_file": row["ez_vtk_file"],
                    "table10_csv": _rel_from_root(paths["table10_csv"]),
                    "run_log": _rel_from_root(paths["run_log"]),
                    "run_timing_csv": _rel_from_root(paths["run_timing_csv"]),
                }
            )

    return output_rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "backend",
        "section",
        "case",
        "priority_point",
        "mode",
        "mesh_label",
        "nx",
        "ny",
        "mesh_cells_expected",
        "mesh_nodes",
        "mesh_tris",
        "preview_d_over_a",
        "fem",
        "ref_primary",
        "ref_secondary",
        "err_primary_pct",
        "err_secondary_pct",
        "selected_rank",
        "selected_eig_index",
        "ez_ratio",
        "match_status",
        "field_status",
        "candidate_count",
        "candidate_min",
        "candidate_max",
        "closest_candidate_to_analytic",
        "closest_candidate_err_abs_analytic",
        "closest_candidate_to_helmvec3",
        "closest_candidate_err_abs_helmvec3",
        "selected_candidate_abs_delta_in_list",
        "candidate_beta_ratios",
        "et_fields_csv_file",
        "ez_fields_csv_file",
        "et_vtk_file",
        "ez_vtk_file",
        "table10_csv",
        "run_log",
        "run_timing_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_md(
    path: Path,
    rows: list[dict[str, str]],
    *,
    meshes: list[MeshSpec],
    csv_path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows_by_point: dict[str, list[dict[str, str]]] = {point.key: [] for point in PRIORITY_POINTS}
    for row in rows:
        rows_by_point[row["priority_point"]].append(row)
    for point_rows in rows_by_point.values():
        point_rows.sort(key=lambda row: (int(row["nx"]), int(row["ny"])))

    lines: list[str] = []
    lines.append("# Figure 13 / Table 10 Convergence Diagnostic")
    lines.append("")
    lines.append(f"- Generated at: `{dt.datetime.now().isoformat(timespec='seconds')}`")
    lines.append(f"- Section: `{SECTION}`")
    lines.append(f"- Case: `{CASE}`")
    lines.append(f"- Backend: `{rows[0]['backend'] if rows else ''}`")
    lines.append(f"- Meshes: `{', '.join(mesh.label for mesh in meshes)}`")
    lines.append("")
    lines.append("This diagnostic focuses only on the three priority points flagged by the scientific gate.")
    lines.append("")

    for point in PRIORITY_POINTS:
        point_rows = rows_by_point[point.key]
        lines.append(f"## {point.key} `{_format_mode(point)}`")
        lines.append("")
        lines.append("| Mesh | FEM | Analytic | HELMVEC3 | Err. ana % | Err. ref % | Closest ana cand | Closest ana | Closest ref cand | Closest ref | Rank | Eig | Ez ratio | #cands |")
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for row in point_rows:
            lines.append(
                "| "
                f"{row['mesh_label']} | "
                f"{row['fem']} | "
                f"{row['ref_primary']} | "
                f"{row['ref_secondary']} | "
                f"{row['err_primary_pct']} | "
                f"{row['err_secondary_pct']} | "
                f"{row['closest_candidate_to_analytic']} | "
                f"{row['closest_candidate_err_abs_analytic']} | "
                f"{row['closest_candidate_to_helmvec3']} | "
                f"{row['closest_candidate_err_abs_helmvec3']} | "
                f"{row['selected_rank']} | "
                f"{row['selected_eig_index']} | "
                f"{row['ez_ratio']} | "
                f"{row['candidate_count']} |"
            )
        if point_rows:
            first = point_rows[0]
            last = point_rows[-1]
            lines.append("")
            lines.append(
                f"- Selected FEM trend: `{first['fem']}` -> `{last['fem']}`."
            )
            lines.append(
                f"- Closest analytic-candidate distance: `{first['closest_candidate_err_abs_analytic']}` -> `{last['closest_candidate_err_abs_analytic']}`."
            )
            lines.append(
                f"- Closest HELMVEC3-candidate distance: `{first['closest_candidate_err_abs_helmvec3']}` -> `{last['closest_candidate_err_abs_helmvec3']}`."
            )
            lines.append(
                f"- Primary error trend: `{first['err_primary_pct']}` -> `{last['err_primary_pct']}`."
            )
            lines.append(
                f"- Secondary error trend: `{first['err_secondary_pct']}` -> `{last['err_secondary_pct']}`."
            )
        lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append(f"- CSV: `[{csv_path.name}]({csv_path.name})`")
    lines.append("")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Run a focused mesh-convergence diagnostic for Figure 13 / Table 10."
    )
    ap.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD_DIR,
        help="Directory containing helmvec3_fig13_rect (default: build).",
    )
    ap.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_OUT_ROOT,
        help="Base output directory used for isolated diagnostic runs (default: out).",
    )
    ap.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/diag_224_table10_convergence.csv"),
        help="Path to the consolidated diagnostic CSV.",
    )
    ap.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/DIAG_224_TABLE10_CONVERGENCE.md"),
        help="Path to the optional Markdown summary.",
    )
    ap.add_argument(
        "--backend",
        default="closed-form",
        help="Assembly backend passed to helmvec3_fig13_rect (default: closed-form).",
    )
    ap.add_argument(
        "--d-over-a-preview",
        default=DEFAULT_PREVIEW_D_OVER_A,
        help="Preview d/a forwarded to helmvec3_fig13_rect (default: 0.20).",
    )
    ap.add_argument(
        "--meshes",
        default=DEFAULT_MESHES,
        help="Comma-separated mesh grid, preserving aspect ratio when desired (default: 10x5,16x8,20x10).",
    )
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    build_dir = _resolve(args.build_dir)
    out_root = _resolve(args.out_root)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)
    meshes = _parse_meshes(args.meshes)

    rows = build_rows(
        meshes=meshes,
        build_dir=build_dir,
        out_root=out_root,
        backend=args.backend,
        preview_d_over_a=args.d_over_a_preview,
    )
    expected_rows = len(meshes) * len(PRIORITY_POINTS)
    if len(rows) != expected_rows:
        raise SystemExit(
            f"Unexpected diagnostic row count: got {len(rows)}, expected {expected_rows}."
        )

    rows.sort(key=lambda row: (row["priority_point"], int(row["nx"]), int(row["ny"])))
    write_csv(out_csv, rows)
    write_md(out_md, rows, meshes=meshes, csv_path=out_csv)

    print(f"[diag] wrote {_rel_from_root(out_csv)}")
    print(f"[diag] wrote {_rel_from_root(out_md)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
