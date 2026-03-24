#!/usr/bin/env python3
"""
Varredura completa de malhas para todos os casos 2D/3D, em ambos os backends.

Objetivos:
- rodar todos os executaveis relevantes em 10 niveis de discretizacao;
- separar saidas por backend e por nivel;
- preservar stdout/log e tempos por execucao;
- gerar CSVs agregados e um README resumo com erros e tempos.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import statistics
import subprocess
import textwrap
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from benchmark_backends import parse_timing
from validate_2d_22 import (
    parse_first_kc_block,
    parse_helmvec2_table,
    parse_helmvec3_table10,
    parse_helmvec3_table9,
    parse_mixed_rect_table,
)
from validate_3d_31 import parse_mode_table


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD = ROOT / "build"
DEFAULT_OUT = ROOT / "out" / "sweeps" / "full_mesh_compare"

TABLE_ROW_RE = re.compile(
    r"^\s*(\d+)\s+\((\d+),(\d+)\)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE-]+)\s+([0-9.+-eE]+)\s*$"
)


RECT_LEVELS = [
    (6, 6),
    (8, 8),
    (10, 10),
    (12, 12),
    (14, 14),
    (16, 16),
    (18, 18),
    (20, 20),
    (24, 24),
    (28, 28),
]

POLAR_LEVELS = [
    (4, 24),
    (5, 28),
    (6, 32),
    (7, 36),
    (8, 40),
    (10, 48),
    (12, 56),
    (14, 64),
    (16, 72),
    (20, 96),
]

RECT_RATIO_LEVELS = [
    (4, 2),
    (6, 3),
    (8, 4),
    (10, 5),
    (12, 6),
    (14, 7),
    (16, 8),
    (18, 9),
    (20, 10),
    (24, 12),
]

HV2_LEVELS = [
    (4, 4),
    (5, 5),
    (6, 6),
    (7, 7),
    (8, 8),
    (10, 10),
    (12, 12),
    (14, 14),
    (16, 16),
    (20, 20),
]

GRID3D_LEVELS = {
    "air": [
        (3, 2, 2), (4, 2, 2), (4, 3, 3), (5, 3, 3), (6, 3, 3),
        (6, 4, 4), (7, 4, 4), (8, 4, 4), (8, 5, 5), (9, 5, 5),
    ],
    "half": [
        (3, 3, 2), (4, 4, 3), (4, 4, 4), (5, 5, 4), (5, 5, 5),
        (6, 5, 4), (6, 6, 5), (7, 6, 5), (7, 7, 6), (8, 7, 6),
    ],
    "cyl": [
        (4, 4, 2), (5, 5, 2), (6, 6, 3), (6, 6, 4), (7, 7, 4),
        (7, 7, 5), (8, 8, 5), (8, 8, 6), (9, 9, 6), (10, 10, 7),
    ],
    "sphere": [
        (4, 3, 3), (5, 4, 4), (5, 5, 4), (6, 5, 5), (6, 6, 5),
        (7, 6, 6), (7, 7, 6), (8, 7, 7), (8, 8, 7), (9, 8, 8),
    ],
}


@dataclass(frozen=True)
class RunResult:
    case_id: str
    section: str
    backend: str
    level: int
    resolution: str
    command: str
    stdout_rel: str
    stderr_rel: str
    wall_ms: float
    assembly_ms: float | None
    solve_ms: float | None
    post_ms: float | None
    total_ms: float | None
    returncode: int


def run_checked(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    stdout_path: Path,
    stderr_path: Path,
) -> tuple[str, str, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    wall_ms = (time.perf_counter() - t0) * 1000.0
    stdout_path.write_text(proc.stdout, encoding="utf-8")
    stderr_path.write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"Falha ao executar: {' '.join(cmd)}\n"
            f"stdout: {stdout_path}\n"
            f"stderr: {stderr_path}"
        )
    return proc.stdout, proc.stderr, wall_ms


def parse_rect_or_polar_table(
    text: str,
    marker: str,
    *,
    section: str,
    case: str,
    geometry: str,
    formulation: str,
    polarization: str,
    index2_name: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    start = text.find(marker)
    if start < 0:
        return rows
    sub = text[start:].splitlines()
    started = False
    for ln in sub:
        if "kc_ana" in ln and "kc_fem" in ln:
            started = True
            continue
        if not started:
            continue
        m = TABLE_ROW_RE.match(ln)
        if not m:
            if rows and not ln.strip():
                break
            continue
        rows.append(
            {
                "section": section,
                "case": case,
                "geometry": geometry,
                "formulation": formulation,
                "polarization": polarization,
                "mode_rank": int(m.group(1)),
                "m": int(m.group(2)),
                "index2_name": index2_name,
                "index2": int(m.group(3)),
                "fem": float(m.group(5)),
                "ref_primary": float(m.group(4)),
                "err_primary_pct": float(m.group(6)),
                "rho_abs": float(m.group(7)),
            }
        )
    return rows


def level_resolution_summary(level: int) -> dict[str, str]:
    rect = RECT_LEVELS[level - 1]
    polar = POLAR_LEVELS[level - 1]
    mixed_rect = RECT_RATIO_LEVELS[level - 1]
    hv2 = HV2_LEVELS[level - 1]
    hv3 = RECT_RATIO_LEVELS[level - 1]
    air = GRID3D_LEVELS["air"][level - 1]
    half = GRID3D_LEVELS["half"][level - 1]
    cyl = GRID3D_LEVELS["cyl"][level - 1]
    sphere = GRID3D_LEVELS["sphere"][level - 1]
    return {
        "helm10_rect_edge_rect": f"{rect[0]}x{rect[1]}",
        "polar_2d": f"nr={polar[0]},nt={polar[1]}",
        "mixed_rect_hv3": f"{mixed_rect[0]}x{mixed_rect[1]}",
        "helmvec2_rect": f"{hv2[0]}x{hv2[1]}",
        "3d_air": f"{air[0]}x{air[1]}x{air[2]}",
        "3d_half": f"{half[0]}x{half[1]}x{half[2]}",
        "3d_cyl": f"{cyl[0]}x{cyl[1]}x{cyl[2]}",
        "3d_sphere": f"{sphere[0]}x{sphere[1]}x{sphere[2]}",
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def mean_abs(values: Iterable[float]) -> float:
    vals = [abs(v) for v in values]
    return statistics.fmean(vals) if vals else math.nan


def fmt(v: float | None) -> str:
    if v is None or math.isnan(v):
        return "-"
    return f"{v:.3f}"


def build_project(build_dir: Path, jobs: int) -> None:
    subprocess.run(
        ["cmake", "-S", str(ROOT), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"],
        check=True,
    )
    subprocess.run(["cmake", "--build", str(build_dir), "-j", str(jobs)], check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description="Varredura completa de malhas para gauss e closed-form.")
    ap.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD)
    ap.add_argument("--out-root", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--mode-export", type=int, default=8)
    ap.add_argument("--levels", type=int, default=10, help="Numero de niveis de malha (maximo 10).")
    ap.add_argument("--clean", action="store_true", help="Remove o diretorio de saida antes da varredura.")
    ap.add_argument("--skip-build", action="store_true")
    ap.add_argument("--skip-images", action="store_true")
    ap.add_argument(
        "--backends",
        default="gauss,closed-form",
        help="Lista separada por virgula. Ex.: gauss,closed-form",
    )
    ap.add_argument(
        "--cases",
        default="all",
        help="Lista separada por virgula de ids de execucao ou 'all'. Ex.: helm10_rect,edge_rect,2d,3d",
    )
    args = ap.parse_args()

    build_dir = args.build_dir.resolve()
    out_root = args.out_root.resolve()
    if args.clean and out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_build:
        build_project(build_dir, args.jobs)

    levels = max(1, min(args.levels, 10))
    backends = [b.strip() for b in args.backends.split(",") if b.strip()]
    valid_backends = {"gauss", "closed-form"}
    for backend in backends:
        if backend not in valid_backends:
            raise SystemExit(f"Backend invalido: {backend}")

    case_tokens = {c.strip().lower() for c in args.cases.split(",") if c.strip()}
    run_all_cases = "all" in case_tokens

    def wants(case_id: str, family: str) -> bool:
        if run_all_cases:
            return True
        cid = case_id.lower()
        fam = family.lower()
        if cid in case_tokens or fam in case_tokens:
            return True
        return any(cid.startswith(token + "_") or token.startswith(cid + "_") for token in case_tokens)

    run_rows: list[dict[str, object]] = []
    mode_rows: list[dict[str, object]] = []
    level_rows: list[dict[str, object]] = []

    for level in range(1, levels + 1):
        rect = RECT_LEVELS[level - 1]
        polar = POLAR_LEVELS[level - 1]
        mixed_rect = RECT_RATIO_LEVELS[level - 1]
        hv2 = HV2_LEVELS[level - 1]
        hv3 = RECT_RATIO_LEVELS[level - 1]
        grids3d = {k: GRID3D_LEVELS[k][level - 1] for k in GRID3D_LEVELS}

        manifest = level_resolution_summary(level)
        for backend in backends:
            level_root = out_root / backend.replace("-", "_") / f"level_{level:02d}"
            logs_dir = level_root / "logs"
            logs_dir.mkdir(parents=True, exist_ok=True)

            env = os.environ.copy()
            env["TP3485_OUT_DIR"] = str(level_root)

            def run_exe(case_id: str, section: str, resolution: str, cmd: list[str]) -> str:
                exe_name = Path(cmd[0]).name
                safe_case = case_id.replace("/", "_")
                stdout_path = logs_dir / f"{safe_case}__{exe_name}.stdout.txt"
                stderr_path = logs_dir / f"{safe_case}__{exe_name}.stderr.txt"
                stdout, _stderr, wall_ms = run_checked(cmd, cwd=build_dir, env=env, stdout_path=stdout_path, stderr_path=stderr_path)
                timing = parse_timing(stdout)
                run_rows.append(
                    {
                        "case_id": case_id,
                        "section": section,
                        "backend": backend,
                        "level": level,
                        "resolution": resolution,
                        "command": " ".join(cmd),
                        "stdout_rel": str(stdout_path.relative_to(out_root)),
                        "stderr_rel": str(stderr_path.relative_to(out_root)),
                        "wall_ms": f"{wall_ms:.3f}",
                        "assembly_ms": "" if math.isnan(timing["assembly_ms"]) else f"{timing['assembly_ms']:.3f}",
                        "solve_ms": "" if math.isnan(timing["solve_ms"]) else f"{timing['solve_ms']:.3f}",
                        "post_ms": "" if math.isnan(timing["post_ms"]) else f"{timing['post_ms']:.3f}",
                        "total_ms": "" if math.isnan(timing["total_ms"]) else f"{timing['total_ms']:.3f}",
                        "returncode": 0,
                    }
                )
                return stdout

            # 2.1 scalar
            if wants("helm10_rect", "2d"):
                stdout = run_exe("helm10_rect", "2.1", f"{rect[0]}x{rect[1]}",
                                 [str(build_dir / "helm10_rect"), str(rect[0]), str(rect[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{rect[0]}x{rect[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 1 (TE)", section="2.1", case="rect_scalar_te", geometry="rect", formulation="scalar", polarization="TE", index2_name="n")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{rect[0]}x{rect[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 1 (TM)", section="2.1", case="rect_scalar_tm", geometry="rect", formulation="scalar", polarization="TM", index2_name="n")
                ]

            if wants("helm10_circle", "2d"):
                stdout = run_exe("helm10_circle", "2.1", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "helm10_circle"), str(polar[0]), str(polar[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 2 (TE)", section="2.1", case="circle_scalar_te", geometry="circle", formulation="scalar", polarization="TE", index2_name="p")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 2 (TM)", section="2.1", case="circle_scalar_tm", geometry="circle", formulation="scalar", polarization="TM", index2_name="p")
                ]

            if wants("helm10_coax", "2d"):
                stdout = run_exe("helm10_coax", "2.1", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "helm10_coax"), str(polar[0]), str(polar[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 3 (TE coax)", section="2.1", case="coax_scalar_te", geometry="coax", formulation="scalar", polarization="TE", index2_name="p")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela 3 (TM coax)", section="2.1", case="coax_scalar_tm", geometry="coax", formulation="scalar", polarization="TM", index2_name="p")
                ]

            # 2.2.1 edge
            if wants("edge_rect", "2d"):
                stdout = run_exe("edge_rect", "2.2.1", f"{rect[0]}x{rect[1]}",
                                 [str(build_dir / "edge_rect"), str(rect[0]), str(rect[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{rect[0]}x{rect[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TE", section="2.2.1", case="rect_edge_te", geometry="rect", formulation="edge", polarization="TE", index2_name="n")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{rect[0]}x{rect[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TM", section="2.2.1", case="rect_edge_tm", geometry="rect", formulation="edge", polarization="TM", index2_name="n")
                ]

            if wants("edge_circle", "2d"):
                stdout = run_exe("edge_circle", "2.2.1", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "edge_circle"), str(polar[0]), str(polar[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TE", section="2.2.1", case="circle_edge_te", geometry="circle", formulation="edge", polarization="TE", index2_name="p")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TM", section="2.2.1", case="circle_edge_tm", geometry="circle", formulation="edge", polarization="TM", index2_name="p")
                ]

            if wants("edge_coax", "2d"):
                stdout = run_exe("edge_coax", "2.2.1", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "edge_coax"), str(polar[0]), str(polar[1]), str(args.mode_export), "--backend", backend])
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TE", section="2.2.1", case="coax_edge_te", geometry="coax", formulation="edge", polarization="TE", index2_name="p")
                ]
                mode_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"nr={polar[0]},nt={polar[1]}"} for row in
                    parse_rect_or_polar_table(stdout, "Tabela (TM", section="2.2.1", case="coax_edge_tm", geometry="coax", formulation="edge", polarization="TM", index2_name="p")
                ]

            # 2.2.2 mixed
            mixed_2d_rows: list[dict[str, object]] = []
            if wants("mixed_rect", "2d"):
                stdout = run_exe("mixed_rect", "2.2.2", f"{mixed_rect[0]}x{mixed_rect[1]}",
                                 [str(build_dir / "mixed_rect"), str(mixed_rect[0]), str(mixed_rect[1]), "--backend", backend])
                mixed_2d_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{mixed_rect[0]}x{mixed_rect[1]}"} for row in
                    parse_mixed_rect_table(stdout, "[E-formulation] TE cutoffs (edge block)", "mixed_rect_E_TE_table")
                ]
                mixed_2d_rows += [
                    {**row, "backend": backend, "level": level, "resolution": f"{mixed_rect[0]}x{mixed_rect[1]}"} for row in
                    parse_mixed_rect_table(stdout, "[E-formulation] TM cutoffs (scalar block)", "mixed_rect_E_TM_table")
                ]
            if wants("mixed_circle", "2d"):
                stdout = run_exe("mixed_circle", "2.2.2", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "mixed_circle"), str(polar[0]), str(polar[1]), "--backend", backend])
                snapshots = [
                    ("mixed_circle_TE_edge", "TE (edge block), first 8 kc:"),
                    ("mixed_circle_TM_scalar", "TM (scalar block), first 8 kc:"),
                ]
                for case_name, label in snapshots:
                    vals = parse_first_kc_block(stdout, label)
                    for i, val in enumerate(vals, start=1):
                        mixed_2d_rows.append({
                            "section": "2.2.2",
                            "case": case_name,
                            "mode": i,
                            "fem": val,
                            "ref_primary": "",
                            "ref_secondary": "",
                            "err_primary_pct": "",
                            "err_secondary_pct": "",
                            "backend": backend,
                            "level": level,
                            "resolution": f"nr={polar[0]},nt={polar[1]}",
                        })
            if wants("mixed_coax", "2d"):
                stdout = run_exe("mixed_coax", "2.2.2", f"nr={polar[0]},nt={polar[1]}",
                                 [str(build_dir / "mixed_coax"), str(polar[0]), str(polar[1]), "--backend", backend])
                snapshots = [
                    ("mixed_coax_TE_edge", "TE (edge block), first 8 kc:"),
                    ("mixed_coax_TM_scalar", "TM (scalar block), first 8 kc:"),
                ]
                for case_name, label in snapshots:
                    vals = parse_first_kc_block(stdout, label)
                    for i, val in enumerate(vals, start=1):
                        mixed_2d_rows.append({
                            "section": "2.2.2",
                            "case": case_name,
                            "mode": i,
                            "fem": val,
                            "ref_primary": "",
                            "ref_secondary": "",
                            "err_primary_pct": "",
                            "err_secondary_pct": "",
                            "backend": backend,
                            "level": level,
                            "resolution": f"nr={polar[0]},nt={polar[1]}",
                        })

            # 2.2.3
            hv2_rows: list[dict[str, object]] = []
            if wants("helmvec2_rect", "2d"):
                stdout = run_exe("helmvec2_rect", "2.2.3", f"{hv2[0]}x{hv2[1]}",
                                 [str(build_dir / "helmvec2_rect"), "10", str(hv2[0]), str(hv2[1]), "--backend", backend])
                hv2_rows += [{**row, "backend": backend, "level": level, "resolution": f"{hv2[0]}x{hv2[1]}"} for row in parse_helmvec2_table(stdout)]

            # 2.2.4
            hv3_rows: list[dict[str, object]] = []
            if wants("helmvec3_rect", "2d"):
                stdout = run_exe("helmvec3_rect", "2.2.4", f"{hv3[0]}x{hv3[1]}",
                                 [str(build_dir / "helmvec3_rect"), "0.20", str(hv3[0]), str(hv3[1]), "--backend", backend])
                hv3_rows += [{**row, "backend": backend, "level": level, "resolution": f"{hv3[0]}x{hv3[1]}"} for row in parse_helmvec3_table9(stdout)]
                hv3_rows += [{**row, "backend": backend, "level": level, "resolution": f"{hv3[0]}x{hv3[1]}"} for row in parse_helmvec3_table10(stdout)]

            if mixed_2d_rows or hv2_rows or hv3_rows:
                csv_22 = level_root / "validation" / "validation_2d_22.csv"
                rows_22 = mixed_2d_rows + hv2_rows + hv3_rows
                write_csv(
                    csv_22,
                    rows_22,
                    ["backend", "level", "resolution", "section", "case", "mode", "fem", "ref_primary", "ref_secondary", "err_primary_pct", "err_secondary_pct"],
                )
                if not args.skip_images:
                    subprocess.run(
                        [
                            "python3",
                            str(ROOT / "scripts" / "plot_validation_2d_22.py"),
                            "--in-csv",
                            str(csv_22),
                            "--out-dir",
                            str(level_root / "img_all" / "validation_2d_22"),
                            "--backend",
                            backend,
                        ],
                        check=True,
                    )

            # Generate 2D images and mode summary for 2.1 + 2.2.1
            if not args.skip_images and any(
                wants(case_name, "2d")
                for case_name in (
                    "helm10_rect",
                    "helm10_circle",
                    "helm10_coax",
                    "edge_rect",
                    "edge_circle",
                    "edge_coax",
                )
            ):
                subprocess.run(
                    [
                        "python3",
                        str(ROOT / "scripts" / "plot_vtk_quiver.py"),
                        "--all-img",
                        "--vtk-root",
                        str(level_root / "2d"),
                        "--out-dir",
                        str(level_root / "img_all"),
                        "--csv",
                        str(level_root / "img_all" / "mode_summary.csv"),
                        "--mode-export",
                        str(args.mode_export),
                        "--max-rank",
                        str(args.mode_export),
                    ],
                    check=True,
                )

            # 3D direct runs and summaries
            rows_3d_modes: list[dict[str, object]] = []
            for solver in ("fem3d0", "fem3d1"):
                for case3d in ("air", "half", "cyl", "sphere"):
                    if not wants(f"{solver}_{case3d}", "3d"):
                        continue
                    exe = f"{solver}_rect"
                    gx, gy, gz = grids3d[case3d]
                    stdout = run_exe(
                        f"{solver}_{case3d}",
                        "3.1",
                        f"{gx}x{gy}x{gz}",
                        [str(build_dir / exe), f"--{case3d}", "--nx", str(gx), "--ny", str(gy), "--nz", str(gz), "--backend", backend],
                    )
                    for row in parse_mode_table(stdout):
                        rows_3d_modes.append(
                            {
                                "section": "3.1",
                                "solver": solver,
                                "backend": backend,
                                "level": level,
                                "case": case3d,
                                "nx": gx,
                                "ny": gy,
                                "nz": gz,
                                **row,
                            }
                        )

            if rows_3d_modes:
                modes_csv = level_root / "validation" / "validation_3d_31_modes.csv"
                write_csv(
                    modes_csv,
                    rows_3d_modes,
                    ["section", "solver", "backend", "level", "case", "nx", "ny", "nz", "mode_idx", "mode", "k0_ana", "k0_fem", "err_ana_pct", "ref_paper", "err_ref_pct"],
                )
                grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
                for row in rows_3d_modes:
                    grouped.setdefault((str(row["solver"]), str(row["case"])), []).append(row)
                summary_rows: list[dict[str, object]] = []
                for (solver, case3d), rows in sorted(grouped.items()):
                    err_ana = [abs(float(r["err_ana_pct"])) for r in rows]
                    err_ref = [abs(float(r["err_ref_pct"])) for r in rows]
                    sample = rows[0]
                    summary_rows.append(
                        {
                            "section": "3.1",
                            "solver": solver,
                            "backend": backend,
                            "level": level,
                            "case": case3d,
                            "nx": sample["nx"],
                            "ny": sample["ny"],
                            "nz": sample["nz"],
                            "n_modes": len(rows),
                            "max_err_ana_pct": max(err_ana),
                            "mean_err_ana_pct": statistics.fmean(err_ana),
                            "max_err_ref_pct": max(err_ref),
                            "mean_err_ref_pct": statistics.fmean(err_ref),
                        }
                    )
                write_csv(
                    level_root / "validation" / "validation_3d_31_summary.csv",
                    summary_rows,
                    ["section", "solver", "backend", "level", "case", "nx", "ny", "nz", "n_modes", "max_err_ana_pct", "mean_err_ana_pct", "max_err_ref_pct", "mean_err_ref_pct"],
                )

            level_rows.append({
                "backend": backend,
                "level": level,
                **manifest,
                "path": str(level_root.relative_to(out_root)),
            })

    # aggregate csvs
    write_csv(
        out_root / "sweep_runs.csv",
        run_rows,
        ["case_id", "section", "backend", "level", "resolution", "command", "stdout_rel", "stderr_rel", "wall_ms", "assembly_ms", "solve_ms", "post_ms", "total_ms", "returncode"],
    )
    write_csv(
        out_root / "sweep_modes.csv",
        mode_rows,
        ["section", "case", "geometry", "formulation", "polarization", "mode_rank", "m", "index2_name", "index2", "fem", "ref_primary", "err_primary_pct", "rho_abs", "backend", "level", "resolution"],
    )
    write_csv(
        out_root / "sweep_levels.csv",
        level_rows,
        ["backend", "level", "helm10_rect_edge_rect", "polar_2d", "mixed_rect_hv3", "helmvec2_rect", "3d_air", "3d_half", "3d_cyl", "3d_sphere", "path"],
    )

    # Root README
    lines: list[str] = []
    lines.append("# Sweep Completo de Malhas")
    lines.append("")
    lines.append(
        f"Esta pasta contem a varredura completa em `{levels}` niveis de discretizacao "
        f"para os backends `{', '.join(backends)}`."
    )
    lines.append("")
    lines.append("## Como foi gerado")
    lines.append("")
    lines.append("```bash")
    cmd_preview = [
        "python3",
        "scripts/run_full_mesh_sweep.py",
        "--out-root",
        str(out_root.relative_to(ROOT)),
        "--mode-export",
        str(args.mode_export),
        "--levels",
        str(levels),
        "--backends",
        ",".join(backends),
    ]
    if args.skip_images:
        cmd_preview.append("--skip-images")
    lines.append(" ".join(cmd_preview))
    lines.append("```")
    lines.append("")
    lines.append("- CSV de tempos por execucao: [sweep_runs.csv](sweep_runs.csv)")
    lines.append("- CSV de modos/erros 2D: [sweep_modes.csv](sweep_modes.csv)")
    lines.append("- CSV de niveis e resolucoes: [sweep_levels.csv](sweep_levels.csv)")
    lines.append("")
    lines.append("## Niveis de malha")
    lines.append("")
    lines.append("| Backend | Nivel | Rect 2D | Polar 2D | Mixed/HV3 | HV2 | Air 3D | Half 3D | Cyl 3D | Sphere 3D | Pasta |")
    lines.append("|---|---:|---|---|---|---|---|---|---|---|---|")
    for row in level_rows:
        lines.append(
            f"| `{row['backend']}` | {row['level']} | `{row['helm10_rect_edge_rect']}` | `{row['polar_2d']}` | "
            f"`{row['mixed_rect_hv3']}` | `{row['helmvec2_rect']}` | `{row['3d_air']}` | `{row['3d_half']}` | "
            f"`{row['3d_cyl']}` | `{row['3d_sphere']}` | [{row['path']}]({row['path']}) |"
        )
    lines.append("")

    # Timing summary
    lines.append("## Resumo de tempos")
    lines.append("")
    lines.append("| Caso | Backend | Menor total_ms | Maior total_ms | Menor wall_ms | Maior wall_ms |")
    lines.append("|---|---|---:|---:|---:|---:|")
    grouped_runs: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in run_rows:
        grouped_runs.setdefault((str(row["case_id"]), str(row["backend"])), []).append(row)
    for (case_id, backend), rows in sorted(grouped_runs.items()):
        total_vals = [float(r["total_ms"]) for r in rows if str(r["total_ms"])]
        wall_vals = [float(r["wall_ms"]) for r in rows]
        lines.append(
            f"| `{case_id}` | `{backend}` | {fmt(min(total_vals) if total_vals else None)} | "
            f"{fmt(max(total_vals) if total_vals else None)} | {fmt(min(wall_vals) if wall_vals else None)} | {fmt(max(wall_vals) if wall_vals else None)} |"
        )
    lines.append("")

    # Error summary for 2D parsed modes
    lines.append("## Resumo de erros 2D com referencia analitica")
    lines.append("")
    lines.append("| Caso | Backend | Melhor media |err| (%) | Pior media |err| (%) |")
    lines.append("|---|---|---:|---:|")
    grouped_modes: dict[tuple[str, str, int], list[dict[str, object]]] = {}
    for row in mode_rows:
        grouped_modes.setdefault((str(row["case"]), str(row["backend"]), int(row["level"])), []).append(row)
    case_backend_values: dict[tuple[str, str], list[float]] = {}
    for (case, backend, level), rows in grouped_modes.items():
        errs = [abs(float(r["err_primary_pct"])) for r in rows if str(r.get("err_primary_pct", "")) not in {"", "None"}]
        if errs:
            case_backend_values.setdefault((case, backend), []).append(statistics.fmean(errs))
    for (case, backend), vals in sorted(case_backend_values.items()):
        lines.append(f"| `{case}` | `{backend}` | {fmt(min(vals))} | {fmt(max(vals))} |")
    lines.append("")

    lines.append("## Observacoes")
    lines.append("")
    lines.append("- Cada nivel/backend tem sua propria pasta com `2d/`, `validation/`, `img_all/` e `logs/`.")
    lines.append("- Os executaveis 2D geram VTK e, quando imagens nao sao puladas, o script tambem gera PNG e `mode_summary.csv` por nivel.")
    lines.append("- As secoes `2.2.2`, `2.2.3` e `2.2.4` geram `validation_2d_22.csv` e figuras de validacao por nivel.")
    lines.append("- As secoes 3D geram `validation_3d_31_modes.csv` e `validation_3d_31_summary.csv` por nivel.")
    lines.append("- Esta varredura e pesada. Se quiser iterar mais rapido, use `--levels 2` ou `--skip-images`.")
    lines.append("")

    (out_root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Sweep salvo em: {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
