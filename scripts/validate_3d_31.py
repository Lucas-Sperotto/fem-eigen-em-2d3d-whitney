#!/usr/bin/env python3
"""
Run and summarize Section 3.1 (3D cavity) validations into CSV files.

Default outputs:
- out/validation/validation_3d_31_modes.csv
- out/validation/validation_3d_31_summary.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD = ROOT / "build"
DEFAULT_OUT_MODES = ROOT / "out" / "validation" / "validation_3d_31_modes.csv"
DEFAULT_OUT_SUMMARY = ROOT / "out" / "validation" / "validation_3d_31_summary.csv"

TABLE_RE = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$"
)


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def run_checked(build_dir: Path, cmd: list[str], verbose: bool) -> str:
    if verbose:
        print(f"Running: {' '.join(cmd)}")
    return subprocess.check_output(cmd, cwd=build_dir, text=True)


def parse_mode_table(text: str) -> list[dict]:
    rows: list[dict] = []
    in_tbl = False
    for ln in text.splitlines():
        if ln.strip().startswith("idx") and "k0_ana" in ln and "k0_fem" in ln:
            in_tbl = True
            continue
        if not in_tbl:
            continue
        if not ln.strip():
            break

        m = TABLE_RE.match(ln)
        if not m:
            continue

        rows.append(
            {
                "mode_idx": int(m.group(1)),
                "mode": m.group(2),
                "k0_ana": float(m.group(3)),
                "k0_fem": float(m.group(4)),
                "err_ana_pct": float(m.group(5)),
                "ref_paper": float(m.group(6)),
                "err_ref_pct": float(m.group(7)),
            }
        )
    return rows


def grids_for(profile: str) -> dict[str, list[tuple[int, int, int]]]:
    # Quick profile is intended for frequent iteration.
    quick = {
        "air": [(5, 3, 3), (6, 3, 3)],
        "half": [(4, 4, 3), (5, 5, 4)],
        "cyl": [(6, 6, 3), (7, 7, 4)],
        "sphere": [(5, 4, 4), (6, 5, 5)],
    }

    # Full profile adds one denser run per case.
    full = {
        "air": quick["air"] + [(7, 4, 4)],
        "half": quick["half"] + [(6, 6, 5)],
        "cyl": quick["cyl"] + [(8, 8, 5)],
        "sphere": quick["sphere"] + [(7, 6, 6)],
    }
    return full if profile == "full" else quick


def case_option(case: str) -> str:
    return {
        "air": "--air",
        "half": "--half",
        "cyl": "--cyl",
        "sphere": "--sphere",
    }[case]


def solver_bin(solver: str) -> str:
    return {
        "fem3d0": "./fem3d0_rect",
        "fem3d1": "./fem3d1_rect",
    }[solver]


def parse_cases_arg(cases_arg: str) -> list[str]:
    selected_cases = [c.strip() for c in cases_arg.split(",") if c.strip()]
    valid_cases = {"air", "half", "cyl", "sphere"}
    bad = [c for c in selected_cases if c not in valid_cases]
    if bad:
        raise SystemExit(f"Invalid cases: {bad}")
    return selected_cases


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate Section 3.1 and save CSV files.")
    ap.add_argument("--profile", choices=["quick", "full"], default="quick")
    ap.add_argument("--solver", choices=["fem3d0", "fem3d1", "both"], default="both")
    ap.add_argument(
        "--cases",
        default="air,half,cyl,sphere",
        help="Comma-separated subset of: air,half,cyl,sphere",
    )
    ap.add_argument("--build-dir", type=Path, default=Path("build"), help="Build directory containing executables.")
    ap.add_argument("--out-modes", type=Path, default=Path("out/validation/validation_3d_31_modes.csv"), help="Output CSV for per-mode rows.")
    ap.add_argument("--out-summary", type=Path, default=Path("out/validation/validation_3d_31_summary.csv"), help="Output CSV for grouped summary rows.")
    ap.add_argument("--verbose", action="store_true", help="Print executed commands.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    build_dir = _resolve(args.build_dir)
    out_modes = _resolve(args.out_modes)
    out_summary = _resolve(args.out_summary)

    if not build_dir.exists():
        raise SystemExit(f"Build directory not found: {build_dir}")

    selected_cases = parse_cases_arg(args.cases)
    selected_solvers = ["fem3d0", "fem3d1"] if args.solver == "both" else [args.solver]
    case_grids = grids_for(args.profile)

    mode_rows: list[dict] = []

    for solver in selected_solvers:
        exe = solver_bin(solver)
        for case in selected_cases:
            opt = case_option(case)
            for nx, ny, nz in case_grids[case]:
                cmd = [exe, opt, "--nx", str(nx), "--ny", str(ny), "--nz", str(nz)]
                out = run_checked(build_dir, cmd, args.verbose)
                parsed = parse_mode_table(out)
                for row in parsed:
                    mode_rows.append(
                        {
                            "section": "3.1",
                            "solver": solver,
                            "case": case,
                            "nx": nx,
                            "ny": ny,
                            "nz": nz,
                            **row,
                        }
                    )

    out_modes.parent.mkdir(parents=True, exist_ok=True)
    with out_modes.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "section",
                "solver",
                "case",
                "nx",
                "ny",
                "nz",
                "mode_idx",
                "mode",
                "k0_ana",
                "k0_fem",
                "err_ana_pct",
                "ref_paper",
                "err_ref_pct",
            ],
        )
        w.writeheader()
        for row in mode_rows:
            w.writerow(row)

    grouped: dict[tuple[str, str, int, int, int], list[dict]] = defaultdict(list)
    for row in mode_rows:
        key = (row["solver"], row["case"], row["nx"], row["ny"], row["nz"])
        grouped[key].append(row)

    summary_rows: list[dict] = []
    for (solver, case, nx, ny, nz), rows in sorted(grouped.items()):
        err_ana = [abs(float(r["err_ana_pct"])) for r in rows]
        err_ref = [abs(float(r["err_ref_pct"])) for r in rows]
        summary_rows.append(
            {
                "section": "3.1",
                "solver": solver,
                "case": case,
                "nx": nx,
                "ny": ny,
                "nz": nz,
                "n_modes": len(rows),
                "max_err_ana_pct": max(err_ana),
                "mean_err_ana_pct": statistics.mean(err_ana),
                "max_err_ref_pct": max(err_ref),
                "mean_err_ref_pct": statistics.mean(err_ref),
            }
        )

    out_summary.parent.mkdir(parents=True, exist_ok=True)
    with out_summary.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "section",
                "solver",
                "case",
                "nx",
                "ny",
                "nz",
                "n_modes",
                "max_err_ana_pct",
                "mean_err_ana_pct",
                "max_err_ref_pct",
                "mean_err_ref_pct",
            ],
        )
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)

    print(f"Saved: {out_modes}")
    print(f"Saved: {out_summary}")


if __name__ == "__main__":
    main()
