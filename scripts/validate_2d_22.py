#!/usr/bin/env python3
"""
Run and summarize Section 2.2 validations into one CSV file.

Default output:
- out/validation/validation_2d_22.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILD = ROOT / "build"
DEFAULT_OUT = ROOT / "out" / "validation" / "validation_2d_22.csv"


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def run_checked(build_dir: Path, cmd: list[str], verbose: bool) -> str:
    if verbose:
        print(f"Running: {' '.join(cmd)}")
    return subprocess.check_output(cmd, cwd=build_dir, text=True)


def parse_helmvec2_table(text: str) -> list[dict]:
    rows: list[dict] = []
    pat = re.compile(r"^\s*(\d+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$")
    in_tbl = False
    for ln in text.splitlines():
        if "k0L(FEM matched)" in ln:
            in_tbl = True
            continue
        if not in_tbl:
            continue
        m = pat.match(ln)
        if not m:
            continue

        mode = int(m.group(1))
        fem = float(m.group(2))
        ref_helmvec2 = float(m.group(3))
        ref_hayata = float(m.group(4))
        err_hvec = 100.0 * (fem - ref_helmvec2) / ref_helmvec2
        err_hay = 100.0 * (fem - ref_hayata) / ref_hayata
        rows.append(
            {
                "section": "2.2.3",
                "case": "Figure11_Table8",
                "mode": mode,
                "fem": fem,
                "ref_primary": ref_helmvec2,
                "ref_secondary": ref_hayata,
                "err_primary_pct": err_hvec,
                "err_secondary_pct": err_hay,
            }
        )
    return rows


def parse_helmvec3_table9(text: str) -> list[dict]:
    rows: list[dict] = []
    pat = re.compile(r"^\s*([0-9.]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$")
    in_tbl = False
    for ln in text.splitlines():
        if "br/lambda0  beta/k0(FEM)" in ln and "Analytic" in ln:
            in_tbl = True
            continue
        if in_tbl and ln.strip().startswith("[2.2.4] beta from given k0 (Figure 13"):
            break
        if not in_tbl:
            continue

        m = pat.match(ln)
        if not m:
            continue

        br_over_lambda = float(m.group(1))
        fem = float(m.group(2))
        ref_ana = float(m.group(3))
        ref_hvec3 = float(m.group(4))
        err_ana = 100.0 * (fem - ref_ana) / ref_ana
        err_hv3 = 100.0 * (fem - ref_hvec3) / ref_hvec3
        rows.append(
            {
                "section": "2.2.4",
                "case": "Figure12_Table9",
                "mode": br_over_lambda,
                "fem": fem,
                "ref_primary": ref_ana,
                "ref_secondary": ref_hvec3,
                "err_primary_pct": err_ana,
                "err_secondary_pct": err_hv3,
            }
        )
    return rows


def parse_helmvec3_table10(text: str) -> list[dict]:
    rows: list[dict] = []
    pat = re.compile(
        r"^\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$"
    )
    in_tbl = False
    for ln in text.splitlines():
        if "d/a  br/lambda0  beta/k0(FEM matched)  Analytical(ref)  HELMVEC3(ref)" in ln:
            in_tbl = True
            continue
        if not in_tbl:
            continue

        m = pat.match(ln)
        if not m:
            continue

        d_over_a = float(m.group(1))
        br_over_lambda = float(m.group(2))
        fem = float(m.group(3))
        ref_ana = float(m.group(4))
        ref_hvec3 = float(m.group(5))
        err_ana = 100.0 * (fem - ref_ana) / ref_ana
        err_hv3 = 100.0 * (fem - ref_hvec3) / ref_hvec3
        rows.append(
            {
                "section": "2.2.4",
                "case": "Figure13_Table10",
                "mode": f"d/a={d_over_a},br/lambda0={br_over_lambda}",
                "fem": fem,
                "ref_primary": ref_ana,
                "ref_secondary": ref_hvec3,
                "err_primary_pct": err_ana,
                "err_secondary_pct": err_hv3,
            }
        )
    return rows


def parse_first_kc_block(text: str, label: str) -> list[float]:
    vals: list[float] = []
    start = text.find(label)
    if start < 0:
        return vals

    sub = text[start:].splitlines()
    for ln in sub[1:]:
        if not ln.strip():
            break
        m = re.search(r"^\s*\d+\s+([0-9.+-eE]+)\s*$", ln)
        if m:
            vals.append(float(m.group(1)))
        else:
            break
    return vals


def parse_mixed_rect_table(text: str, block_title: str, case: str) -> list[dict]:
    rows: list[dict] = []
    start = text.find(block_title)
    if start < 0:
        return rows

    sub = text[start:].splitlines()
    pat = re.compile(r"^\s*(\d+)\s+\((\d+),(\d+)\)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s+([0-9.+-eE]+)\s*$")
    for ln in sub:
        m = pat.match(ln)
        if not m:
            if rows and not ln.strip():
                break
            continue

        mode = int(m.group(1))
        kc_ana = float(m.group(4))
        kc_fem = float(m.group(5))
        err = float(m.group(6))
        rows.append(
            {
                "section": "2.2.2",
                "case": case,
                "mode": mode,
                "fem": kc_fem,
                "ref_primary": kc_ana,
                "ref_secondary": "",
                "err_primary_pct": err,
                "err_secondary_pct": "",
            }
        )
    return rows


def append_snapshot_rows(rows: list[dict], snapshots: Iterable[tuple[str, str, list[float]]]) -> None:
    for sec, case, vals in snapshots:
        for i, val in enumerate(vals, start=1):
            rows.append(
                {
                    "section": sec,
                    "case": case,
                    "mode": i,
                    "fem": val,
                    "ref_primary": "",
                    "ref_secondary": "",
                    "err_primary_pct": "",
                    "err_secondary_pct": "",
                }
            )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate Section 2.2 and save CSV.")
    ap.add_argument("--build-dir", type=Path, default=Path("build"), help="Build directory containing executables.")
    ap.add_argument("--out-csv", type=Path, default=Path("out/validation/validation_2d_22.csv"), help="Output CSV path.")
    ap.add_argument("--rect-nx", type=int, default=12)
    ap.add_argument("--rect-ny", type=int, default=6)
    ap.add_argument("--circle-nr", type=int, default=10)
    ap.add_argument("--circle-nt", type=int, default=48)
    ap.add_argument("--coax-nr", type=int, default=10)
    ap.add_argument("--coax-nt", type=int, default=48)
    ap.add_argument("--beta", type=float, default=10.0, help="beta*L for helmvec2_rect.")
    ap.add_argument("--hv2-nx", type=int, default=6)
    ap.add_argument("--hv2-ny", type=int, default=6)
    ap.add_argument("--d-over-a", type=float, default=0.20, help="Figure 13 preview parameter for helmvec3_rect.")
    ap.add_argument("--hv3-nx", type=int, default=10)
    ap.add_argument("--hv3-ny", type=int, default=5)
    ap.add_argument("--verbose", action="store_true", help="Print executed commands.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    build_dir = _resolve(args.build_dir)
    out_csv = _resolve(args.out_csv)

    if not build_dir.exists():
        raise SystemExit(f"Build directory not found: {build_dir}")

    out_mixed_rect = run_checked(build_dir, ["./mixed_rect", str(args.rect_nx), str(args.rect_ny)], args.verbose)
    out_mixed_circle = run_checked(build_dir, ["./mixed_circle", str(args.circle_nr), str(args.circle_nt)], args.verbose)
    out_mixed_coax = run_checked(build_dir, ["./mixed_coax", str(args.coax_nr), str(args.coax_nt)], args.verbose)
    out_hv2 = run_checked(
        build_dir,
        ["./helmvec2_rect", f"{args.beta}", str(args.hv2_nx), str(args.hv2_ny)],
        args.verbose,
    )
    out_hv3 = run_checked(
        build_dir,
        ["./helmvec3_rect", f"{args.d_over_a}", str(args.hv3_nx), str(args.hv3_ny)],
        args.verbose,
    )

    rows: list[dict] = []
    rows.extend(parse_helmvec2_table(out_hv2))
    rows.extend(parse_helmvec3_table9(out_hv3))
    rows.extend(parse_helmvec3_table10(out_hv3))
    rows.extend(
        parse_mixed_rect_table(
            out_mixed_rect,
            "[E-formulation] TE cutoffs (edge block)",
            "mixed_rect_E_TE_table",
        )
    )
    rows.extend(
        parse_mixed_rect_table(
            out_mixed_rect,
            "[E-formulation] TM cutoffs (scalar block)",
            "mixed_rect_E_TM_table",
        )
    )

    snapshots = [
        ("2.2.2", "mixed_circle_TE_edge", parse_first_kc_block(out_mixed_circle, "TE (edge block), first 8 kc:")),
        ("2.2.2", "mixed_circle_TM_scalar", parse_first_kc_block(out_mixed_circle, "TM (scalar block), first 8 kc:")),
        ("2.2.2", "mixed_coax_TE_edge", parse_first_kc_block(out_mixed_coax, "TE (edge block), first 8 kc:")),
        ("2.2.2", "mixed_coax_TM_scalar", parse_first_kc_block(out_mixed_coax, "TM (scalar block), first 8 kc:")),
    ]
    append_snapshot_rows(rows, snapshots)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "section",
                "case",
                "mode",
                "fem",
                "ref_primary",
                "ref_secondary",
                "err_primary_pct",
                "err_secondary_pct",
            ],
        )
        w.writeheader()
        for row in rows:
            w.writerow(row)

    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()
