#!/usr/bin/env python3
"""
Generate summary images for Section 2.2 cases using validation_2d_22.csv.

Covers:
- 2.2.2 (mixed cutoff snapshots and rect comparisons)
- 2.2.3 (Figure 11 / Table 8)
- 2.2.4 (Figure 12 / Table 9 and Figure 13 / Table 10)
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _to_float(s: str) -> Optional[float]:
    if s is None:
        return None
    t = s.strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def _read_rows(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")


def _plot_223_table8(rows: List[Dict[str, str]], out_dir: Path, dpi: int) -> None:
    data = [r for r in rows if r["section"] == "2.2.3" and r["case"] == "Figure11_Table8"]
    if not data:
        return
    data.sort(key=lambda r: int(float(r["mode"])))

    x = [int(float(r["mode"])) for r in data]
    fem = [_to_float(r["fem"]) for r in data]
    ref1 = [_to_float(r["ref_primary"]) for r in data]
    ref2 = [_to_float(r["ref_secondary"]) for r in data]
    err1 = [abs(_to_float(r["err_primary_pct"]) or 0.0) for r in data]
    err2 = [abs(_to_float(r["err_secondary_pct"]) or 0.0) for r in data]

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.plot(x, fem, "o-", label="FEM")
    ax.plot(x, ref1, "s--", label="HELMVEC2 (ref)")
    ax.plot(x, ref2, "d-.", label="Hayata (ref)")
    ax.set_xlabel("Mode")
    ax.set_ylabel("k0 L")
    ax.set_title("2.2.3 Figure 11 / Table 8: k0L by mode")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, out_dir / "2_2_3_table8_k0L.png", dpi)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.plot(x, err1, "s--", label="|err| vs HELMVEC2 (%)")
    ax.plot(x, err2, "d-.", label="|err| vs Hayata (%)")
    ax.set_xlabel("Mode")
    ax.set_ylabel("Absolute relative error (%)")
    ax.set_title("2.2.3 Figure 11 / Table 8: error by mode")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, out_dir / "2_2_3_table8_error_pct.png", dpi)


def _plot_224_table9(rows: List[Dict[str, str]], out_dir: Path, dpi: int) -> None:
    data = [r for r in rows if r["section"] == "2.2.4" and r["case"] == "Figure12_Table9"]
    if not data:
        return
    data.sort(key=lambda r: float(r["mode"]))

    x = [float(r["mode"]) for r in data]
    fem = [_to_float(r["fem"]) for r in data]
    ana = [_to_float(r["ref_primary"]) for r in data]
    ref = [_to_float(r["ref_secondary"]) for r in data]
    err_ana = [abs(_to_float(r["err_primary_pct"]) or 0.0) for r in data]
    err_ref = [abs(_to_float(r["err_secondary_pct"]) or 0.0) for r in data]

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.plot(x, fem, "o-", label="FEM")
    ax.plot(x, ana, "s--", label="Analytical (ref)")
    ax.plot(x, ref, "d-.", label="HELMVEC3 (ref)")
    ax.set_xlabel("br / lambda0")
    ax.set_ylabel("beta / k0")
    ax.set_title("2.2.4 Figure 12 / Table 9: dispersion points")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, out_dir / "2_2_4_table9_beta_over_k0.png", dpi)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    ax.plot(x, err_ana, "s--", label="|err| vs Analytical (%)")
    ax.plot(x, err_ref, "d-.", label="|err| vs HELMVEC3 (%)")
    ax.set_xlabel("br / lambda0")
    ax.set_ylabel("Absolute relative error (%)")
    ax.set_title("2.2.4 Figure 12 / Table 9: error")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, out_dir / "2_2_4_table9_error_pct.png", dpi)


def _parse_mode_224_10(mode_s: str) -> Optional[Tuple[float, float]]:
    m = re.match(r"d/a=([0-9.]+),br/lambda0=([0-9.]+)", mode_s.strip())
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def _plot_224_table10(rows: List[Dict[str, str]], out_dir: Path, dpi: int) -> None:
    data = [r for r in rows if r["section"] == "2.2.4" and r["case"] == "Figure13_Table10"]
    if not data:
        return

    groups: Dict[float, List[Tuple[float, float, float, float]]] = {}
    for r in data:
        parsed = _parse_mode_224_10(r["mode"])
        if parsed is None:
            continue
        d_over_a, br_over_l = parsed
        groups.setdefault(d_over_a, []).append(
            (
                br_over_l,
                _to_float(r["fem"]) or float("nan"),
                _to_float(r["ref_primary"]) or float("nan"),
                _to_float(r["ref_secondary"]) or float("nan"),
            )
        )

    if not groups:
        return

    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    for d in sorted(groups.keys()):
        pts = sorted(groups[d], key=lambda t: t[0])
        x = [t[0] for t in pts]
        y = [t[1] for t in pts]
        ax.plot(x, y, marker="o", label=f"FEM d/a={d:g}")
    ax.set_xlabel("br / lambda0")
    ax.set_ylabel("beta / k0")
    ax.set_title("2.2.4 Figure 13 / Table 10: FEM branches by d/a")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    _save(fig, out_dir / "2_2_4_table10_fem_branches.png", dpi)

    # Per d/a comparison with references.
    for d in sorted(groups.keys()):
        pts = sorted(groups[d], key=lambda t: t[0])
        x = [t[0] for t in pts]
        fem = [t[1] for t in pts]
        ana = [t[2] for t in pts]
        ref = [t[3] for t in pts]

        fig, ax = plt.subplots(figsize=(8.8, 4.8))
        ax.plot(x, fem, "o-", label="FEM")
        ax.plot(x, ana, "s--", label="Analytical (ref)")
        ax.plot(x, ref, "d-.", label="HELMVEC3 (ref)")
        ax.set_xlabel("br / lambda0")
        ax.set_ylabel("beta / k0")
        ax.set_title(f"2.2.4 Figure 13 / Table 10: d/a={d:g}")
        ax.grid(True, alpha=0.3)
        ax.legend()
        safe = str(d).replace(".", "p")
        _save(fig, out_dir / f"2_2_4_table10_da_{safe}.png", dpi)


def _plot_222(rows: List[Dict[str, str]], out_dir: Path, dpi: int) -> None:
    # Rectangular mixed case with analytical reference.
    te = [r for r in rows if r["section"] == "2.2.2" and r["case"] == "mixed_rect_E_TE_table"]
    tm = [r for r in rows if r["section"] == "2.2.2" and r["case"] == "mixed_rect_E_TM_table"]
    if te or tm:
        fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.8), sharex=False)
        for ax, data, title in [
            (axs[0], te, "TE (edge block)"),
            (axs[1], tm, "TM (scalar block)"),
        ]:
            data = sorted(data, key=lambda r: int(float(r["mode"])))
            x = [int(float(r["mode"])) for r in data]
            fem = [_to_float(r["fem"]) for r in data]
            ref = [_to_float(r["ref_primary"]) for r in data]
            ax.plot(x, fem, "o-", label="FEM")
            if any(v is not None for v in ref):
                ax.plot(x, ref, "s--", label="Analytical (ref)")
            ax.set_xlabel("Mode rank")
            ax.set_ylabel("kc")
            ax.set_title(title)
            ax.grid(True, alpha=0.3)
            ax.legend()
        fig.suptitle("2.2.2 Mixed rectangular cutoff comparison")
        _save(fig, out_dir / "2_2_2_mixed_rect_cutoff.png", dpi)

    # Circle/coax snapshots (no analytical columns in csv).
    snap_cases = [
        ("mixed_circle_TE_edge", "Circle TE edge"),
        ("mixed_circle_TM_scalar", "Circle TM scalar"),
        ("mixed_coax_TE_edge", "Coax TE edge"),
        ("mixed_coax_TM_scalar", "Coax TM scalar"),
    ]
    fig, axs = plt.subplots(2, 2, figsize=(12.0, 8.0))
    used_any = False
    for ax, (case_name, title) in zip(axs.flatten(), snap_cases):
        data = [r for r in rows if r["section"] == "2.2.2" and r["case"] == case_name]
        data = sorted(data, key=lambda r: int(float(r["mode"])))
        if not data:
            ax.set_visible(False)
            continue
        used_any = True
        x = [int(float(r["mode"])) for r in data]
        fem = [_to_float(r["fem"]) for r in data]
        ax.plot(x, fem, "o-")
        ax.set_xlabel("Mode rank")
        ax.set_ylabel("kc")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
    if used_any:
        fig.suptitle("2.2.2 Mixed circle/coax snapshots")
        _save(fig, out_dir / "2_2_2_mixed_circle_coax_snapshots.png", dpi)
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot images for 2.2.2/2.2.3/2.2.4 from validation_2d_22.csv")
    ap.add_argument("--in-csv", type=Path, default=Path("out/validation/validation_2d_22.csv"))
    ap.add_argument("--out-dir", type=Path, default=Path("out/img_all/validation_2d_22"))
    ap.add_argument("--dpi", type=int, default=210)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    in_csv = _resolve(args.in_csv)
    out_dir = _resolve(args.out_dir)

    if not in_csv.exists():
        raise SystemExit(f"Input CSV not found: {in_csv}")

    rows = _read_rows(in_csv)
    if not rows:
        raise SystemExit(f"Input CSV has no rows: {in_csv}")

    _plot_222(rows, out_dir / "2.2.2", args.dpi)
    _plot_223_table8(rows, out_dir / "2.2.3", args.dpi)
    _plot_224_table9(rows, out_dir / "2.2.4", args.dpi)
    _plot_224_table10(rows, out_dir / "2.2.4", args.dpi)


if __name__ == "__main__":
    main()

