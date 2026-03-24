#!/usr/bin/env python3
"""
Benchmark simples para comparar os backends "gauss" e "closed-form".

O script executa uma bateria de casos representativos, mede:
- tempo de parede externo (wall_ms),
- tempos internos reportados pelos executaveis (assembly/solve/post/total),
e grava:
- CSV detalhado por repeticao;
- CSV agregado por caso/backend;
- resumo em Markdown com speedup.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


TIMING_RE = re.compile(
    r"\[timing\]\s+label=(?P<label>\S+)\s+"
    r"assembly_ms=(?P<assembly>[0-9.]+)\s+"
    r"solve_ms=(?P<solve>[0-9.]+)\s+"
    r"post_ms=(?P<post>[0-9.]+)\s+"
    r"total_ms=(?P<total>[0-9.]+)"
)


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    section: str
    exe: str
    args: tuple[str, ...]
    note: str


CORE_CASES = [
    CaseSpec("helm10_rect", "2.1", "helm10_rect", ("14", "14", "0"), "retangular escalar, sem exportacao"),
    CaseSpec("edge_rect", "2.2.1", "edge_rect", ("14", "14", "0"), "retangular vetorial, sem exportacao"),
    CaseSpec("mixed_rect", "2.2.2", "mixed_rect", ("12", "6"), "sistema misto retangular"),
    CaseSpec("helmvec2_rect", "2.2.3", "helmvec2_rect", ("10", "6", "6", "0"), "k0 dado beta, Figura 11"),
    CaseSpec("helmvec3_rect", "2.2.4", "helmvec3_rect", ("0.20", "10", "5", "0"), "beta dado k0, Fig. 12/13"),
    CaseSpec("fem3d0_air", "3.1", "fem3d0_rect", ("--air", "--nx", "8", "--ny", "4", "--nz", "6"), "cavidade retangular 3D densa"),
    CaseSpec("fem3d1_air", "3.1", "fem3d1_rect", ("--air", "--nx", "8", "--ny", "4", "--nz", "6"), "cavidade retangular 3D esparsa"),
]

ALL_CASES = CORE_CASES + [
    CaseSpec("helm10_circle", "2.1", "helm10_circle", ("16", "72", "0"), "circular escalar, sem exportacao"),
    CaseSpec("helm10_coax", "2.1", "helm10_coax", ("16", "72", "0"), "coaxial escalar, sem exportacao"),
    CaseSpec("edge_circle", "2.2.1", "edge_circle", ("14", "72", "0"), "circular vetorial, sem exportacao"),
    CaseSpec("edge_coax", "2.2.1", "edge_coax", ("16", "96", "0"), "coaxial vetorial, sem exportacao"),
    CaseSpec("mixed_circle", "2.2.2", "mixed_circle", ("14", "72"), "sistema misto circular"),
    CaseSpec("mixed_coax", "2.2.2", "mixed_coax", ("14", "72"), "sistema misto coaxial"),
    CaseSpec("fem3d0_half", "3.1", "fem3d0_rect", ("--half", "--nx", "8", "--ny", "6", "--nz", "6"), "cavidade meia-preenchida 3D densa"),
    CaseSpec("fem3d1_half", "3.1", "fem3d1_rect", ("--half", "--nx", "8", "--ny", "6", "--nz", "6"), "cavidade meia-preenchida 3D esparsa"),
]


def pick_cases(suite: str) -> list[CaseSpec]:
    if suite == "core":
        return CORE_CASES
    if suite == "all":
        return ALL_CASES
    raise ValueError(f"suite invalida: {suite}")


def mean_or_nan(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def parse_timing(stdout: str) -> dict[str, float]:
    match = TIMING_RE.search(stdout)
    if not match:
        return {
            "assembly_ms": math.nan,
            "solve_ms": math.nan,
            "post_ms": math.nan,
            "total_ms": math.nan,
        }
    return {
        "assembly_ms": float(match.group("assembly")),
        "solve_ms": float(match.group("solve")),
        "post_ms": float(match.group("post")),
        "total_ms": float(match.group("total")),
    }


def run_case(build_dir: Path, case: CaseSpec, backend: str) -> dict[str, object]:
    exe = build_dir / case.exe
    if not exe.exists():
        raise FileNotFoundError(f"Executavel nao encontrado: {exe}")

    cmd = [str(exe), *case.args, "--backend", backend]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall_ms = (time.perf_counter() - t0) * 1000.0
    timing = parse_timing(proc.stdout)
    return {
        "command": " ".join(cmd),
        "returncode": proc.returncode,
        "wall_ms": wall_ms,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        **timing,
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description="Benchmark gauss vs closed-form.")
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--out-dir", default="out/benchmark")
    ap.add_argument("--suite", choices=("core", "all"), default="core")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--keep-logs", action="store_true", help="Salva stdout/stderr bruto de cada execucao.")
    args = ap.parse_args()

    root_dir = Path(__file__).resolve().parents[1]
    build_dir = (root_dir / args.build_dir).resolve()
    out_dir = (root_dir / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = out_dir / "raw"
    if args.keep_logs:
      raw_dir.mkdir(parents=True, exist_ok=True)

    cases = pick_cases(args.suite)
    detail_rows: list[dict[str, object]] = []

    for case in cases:
        for backend in ("gauss", "closed-form"):
            for repeat in range(1, args.repeats + 1):
                result = run_case(build_dir, case, backend)
                row = {
                    "section": case.section,
                    "case_id": case.case_id,
                    "backend": backend,
                    "repeat": repeat,
                    "wall_ms": f"{result['wall_ms']:.3f}",
                    "assembly_ms": f"{result['assembly_ms']:.3f}" if not math.isnan(result["assembly_ms"]) else "",
                    "solve_ms": f"{result['solve_ms']:.3f}" if not math.isnan(result["solve_ms"]) else "",
                    "post_ms": f"{result['post_ms']:.3f}" if not math.isnan(result["post_ms"]) else "",
                    "total_ms": f"{result['total_ms']:.3f}" if not math.isnan(result["total_ms"]) else "",
                    "returncode": result["returncode"],
                    "command": result["command"],
                    "note": case.note,
                }
                detail_rows.append(row)

                print(
                    f"[bench] {case.case_id:<14} backend={backend:<11} rep={repeat} "
                    f"wall_ms={result['wall_ms']:.3f} rc={result['returncode']}",
                    flush=True,
                )

                if args.keep_logs:
                    stem = f"{case.case_id}__{backend.replace('-', '_')}__rep{repeat}"
                    (raw_dir / f"{stem}.stdout.txt").write_text(str(result["stdout"]), encoding="utf-8")
                    (raw_dir / f"{stem}.stderr.txt").write_text(str(result["stderr"]), encoding="utf-8")

                if result["returncode"] != 0:
                    raise RuntimeError(f"Falha ao executar: {result['command']}")

    detail_csv = out_dir / "backend_benchmark_detail.csv"
    write_csv(
        detail_csv,
        detail_rows,
        [
            "section",
            "case_id",
            "backend",
            "repeat",
            "wall_ms",
            "assembly_ms",
            "solve_ms",
            "post_ms",
            "total_ms",
            "returncode",
            "command",
            "note",
        ],
    )

    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in detail_rows:
        grouped.setdefault((str(row["case_id"]), str(row["backend"])), []).append(row)

    aggregate_rows: list[dict[str, object]] = []
    for case in cases:
        for backend in ("gauss", "closed-form"):
            rows = grouped[(case.case_id, backend)]
            wall_vals = [float(r["wall_ms"]) for r in rows]
            assembly_vals = [float(r["assembly_ms"]) for r in rows if str(r["assembly_ms"])]
            solve_vals = [float(r["solve_ms"]) for r in rows if str(r["solve_ms"])]
            post_vals = [float(r["post_ms"]) for r in rows if str(r["post_ms"])]
            total_vals = [float(r["total_ms"]) for r in rows if str(r["total_ms"])]
            aggregate_rows.append(
                {
                    "section": case.section,
                    "case_id": case.case_id,
                    "backend": backend,
                    "repeats": len(rows),
                    "wall_ms_mean": f"{mean_or_nan(wall_vals):.3f}",
                    "wall_ms_min": f"{min(wall_vals):.3f}",
                    "assembly_ms_mean": f"{mean_or_nan(assembly_vals):.3f}" if assembly_vals else "",
                    "solve_ms_mean": f"{mean_or_nan(solve_vals):.3f}" if solve_vals else "",
                    "post_ms_mean": f"{mean_or_nan(post_vals):.3f}" if post_vals else "",
                    "total_ms_mean": f"{mean_or_nan(total_vals):.3f}" if total_vals else "",
                    "note": case.note,
                }
            )

    aggregate_csv = out_dir / "backend_benchmark_summary.csv"
    write_csv(
        aggregate_csv,
        aggregate_rows,
        [
            "section",
            "case_id",
            "backend",
            "repeats",
            "wall_ms_mean",
            "wall_ms_min",
            "assembly_ms_mean",
            "solve_ms_mean",
            "post_ms_mean",
            "total_ms_mean",
            "note",
        ],
    )

    md_lines = [
        "# Benchmark de Backends",
        "",
        f"- Suite: `{args.suite}`",
        f"- Repeticoes por caso/backend: `{args.repeats}`",
        f"- CSV detalhado: [{detail_csv.name}]({detail_csv.name})",
        f"- CSV agregado: [{aggregate_csv.name}]({aggregate_csv.name})",
        "",
        "| Secao | Caso | Wall Gauss (ms) | Wall Closed (ms) | Speedup Closed/Gauss | Assembly Gauss (ms) | Assembly Closed (ms) |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]

    for case in cases:
        gauss = next(r for r in aggregate_rows if r["case_id"] == case.case_id and r["backend"] == "gauss")
        closed = next(r for r in aggregate_rows if r["case_id"] == case.case_id and r["backend"] == "closed-form")
        wall_gauss = float(gauss["wall_ms_mean"])
        wall_closed = float(closed["wall_ms_mean"])
        speedup = wall_gauss / wall_closed if wall_closed > 0.0 else math.nan

        def fmt_maybe(v: str) -> str:
            return v if v else "-"

        md_lines.append(
            f"| {case.section} | `{case.case_id}` | {wall_gauss:.3f} | {wall_closed:.3f} | {speedup:.3f}x | "
            f"{fmt_maybe(str(gauss['assembly_ms_mean']))} | {fmt_maybe(str(closed['assembly_ms_mean']))} |"
        )

    report_md = out_dir / "BACKEND_BENCHMARK.md"
    report_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"\n[bench] CSV detalhado : {detail_csv}")
    print(f"[bench] CSV agregado : {aggregate_csv}")
    print(f"[bench] Relatorio MD : {report_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
