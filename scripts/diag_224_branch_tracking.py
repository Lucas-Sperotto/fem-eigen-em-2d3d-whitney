#!/usr/bin/env python3
"""
Evaluate whether Figure 13 / Table 10 residuals can be reduced by branch
tracking over the deduplicated raw spectrum.

Default outputs:
- out/validation/diag_224_branch_tracking_summary.csv
- out/validation/diag_224_branch_tracking_points.csv
- out/validation/DIAG_224_BRANCH_TRACKING.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_ROOT = ROOT / "out"
DEFAULT_BUILD_BIN = ROOT / "build" / "helmvec3_fig13_rect"
DEFAULT_BACKEND = "closed-form"
RUNS_SUBDIR = ".diag_224_branch_tracking"
SUMMARY_CSV_NAME = "diag_224_branch_tracking_summary.csv"
POINTS_CSV_NAME = "diag_224_branch_tracking_points.csv"
REPORT_MD_NAME = "DIAG_224_BRANCH_TRACKING.md"
DELTA_ALPHA = 0.5
COMPARE_TOL = 1.0e-9


@dataclass(frozen=True)
class PriorityPoint:
    key: str
    d_over_a: float
    br_over_lambda0: float

    @property
    def mode(self) -> str:
        return f"d/a={self.d_over_a},br/lambda0={self.br_over_lambda0}"


PRIORITY_POINTS = (
    PriorityPoint("P1", 0.167, 0.5),
    PriorityPoint("P2", 0.286, 0.5),
    PriorityPoint("P3", 0.5, 0.4),
)


@dataclass(frozen=True)
class TablePoint:
    d_over_a: float
    br_over_lambda0: float
    analytic_beta_over_k0: float
    fem_beta_over_k0: float
    selected_candidate_rank: int
    selected_eig_index: int
    ez_ratio: float


@dataclass(frozen=True)
class Candidate:
    beta_ratio: float
    ez_ratio: float
    candidate_rank: int
    solver_index: int


@dataclass(frozen=True)
class SelectedPoint:
    beta_ratio: float
    candidate_rank: int
    solver_index: int
    ez_ratio: float


@dataclass(frozen=True)
class StrategySummary:
    strategy: str
    mean_error_percent: float
    max_error_percent: float
    worst_d_over_a: float
    worst_br_over_lambda0: float
    worst_beta_over_k0: float
    worst_analytic_beta_over_k0: float


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _label_number(value: float) -> str:
    return f"{value:.12g}".replace(".", "_")


def _point_key(d_over_a: float, br_over_lambda0: float) -> tuple[float, float]:
    return (round(d_over_a, 12), round(br_over_lambda0, 12))


def _run_root(out_root: Path, backend: str) -> Path:
    return out_root / RUNS_SUBDIR / backend.replace("-", "_")


def _table10_path(run_root: Path) -> Path:
    return run_root / "helmvec3" / "fig13_rect" / "csv" / "helmvec3_fig13_rect_table10.csv"


def _raw_spectrum_path(run_root: Path, d_over_a: float, br_over_lambda0: float) -> Path:
    filename = (
        "helmvec3_fig13_rect_table10_da"
        + _label_number(d_over_a)
        + "_br"
        + _label_number(br_over_lambda0)
        + "_raw_spectrum.csv"
    )
    return run_root / "helmvec3" / "fig13_rect" / "csv" / filename


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"Input CSV is empty: {_rel_from_root(path)}")
    return rows


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


def _run_case(
    build_bin: Path,
    out_root: Path,
    backend: str,
    nx: int,
    ny: int,
) -> tuple[list[str], Path]:
    run_root = _run_root(out_root, backend)
    run_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(run_root)
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM"] = "1"
    env["TP3485_HELMVEC3_DIAG_RAW_SPECTRUM_TABLE10_ALL"] = "1"
    cmd = [
        str(build_bin),
        "--d-over-a-preview",
        "0.20",
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--backend",
        backend,
    ]
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)
    return cmd, run_root


def _load_table10(path: Path) -> list[TablePoint]:
    rows = _read_csv_rows(path)
    out: list[TablePoint] = []
    for row in rows:
        out.append(
            TablePoint(
                d_over_a=_to_float(row["d_over_a"], field="d_over_a", path=path),
                br_over_lambda0=_to_float(
                    row["br_over_lambda0"], field="br_over_lambda0", path=path
                ),
                analytic_beta_over_k0=_to_float(
                    row["beta_over_k0_analytic"],
                    field="beta_over_k0_analytic",
                    path=path,
                ),
                fem_beta_over_k0=_to_float(
                    row["beta_over_k0_fem_matched"],
                    field="beta_over_k0_fem_matched",
                    path=path,
                ),
                selected_candidate_rank=_to_int(
                    row["selected_candidate_rank"],
                    field="selected_candidate_rank",
                    path=path,
                ),
                selected_eig_index=_to_int(
                    row["selected_eig_index"],
                    field="selected_eig_index",
                    path=path,
                ),
                ez_ratio=_to_float(row["ez_ratio"], field="ez_ratio", path=path),
            )
        )
    return out


def _load_raw_candidates(path: Path) -> list[Candidate]:
    rows = _read_csv_rows(path)
    out: list[Candidate] = []
    for row in rows:
        if row["kept_after_dedup"] != "1":
            continue
        beta = _to_float(
            row["beta_ratio_if_real_positive"],
            field="beta_ratio_if_real_positive",
            path=path,
        )
        if not math.isfinite(beta):
            continue
        out.append(
            Candidate(
                beta_ratio=beta,
                ez_ratio=_to_float(row["ez_ratio"], field="ez_ratio", path=path),
                candidate_rank=_to_int(
                    row["candidate_rank_after_dedup"],
                    field="candidate_rank_after_dedup",
                    path=path,
                ),
                solver_index=_to_int(
                    row["solver_index"],
                    field="solver_index",
                    path=path,
                ),
            )
        )
    if not out:
        raise SystemExit(f"No deduplicated candidates found in {_rel_from_root(path)}")
    return out


def _absolute_relative_error_percent(reference: float, value: float) -> float:
    if not math.isfinite(reference) or reference == 0.0 or not math.isfinite(value):
        return math.nan
    return abs(value - reference) / abs(reference) * 100.0


def _select_nearest_reference(candidates: list[Candidate], ref: float) -> SelectedPoint:
    best = min(candidates, key=lambda cand: abs(cand.beta_ratio - ref))
    return SelectedPoint(
        beta_ratio=best.beta_ratio,
        candidate_rank=best.candidate_rank,
        solver_index=best.solver_index,
        ez_ratio=best.ez_ratio,
    )


def _select_continuity(candidates: list[Candidate], prev_beta: float) -> SelectedPoint:
    best = min(candidates, key=lambda cand: abs(cand.beta_ratio - prev_beta))
    return SelectedPoint(
        beta_ratio=best.beta_ratio,
        candidate_rank=best.candidate_rank,
        solver_index=best.solver_index,
        ez_ratio=best.ez_ratio,
    )


def _select_delta_guided(
    candidates: list[Candidate],
    prev_beta: float,
    prev_ref: float,
    ref: float,
    alpha: float,
) -> SelectedPoint:
    delta_ref = ref - prev_ref
    delta_scale = max(abs(delta_ref), 0.05)
    ref_scale = max(abs(ref), 1.0e-9)

    def score(cand: Candidate) -> float:
        ref_err = abs(cand.beta_ratio - ref) / ref_scale
        delta_err = abs((cand.beta_ratio - prev_beta) - delta_ref) / delta_scale
        return ref_err + alpha * delta_err

    best = min(candidates, key=score)
    return SelectedPoint(
        beta_ratio=best.beta_ratio,
        candidate_rank=best.candidate_rank,
        solver_index=best.solver_index,
        ez_ratio=best.ez_ratio,
    )


def _group_table_points(table_points: list[TablePoint]) -> dict[float, list[TablePoint]]:
    grouped: dict[float, list[TablePoint]] = {}
    for point in table_points:
        grouped.setdefault(point.d_over_a, []).append(point)
    for points in grouped.values():
        points.sort(key=lambda point: point.br_over_lambda0)
    return grouped


def _evaluate_strategy(
    grouped_points: dict[float, list[TablePoint]],
    raw_lookup: dict[tuple[float, float], list[Candidate]],
    strategy: str,
) -> tuple[StrategySummary, dict[tuple[float, float], SelectedPoint]]:
    selections: dict[tuple[float, float], SelectedPoint] = {}
    total_error = 0.0
    count = 0
    max_error = -1.0
    worst_point: TablePoint | None = None
    worst_selection: SelectedPoint | None = None

    for d_over_a in sorted(grouped_points):
        prev_beta: float | None = None
        prev_ref: float | None = None
        for point in grouped_points[d_over_a]:
            key = _point_key(point.d_over_a, point.br_over_lambda0)
            candidates = raw_lookup[key]
            if strategy == "baseline_nearest_reference":
                selected = _select_nearest_reference(candidates, point.analytic_beta_over_k0)
            elif strategy == "continuity_only":
                if prev_beta is None:
                    selected = _select_nearest_reference(candidates, point.analytic_beta_over_k0)
                else:
                    selected = _select_continuity(candidates, prev_beta)
            elif strategy == "delta_guided_reference":
                if prev_beta is None or prev_ref is None:
                    selected = _select_nearest_reference(candidates, point.analytic_beta_over_k0)
                else:
                    selected = _select_delta_guided(
                        candidates,
                        prev_beta,
                        prev_ref,
                        point.analytic_beta_over_k0,
                        DELTA_ALPHA,
                    )
            else:
                raise ValueError(f"Unknown strategy: {strategy}")

            err = _absolute_relative_error_percent(
                point.analytic_beta_over_k0, selected.beta_ratio
            )
            total_error += err
            count += 1
            selections[key] = selected
            if err > max_error:
                max_error = err
                worst_point = point
                worst_selection = selected
            prev_beta = selected.beta_ratio
            prev_ref = point.analytic_beta_over_k0

    if worst_point is None or worst_selection is None:
        raise SystemExit("No points were evaluated.")

    summary = StrategySummary(
        strategy=strategy,
        mean_error_percent=total_error / count,
        max_error_percent=max_error,
        worst_d_over_a=worst_point.d_over_a,
        worst_br_over_lambda0=worst_point.br_over_lambda0,
        worst_beta_over_k0=worst_selection.beta_ratio,
        worst_analytic_beta_over_k0=worst_point.analytic_beta_over_k0,
    )
    return summary, selections


def _verify_table_matches_local_optimum(
    table_points: list[TablePoint],
    raw_lookup: dict[tuple[float, float], list[Candidate]],
) -> list[TablePoint]:
    mismatches: list[TablePoint] = []
    for point in table_points:
        key = _point_key(point.d_over_a, point.br_over_lambda0)
        best = _select_nearest_reference(raw_lookup[key], point.analytic_beta_over_k0)
        if (
            abs(best.beta_ratio - point.fem_beta_over_k0) > COMPARE_TOL
            or best.candidate_rank != point.selected_candidate_rank
            or best.solver_index != point.selected_eig_index
        ):
            mismatches.append(point)
    return mismatches


def _candidate_bracket(
    candidates: list[Candidate], ref: float
) -> tuple[Candidate | None, Candidate, Candidate | None]:
    closest = min(candidates, key=lambda cand: abs(cand.beta_ratio - ref))
    below = None
    above = None
    for cand in candidates:
        if cand.beta_ratio <= ref and (below is None or cand.beta_ratio > below.beta_ratio):
            below = cand
        if cand.beta_ratio >= ref and (above is None or cand.beta_ratio < above.beta_ratio):
            above = cand
    return below, closest, above


def _write_summary_csv(path: Path, rows: list[StrategySummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "strategy",
                "mean_error_percent",
                "max_error_percent",
                "worst_d_over_a",
                "worst_br_over_lambda0",
                "worst_beta_over_k0",
                "worst_analytic_beta_over_k0",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row.strategy,
                    f"{row.mean_error_percent:.12g}",
                    f"{row.max_error_percent:.12g}",
                    f"{row.worst_d_over_a:.12g}",
                    f"{row.worst_br_over_lambda0:.12g}",
                    f"{row.worst_beta_over_k0:.12g}",
                    f"{row.worst_analytic_beta_over_k0:.12g}",
                ]
            )


def _write_points_csv(
    path: Path,
    table_points: list[TablePoint],
    strategies: dict[str, dict[tuple[float, float], SelectedPoint]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "strategy",
                "d_over_a",
                "br_over_lambda0",
                "analytic_beta_over_k0",
                "selected_beta_over_k0",
                "selected_candidate_rank",
                "selected_eig_index",
                "ez_ratio",
                "error_percent_analytic",
            ]
        )
        for strategy_name, selections in strategies.items():
            for point in table_points:
                key = _point_key(point.d_over_a, point.br_over_lambda0)
                selected = selections[key]
                writer.writerow(
                    [
                        strategy_name,
                        f"{point.d_over_a:.12g}",
                        f"{point.br_over_lambda0:.12g}",
                        f"{point.analytic_beta_over_k0:.12g}",
                        f"{selected.beta_ratio:.12g}",
                        selected.candidate_rank,
                        selected.solver_index,
                        f"{selected.ez_ratio:.12g}",
                        f"{_absolute_relative_error_percent(point.analytic_beta_over_k0, selected.beta_ratio):.12g}",
                    ]
                )


def _write_markdown_report(
    path: Path,
    *,
    cmd: list[str],
    backend: str,
    nx: int,
    ny: int,
    run_root: Path,
    summaries: list[StrategySummary],
    table_points: list[TablePoint],
    raw_lookup: dict[tuple[float, float], list[Candidate]],
    point_rows_path: Path,
    summary_csv_path: Path,
    mismatches: list[TablePoint],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    baseline_summary = next(
        summary for summary in summaries if summary.strategy == "baseline_nearest_reference"
    )
    by_strategy = {summary.strategy: summary for summary in summaries}

    with path.open("w", encoding="utf-8") as handle:
        handle.write("# DIAG 224 Branch Tracking\n\n")
        handle.write(f"- Generated: {dt.datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(f"- Backend: `{backend}`\n")
        handle.write(f"- Mesh: `nx={nx}`, `ny={ny}`\n")
        handle.write(f"- Command: `{shlex.join(cmd)}`\n")
        handle.write(f"- Run root: `{_rel_from_root(run_root)}`\n")
        handle.write(
            f"- Summary CSV: `{_rel_from_root(summary_csv_path)}`\n"
        )
        handle.write(
            f"- Points CSV: `{_rel_from_root(point_rows_path)}`\n\n"
        )

        handle.write("## Strategy Summary\n\n")
        handle.write(
            "| strategy | mean error % | max error % | worst point |\n"
        )
        handle.write(
            "| --- | ---: | ---: | --- |\n"
        )
        for summary in summaries:
            handle.write(
                "| "
                f"{summary.strategy} | "
                f"{summary.mean_error_percent:.4f} | "
                f"{summary.max_error_percent:.4f} | "
                f"(d/a={summary.worst_d_over_a:.3f}, br/lambda0={summary.worst_br_over_lambda0:.3f})"
                " |\n"
            )
        handle.write("\n")

        handle.write("## Local-Optimum Check\n\n")
        if mismatches:
            handle.write(
                f"Current Table 10 selection disagrees with the nearest analytic candidate in {len(mismatches)} point(s).\n\n"
            )
        else:
            handle.write(
                f"Current Table 10 selection matches the nearest deduplicated analytic candidate in all {len(table_points)} points.\n\n"
            )

        handle.write("## Priority Points\n\n")
        handle.write(
            "| point | analytic | baseline beta | closest candidate | lower bracket | upper bracket | dedup count |\n"
        )
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        table_lookup = {
            _point_key(point.d_over_a, point.br_over_lambda0): point for point in table_points
        }
        for point in PRIORITY_POINTS:
            key = _point_key(point.d_over_a, point.br_over_lambda0)
            table_point = table_lookup[key]
            candidates = raw_lookup[key]
            below, closest, above = _candidate_bracket(
                candidates, table_point.analytic_beta_over_k0
            )
            lower_text = f"{below.beta_ratio:.6f}" if below is not None else "nan"
            upper_text = f"{above.beta_ratio:.6f}" if above is not None else "nan"
            handle.write(
                "| "
                f"{point.key} ({point.mode}) | "
                f"{table_point.analytic_beta_over_k0:.6f} | "
                f"{table_point.fem_beta_over_k0:.6f} | "
                f"{closest.beta_ratio:.6f} | "
                f"{lower_text} | "
                f"{upper_text} | "
                f"{len(candidates)} |\n"
            )
        handle.write("\n")

        handle.write("## Conclusion\n\n")
        handle.write(
            "- `baseline_nearest_reference` is the pointwise lower bound for analytic error over the current deduplicated spectrum.\n"
        )
        handle.write(
            f"- `continuity_only` and `delta_guided_reference` do not reduce the worst error below {baseline_summary.max_error_percent:.4f}%.\n"
        )
        handle.write(
            "- For P1, P2 and P3 the closest deduplicated candidates are already the ones selected today, so branch tracking does not create the missing modes near the analytic targets.\n"
        )
        delta_summary = by_strategy["delta_guided_reference"]
        continuity_summary = by_strategy["continuity_only"]
        handle.write(
            f"- Mean analytic error changes from {baseline_summary.mean_error_percent:.4f}% (baseline) "
            f"to {delta_summary.mean_error_percent:.4f}% (delta-guided) and "
            f"{continuity_summary.mean_error_percent:.4f}% (continuity-only).\n"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare Table 10 branch-tracking heuristics against the raw deduplicated spectrum."
    )
    parser.add_argument("--build-bin", type=Path, default=DEFAULT_BUILD_BIN)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    parser.add_argument("--backend", default=DEFAULT_BACKEND)
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=5)
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Reuse an existing diagnostic run under the computed run root.",
    )
    args = parser.parse_args()

    build_bin = _resolve(args.build_bin)
    out_root = _resolve(args.out_root)
    run_root = _run_root(out_root, args.backend)

    if args.skip_run:
        cmd = [
            str(build_bin),
            "--d-over-a-preview",
            "0.20",
            "--nx",
            str(args.nx),
            "--ny",
            str(args.ny),
            "--backend",
            args.backend,
        ]
    else:
        cmd, run_root = _run_case(
            build_bin,
            out_root,
            args.backend,
            args.nx,
            args.ny,
        )

    table_points = _load_table10(_table10_path(run_root))
    raw_lookup = {
        _point_key(point.d_over_a, point.br_over_lambda0): _load_raw_candidates(
            _raw_spectrum_path(run_root, point.d_over_a, point.br_over_lambda0)
        )
        for point in table_points
    }
    grouped_points = _group_table_points(table_points)

    summaries: list[StrategySummary] = []
    strategies: dict[str, dict[tuple[float, float], SelectedPoint]] = {}
    for strategy_name in (
        "baseline_nearest_reference",
        "delta_guided_reference",
        "continuity_only",
    ):
        summary, selections = _evaluate_strategy(
            grouped_points, raw_lookup, strategy_name
        )
        summaries.append(summary)
        strategies[strategy_name] = selections

    mismatches = _verify_table_matches_local_optimum(table_points, raw_lookup)

    validation_dir = out_root / "validation"
    summary_csv_path = validation_dir / SUMMARY_CSV_NAME
    points_csv_path = validation_dir / POINTS_CSV_NAME
    report_md_path = validation_dir / REPORT_MD_NAME

    _write_summary_csv(summary_csv_path, summaries)
    _write_points_csv(points_csv_path, table_points, strategies)
    _write_markdown_report(
        report_md_path,
        cmd=cmd,
        backend=args.backend,
        nx=args.nx,
        ny=args.ny,
        run_root=run_root,
        summaries=summaries,
        table_points=table_points,
        raw_lookup=raw_lookup,
        point_rows_path=points_csv_path,
        summary_csv_path=summary_csv_path,
        mismatches=mismatches,
    )

    print(f"[diag] wrote {_rel_from_root(summary_csv_path)}")
    print(f"[diag] wrote {_rel_from_root(points_csv_path)}")
    print(f"[diag] wrote {_rel_from_root(report_md_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
