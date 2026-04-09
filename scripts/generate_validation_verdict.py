#!/usr/bin/env python3
"""
Generate a scientific validation verdict from the master validation index and a threshold policy.

Default outputs:
- out/validation/validation_verdict.csv
- out/validation/VALIDATION_VERDICT.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POLICY_KEY_FIELDS = ("dimension", "section", "article_ref")


def _resolve(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def _rel_from_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path)


def _to_float(raw: str) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _fmt_float(raw: str) -> str:
    value = _to_float(raw)
    return "" if value is None else f"{value:.12g}"


def _split_pipe(raw: str) -> list[str]:
    return [part.strip() for part in str(raw or "").split("|") if part.strip()]


def _read_csv_required(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"{label} not found: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"{label} is empty: {path}")
    return rows


def _key(row: dict[str, str]) -> tuple[str, str, str]:
    return tuple(str(row.get(field, "")).strip() for field in POLICY_KEY_FIELDS)


def _validate_threshold_pair(row: dict[str, str], lower_field: str, upper_field: str) -> None:
    lower = _to_float(row.get(lower_field, ""))
    upper = _to_float(row.get(upper_field, ""))
    if lower is None and upper is None:
        return
    if lower is None or upper is None:
        raise SystemExit(
            f"Policy row {_key(row)} must define both {lower_field} and {upper_field}, or leave both blank."
        )
    if lower > upper:
        raise SystemExit(f"Policy row {_key(row)} has {lower_field} > {upper_field}.")


def _load_policy(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    rows = _read_csv_required(path, "Validation policy CSV")
    required_fields = {
        "policy_name",
        "dimension",
        "section",
        "article_ref",
        "required_status",
        "required_solvers",
        "primary_pass_max_pct",
        "primary_warn_max_pct",
        "secondary_pass_max_pct",
        "secondary_warn_max_pct",
        "notes",
    }
    missing = sorted(required_fields - set(rows[0].keys()))
    if missing:
        raise SystemExit(f"Validation policy CSV missing columns: {missing}")

    out: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = _key(row)
        if key in out:
            raise SystemExit(f"Duplicate policy row for key {key}")
        _validate_threshold_pair(row, "primary_pass_max_pct", "primary_warn_max_pct")
        _validate_threshold_pair(row, "secondary_pass_max_pct", "secondary_warn_max_pct")
        out[key] = row
    return out


def _load_master_index(path: Path) -> list[dict[str, str]]:
    rows = _read_csv_required(path, "Master validation index CSV")
    required_fields = {
        "dimension",
        "section",
        "article_ref",
        "family",
        "executables",
        "source_index_file",
        "aggregate_files",
        "coverage_key",
        "primary_metric",
        "secondary_metric",
        "status",
        "data_rows",
        "summary_rows",
        "solvers_present",
        "backends_present",
        "max_err_primary_pct",
        "max_err_secondary_pct",
        "notes",
    }
    missing = sorted(required_fields - set(rows[0].keys()))
    if missing:
        raise SystemExit(f"Master validation index CSV missing columns: {missing}")
    return rows


def _validate_coverage(
    master_rows: list[dict[str, str]],
    policy_map: dict[tuple[str, str, str], dict[str, str]],
) -> None:
    master_keys = {_key(row) for row in master_rows}
    policy_keys = set(policy_map.keys())
    missing_policy = sorted(master_keys - policy_keys)
    extra_policy = sorted(policy_keys - master_keys)
    if missing_policy or extra_policy:
        chunks: list[str] = []
        if missing_policy:
            chunks.append(f"missing policy rows for {missing_policy}")
        if extra_policy:
            chunks.append(f"policy rows without master index entry {extra_policy}")
        raise SystemExit("Validation policy coverage mismatch: " + "; ".join(chunks))


def _evaluate_metric(
    *,
    label: str,
    observed_raw: str,
    pass_raw: str,
    warn_raw: str,
) -> tuple[list[str], list[str], list[str]]:
    fail_reasons: list[str] = []
    warn_reasons: list[str] = []
    pass_notes: list[str] = []

    pass_limit = _to_float(pass_raw)
    warn_limit = _to_float(warn_raw)
    if pass_limit is None and warn_limit is None:
        return fail_reasons, warn_reasons, pass_notes

    observed = _to_float(observed_raw)
    if observed is None:
        fail_reasons.append(f"{label} metric missing")
        return fail_reasons, warn_reasons, pass_notes

    if observed > warn_limit:
        fail_reasons.append(f"{label}={observed:.12g} > warn {warn_limit:.12g}")
    elif observed > pass_limit:
        warn_reasons.append(f"{label}={observed:.12g} > pass {pass_limit:.12g} and <= warn {warn_limit:.12g}")
    else:
        pass_notes.append(f"{label}={observed:.12g} <= pass {pass_limit:.12g}")

    return fail_reasons, warn_reasons, pass_notes


def _evaluate_row(
    index_row: dict[str, str],
    policy_row: dict[str, str],
    *,
    policy_path: Path,
) -> dict[str, str]:
    status = str(index_row.get("status", "")).strip()
    required_status = str(policy_row.get("required_status", "")).strip() or "present"

    notes: list[str] = []
    if index_row.get("notes", "").strip():
        notes.append(index_row["notes"].strip())
    if policy_row.get("notes", "").strip():
        notes.append(policy_row["notes"].strip())

    base = {
        "dimension": index_row.get("dimension", ""),
        "section": index_row.get("section", ""),
        "article_ref": index_row.get("article_ref", ""),
        "family": index_row.get("family", ""),
        "executables": index_row.get("executables", ""),
        "source_index_file": index_row.get("source_index_file", ""),
        "aggregate_files": index_row.get("aggregate_files", ""),
        "coverage_key": index_row.get("coverage_key", ""),
        "primary_metric": index_row.get("primary_metric", ""),
        "secondary_metric": index_row.get("secondary_metric", ""),
        "status": status,
        "policy_name": policy_row.get("policy_name", ""),
        "policy_file": _rel_from_root(policy_path),
        "primary_pass_max_pct": _fmt_float(policy_row.get("primary_pass_max_pct", "")),
        "primary_warn_max_pct": _fmt_float(policy_row.get("primary_warn_max_pct", "")),
        "secondary_pass_max_pct": _fmt_float(policy_row.get("secondary_pass_max_pct", "")),
        "secondary_warn_max_pct": _fmt_float(policy_row.get("secondary_warn_max_pct", "")),
        "required_status": required_status,
        "required_solvers": "|".join(_split_pipe(policy_row.get("required_solvers", ""))),
        "observed_solvers": "|".join(_split_pipe(index_row.get("solvers_present", ""))),
        "backends_present": index_row.get("backends_present", ""),
        "max_err_primary_pct": _fmt_float(index_row.get("max_err_primary_pct", "")),
        "max_err_secondary_pct": _fmt_float(index_row.get("max_err_secondary_pct", "")),
        "data_rows": index_row.get("data_rows", ""),
        "summary_rows": index_row.get("summary_rows", ""),
        "notes": " | ".join(notes),
    }

    if status != required_status:
        return {
            **base,
            "scientific_verdict": "MISSING",
            "rule_applied": f"status={status or '<blank>'} != required_status={required_status}",
        }

    fail_reasons: list[str] = []
    warn_reasons: list[str] = []
    pass_notes: list[str] = []

    p_fail, p_warn, p_pass = _evaluate_metric(
        label="primary",
        observed_raw=index_row.get("max_err_primary_pct", ""),
        pass_raw=policy_row.get("primary_pass_max_pct", ""),
        warn_raw=policy_row.get("primary_warn_max_pct", ""),
    )
    fail_reasons.extend(p_fail)
    warn_reasons.extend(p_warn)
    pass_notes.extend(p_pass)

    s_fail, s_warn, s_pass = _evaluate_metric(
        label="secondary",
        observed_raw=index_row.get("max_err_secondary_pct", ""),
        pass_raw=policy_row.get("secondary_pass_max_pct", ""),
        warn_raw=policy_row.get("secondary_warn_max_pct", ""),
    )
    fail_reasons.extend(s_fail)
    warn_reasons.extend(s_warn)
    pass_notes.extend(s_pass)

    required_solvers = set(_split_pipe(policy_row.get("required_solvers", "")))
    observed_solvers = set(_split_pipe(index_row.get("solvers_present", "")))
    if required_solvers:
        missing_solvers = sorted(required_solvers - observed_solvers)
        if missing_solvers:
            warn_reasons.append(f"missing solvers={'|'.join(missing_solvers)}")
        else:
            pass_notes.append(f"solvers cover {'|'.join(sorted(required_solvers))}")

    if fail_reasons:
        verdict = "FAIL"
        rule = "; ".join(fail_reasons + warn_reasons)
    elif warn_reasons:
        verdict = "WARN"
        rule = "; ".join(warn_reasons + pass_notes)
    else:
        verdict = "PASS"
        rule = "; ".join(pass_notes or ["status present and metrics within pass limits"])

    return {
        **base,
        "scientific_verdict": verdict,
        "rule_applied": rule,
    }


def build_verdict_rows(master_index: Path, policy_csv: Path) -> list[dict[str, str]]:
    master_rows = _load_master_index(master_index)
    policy_map = _load_policy(policy_csv)
    _validate_coverage(master_rows, policy_map)

    out_rows: list[dict[str, str]] = []
    for row in master_rows:
        out_rows.append(_evaluate_row(row, policy_map[_key(row)], policy_path=policy_csv))
    return out_rows


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "dimension",
        "section",
        "article_ref",
        "family",
        "executables",
        "source_index_file",
        "aggregate_files",
        "coverage_key",
        "status",
        "scientific_verdict",
        "rule_applied",
        "policy_name",
        "policy_file",
        "primary_metric",
        "max_err_primary_pct",
        "primary_pass_max_pct",
        "primary_warn_max_pct",
        "secondary_metric",
        "max_err_secondary_pct",
        "secondary_pass_max_pct",
        "secondary_warn_max_pct",
        "required_status",
        "required_solvers",
        "observed_solvers",
        "backends_present",
        "data_rows",
        "summary_rows",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_md(path: Path, rows: list[dict[str, str]]) -> None:
    now = dt.datetime.now().isoformat(timespec="seconds")
    counts = Counter(row["scientific_verdict"] for row in rows)
    policy_file = rows[0]["policy_file"] if rows else ""
    source_files: list[str] = []
    for row in rows:
        source = row["source_index_file"]
        if source and source not in source_files:
            source_files.append(source)

    lines: list[str] = []
    lines.append("# Veredito de Validacao Cientifica")
    lines.append("")
    lines.append(f"Gerado em: `{now}`")
    lines.append("")
    lines.append("Este veredito aplica uma politica CSV auditavel sobre o indice mestre de validacao.")
    lines.append("Cada item do artigo recebe `PASS`, `WARN`, `FAIL` ou `MISSING` sem recalcular solvers.")
    lines.append("")
    lines.append("## Resumo")
    lines.append("")
    lines.append(f"- PASS: `{counts.get('PASS', 0)}`")
    lines.append(f"- WARN: `{counts.get('WARN', 0)}`")
    lines.append(f"- FAIL: `{counts.get('FAIL', 0)}`")
    lines.append(f"- MISSING: `{counts.get('MISSING', 0)}`")
    lines.append("")
    lines.append("| Dim | Secao | Tabela/Figura | Familia | Veredito | Status | Max err prim. (%) | Max err sec. (%) | Regra aplicada |")
    lines.append("|---|---|---|---|---|---|---:|---:|---|")
    for row in rows:
        lines.append(
            f"| `{row['dimension']}` | `{row['section']}` | `{row['article_ref']}` | `{row['family']}` | "
            f"`{row['scientific_verdict']}` | `{row['status']}` | `{row['max_err_primary_pct']}` | "
            f"`{row['max_err_secondary_pct']}` | `{row['rule_applied']}` |"
        )
    lines.append("")
    lines.append("## Fontes")
    lines.append("")
    if policy_file:
        lines.append(f"- Politica de thresholds: `{policy_file}`")
    for source in source_files:
        lines.append(f"- Indice de entrada: `{source}`")
    lines.append("")
    lines.append("## Observacoes")
    lines.append("")
    non_pass = [row for row in rows if row["scientific_verdict"] != "PASS"]
    if non_pass:
        for row in non_pass:
            lines.append(
                f"- `{row['dimension']} {row['article_ref']}` ficou `{row['scientific_verdict']}` "
                f"porque `{row['rule_applied']}`."
            )
    else:
        lines.append("- Todos os itens avaliados ficaram dentro das faixas de aprovacao desta politica.")
    lines.append("")
    lines.append("- CSV maquina-legivel: `[validation_verdict.csv](validation_verdict.csv)`")

    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a scientific validation verdict from the master index and a threshold policy.")
    parser.add_argument(
        "--master-index",
        type=Path,
        default=Path("out/validation/validation_master_index.csv"),
        help="Input master validation index CSV.",
    )
    parser.add_argument(
        "--policy-csv",
        type=Path,
        default=Path("docs/validation_thresholds.csv"),
        help="Threshold policy CSV.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("out/validation/validation_verdict.csv"),
        help="Output CSV path for the scientific validation verdict.",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("out/validation/VALIDATION_VERDICT.md"),
        help="Output Markdown path for the scientific validation verdict.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    master_index = _resolve(args.master_index)
    policy_csv = _resolve(args.policy_csv)
    out_csv = _resolve(args.out_csv)
    out_md = _resolve(args.out_md)

    rows = build_verdict_rows(master_index, policy_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_csv, rows)
    _write_md(out_md, rows)

    print(f"Saved: {_rel_from_root(out_csv)}")
    print(f"Saved: {_rel_from_root(out_md)}")


if __name__ == "__main__":
    main()
