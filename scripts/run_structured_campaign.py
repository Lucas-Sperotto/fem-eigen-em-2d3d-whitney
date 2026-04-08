#!/usr/bin/env python3
"""
Create and execute a structured campaign tree under out/.

Implemented scope:
- resolve selected cases and backends;
- generate or sanitize the execution name;
- create the execution root inside out/;
- scaffold backend/section/case/configuration directories;
- isolate TP3485_OUT_DIR per configuration;
- execute each configuration and save stdout/stderr logs;
- normalize 2D field VTK files into vtk/;
- generate PNG images and CSV summaries inside each configuration folder;
- write manifest.json, RUN.md, CASE.json, config.json and per-config result.json.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import time
import traceback
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

from benchmark_backends import parse_timing
from plot_vtk_quiver import _parse_dims as parse_dims_2d
from plot_vtk_quiver import _parse_mode_table as parse_mode_table_2d
from plot_vtk_quiver import _plot_all_images as plot_all_field_images
from plot_vtk_quiver import _write_csv as write_mode_summary_csv
from validate_2d_22 import parse_first_kc_block
from validate_2d_22 import parse_helmvec2_table
from validate_2d_22 import parse_helmvec3_table9
from validate_2d_22 import parse_helmvec3_table10
from validate_2d_22 import parse_mixed_rect_table
from validate_3d_31 import parse_mode_table as parse_mode_table_3d


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_BASE = ROOT / "out"
DEFAULT_BUILD = ROOT / "build"
ARTIFACT_SUBDIRS = ("vtk", "img", "csv", "logs", "meta")
SECTION_ORDER = {"2.1": 1, "2.2.1": 2, "2.2.2": 3, "2.2.3": 4, "2.2.4": 5, "3.1": 6}
FIELD_CASE_IDS = {
    "helm10_rect",
    "helm10_circle",
    "helm10_coax",
    "edge_rect",
    "edge_circle",
    "edge_coax",
}
VALIDATION_2D_CASE_IDS = {
    "mixed_rect",
    "mixed_circle",
    "mixed_coax",
    "helmvec2_rect",
    "helmvec3_fig12_rect",
    "helmvec3_fig13_rect",
}


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    section: str
    family: str
    program: str
    label: str


@dataclass(frozen=True)
class ConfigSpec:
    config_id: str
    label: str
    purposes: tuple[str, ...]
    params: dict[str, object]


@dataclass(frozen=True)
class ConfigEntry:
    backend: str
    backend_dir: str
    case: CaseSpec
    config: ConfigSpec

    @property
    def rel_dir(self) -> Path:
        return Path(self.backend_dir) / self.case.section / self.case.case_id / self.config.config_id


@dataclass(frozen=True)
class CommandResult:
    step_name: str
    command: list[str]
    wall_ms: float
    returncode: int
    stdout_path: Path
    stderr_path: Path
    stdout_text: str
    stderr_text: str
    timing: dict[str, float]


CASE_SPECS: list[CaseSpec] = [
    CaseSpec("helm10_rect", "2.1", "2d", "helm10_rect", "Scalar rectangular guide"),
    CaseSpec("helm10_circle", "2.1", "2d", "helm10_circle", "Scalar circular guide"),
    CaseSpec("helm10_coax", "2.1", "2d", "helm10_coax", "Scalar coaxial guide"),
    CaseSpec("edge_rect", "2.2.1", "2d", "edge_rect", "Edge rectangular guide"),
    CaseSpec("edge_circle", "2.2.1", "2d", "edge_circle", "Edge circular guide"),
    CaseSpec("edge_coax", "2.2.1", "2d", "edge_coax", "Edge coaxial guide"),
    CaseSpec("mixed_rect", "2.2.2", "2d", "mixed_rect", "Mixed rectangular cutoff validation"),
    CaseSpec("mixed_circle", "2.2.2", "2d", "mixed_circle", "Mixed circular cutoff snapshots"),
    CaseSpec("mixed_coax", "2.2.2", "2d", "mixed_coax", "Mixed coaxial cutoff snapshots"),
    CaseSpec("helmvec2_rect", "2.2.3", "2d", "helmvec2_rect", "Figure 11 / Table 8 validation"),
    CaseSpec("helmvec3_fig12_rect", "2.2.4", "2d", "helmvec3_fig12_rect", "Figure 12 / Table 9 validation"),
    CaseSpec("helmvec3_fig13_rect", "2.2.4", "2d", "helmvec3_fig13_rect", "Figure 13 / Table 10 validation"),
    CaseSpec("fem3d0_air", "3.1", "3d", "fem3d0_air", "3D air cavity with fem3d0"),
    CaseSpec("fem3d1_air", "3.1", "3d", "fem3d1_air", "3D air cavity with fem3d1"),
    CaseSpec("fem3d0_half", "3.1", "3d", "fem3d0_half", "3D half-filled cavity with fem3d0"),
    CaseSpec("fem3d1_half", "3.1", "3d", "fem3d1_half", "3D half-filled cavity with fem3d1"),
    CaseSpec("fem3d0_cyl", "3.1", "3d", "fem3d0_cyl", "3D cylindrical cavity with fem3d0"),
    CaseSpec("fem3d1_cyl", "3.1", "3d", "fem3d1_cyl", "3D cylindrical cavity with fem3d1"),
    CaseSpec("fem3d0_sphere", "3.1", "3d", "fem3d0_sphere", "3D spherical cavity with fem3d0"),
    CaseSpec("fem3d1_sphere", "3.1", "3d", "fem3d1_sphere", "3D spherical cavity with fem3d1"),
]

CASE_BY_ID = {spec.case_id: spec for spec in CASE_SPECS}
CASE_ORDER = {spec.case_id: idx for idx, spec in enumerate(CASE_SPECS)}


def normalize_token(value: str) -> str:
    token = value.strip().lower()
    token = token.replace(".", "")
    token = token.replace("_", "")
    token = token.replace("-", "")
    token = token.replace("/", "")
    token = token.replace(":", "")
    token = token.replace(" ", "")
    return token


def backend_dirname(backend: str) -> str:
    return backend.replace("-", "_")


def case_ids_for_family(family: str) -> list[str]:
    return [spec.case_id for spec in CASE_SPECS if spec.family == family]


def case_ids_for_section(section: str) -> list[str]:
    return [spec.case_id for spec in CASE_SPECS if spec.section == section]


ALIAS_GROUPS: dict[str, list[str]] = {
    "all": [spec.case_id for spec in CASE_SPECS],
    "2": case_ids_for_family("2d"),
    "2d": case_ids_for_family("2d"),
    "sec2": case_ids_for_family("2d"),
    "section2": case_ids_for_family("2d"),
    "21": case_ids_for_section("2.1"),
    "sec21": case_ids_for_section("2.1"),
    "section21": case_ids_for_section("2.1"),
    "scalar": case_ids_for_section("2.1"),
    "helm10": case_ids_for_section("2.1"),
    "table1": case_ids_for_section("2.1"),
    "table2": case_ids_for_section("2.1"),
    "table3": case_ids_for_section("2.1"),
    "tabela1": case_ids_for_section("2.1"),
    "tabela2": case_ids_for_section("2.1"),
    "tabela3": case_ids_for_section("2.1"),
    "221": case_ids_for_section("2.2.1"),
    "sec221": case_ids_for_section("2.2.1"),
    "section221": case_ids_for_section("2.2.1"),
    "edge": case_ids_for_section("2.2.1"),
    "222": case_ids_for_section("2.2.2"),
    "sec222": case_ids_for_section("2.2.2"),
    "section222": case_ids_for_section("2.2.2"),
    "mixed": case_ids_for_section("2.2.2"),
    "coupled": case_ids_for_section("2.2.2"),
    "acoplado": case_ids_for_section("2.2.2"),
    "223": case_ids_for_section("2.2.3"),
    "sec223": case_ids_for_section("2.2.3"),
    "section223": case_ids_for_section("2.2.3"),
    "fig11": case_ids_for_section("2.2.3"),
    "figure11": case_ids_for_section("2.2.3"),
    "table8": case_ids_for_section("2.2.3"),
    "tabela8": case_ids_for_section("2.2.3"),
    "helmvec2": case_ids_for_section("2.2.3"),
    "224": case_ids_for_section("2.2.4"),
    "sec224": case_ids_for_section("2.2.4"),
    "section224": case_ids_for_section("2.2.4"),
    "fig12": case_ids_for_section("2.2.4"),
    "figure12": case_ids_for_section("2.2.4"),
    "fig13": case_ids_for_section("2.2.4"),
    "figure13": case_ids_for_section("2.2.4"),
    "table9": case_ids_for_section("2.2.4"),
    "table10": case_ids_for_section("2.2.4"),
    "tabela9": case_ids_for_section("2.2.4"),
    "tabela10": case_ids_for_section("2.2.4"),
    "helmvec3": case_ids_for_section("2.2.4"),
    "31": case_ids_for_section("3.1"),
    "sec31": case_ids_for_section("3.1"),
    "section31": case_ids_for_section("3.1"),
    "3d": case_ids_for_section("3.1"),
    "fem3d": case_ids_for_section("3.1"),
    "all3d": case_ids_for_section("3.1"),
    "table12": ["fem3d0_air", "fem3d1_air"],
    "tabela12": ["fem3d0_air", "fem3d1_air"],
    "fig15": ["fem3d0_air", "fem3d1_air"],
    "figure15": ["fem3d0_air", "fem3d1_air"],
    "air": ["fem3d0_air", "fem3d1_air"],
    "table13": ["fem3d0_half", "fem3d1_half"],
    "tabela13": ["fem3d0_half", "fem3d1_half"],
    "fig16": ["fem3d0_half", "fem3d1_half"],
    "figure16": ["fem3d0_half", "fem3d1_half"],
    "half": ["fem3d0_half", "fem3d1_half"],
    "halffilled": ["fem3d0_half", "fem3d1_half"],
    "table14": ["fem3d0_cyl", "fem3d1_cyl"],
    "tabela14": ["fem3d0_cyl", "fem3d1_cyl"],
    "fig17": ["fem3d0_cyl", "fem3d1_cyl"],
    "figure17": ["fem3d0_cyl", "fem3d1_cyl"],
    "cyl": ["fem3d0_cyl", "fem3d1_cyl"],
    "cylinder": ["fem3d0_cyl", "fem3d1_cyl"],
    "cil": ["fem3d0_cyl", "fem3d1_cyl"],
    "cilindro": ["fem3d0_cyl", "fem3d1_cyl"],
    "table15": ["fem3d0_sphere", "fem3d1_sphere"],
    "tabela15": ["fem3d0_sphere", "fem3d1_sphere"],
    "sphere": ["fem3d0_sphere", "fem3d1_sphere"],
    "spherical": ["fem3d0_sphere", "fem3d1_sphere"],
    "esfera": ["fem3d0_sphere", "fem3d1_sphere"],
    "fem3d0": [spec.case_id for spec in CASE_SPECS if spec.case_id.startswith("fem3d0_")],
    "fem3d1": [spec.case_id for spec in CASE_SPECS if spec.case_id.startswith("fem3d1_")],
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Create and execute a structured execution tree inside out/."
    )
    ap.add_argument(
        "--run-name",
        help="Execution folder name under out/. Auto-generated when omitted.",
    )
    ap.add_argument(
        "--out-base",
        type=Path,
        default=DEFAULT_OUT_BASE,
        help="Base directory that contains execution folders (default: out).",
    )
    ap.add_argument(
        "--build-dir",
        type=Path,
        default=DEFAULT_BUILD,
        help="Build directory containing executables.",
    )
    ap.add_argument(
        "--backend",
        choices=["gauss", "closed-form", "both"],
        default="both",
        help="Backends executed for this campaign.",
    )
    ap.add_argument(
        "--case",
        action="append",
        default=[],
        help="Case or section selector. Repeatable. Defaults to all cases.",
    )
    ap.add_argument(
        "--profile",
        choices=["quick", "full"],
        default="quick",
        help="3D profile used to scaffold per-grid configuration folders.",
    )
    ap.add_argument(
        "--article-only",
        action="store_true",
        help="Use only the article baseline configuration for each selected case.",
    )
    ap.add_argument(
        "--node-scale",
        type=float,
        default=1.0,
        help="Approximate nodal-count multiplier relative to each article baseline mesh.",
    )
    ap.add_argument(
        "--mode-export",
        type=int,
        default=8,
        help="Number of 2D modes exported in field-producing configurations.",
    )
    ap.add_argument(
        "--skip-images",
        action="store_true",
        help="Skip PNG generation while still preserving raw VTK/CSV outputs.",
    )
    ap.add_argument(
        "--skip-validate",
        action="store_true",
        help="Skip validation CSV/image generation for sections 2.2.2+ and 3.1.",
    )
    ap.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip CMake configure/build before running selected configurations.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the execution plan without creating files or directories.",
    )
    return ap.parse_args()


def sanitize_run_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    clean = re.sub(r"_+", "_", clean)
    clean = clean.strip("._-")
    if not clean:
        raise SystemExit("Invalid --run-name: expected at least one alphanumeric character.")
    return clean


def auto_run_name(now: dt.datetime) -> str:
    return now.strftime("run_%Y%m%d_%H%M%S")


def unique_auto_run_name(out_base: Path, now: dt.datetime) -> str:
    base = auto_run_name(now)
    candidate = base
    suffix = 1
    while (out_base / candidate).exists():
        suffix += 1
        candidate = f"{base}_{suffix:02d}"
    return candidate


def resolve_backends(raw_backend: str) -> list[str]:
    if raw_backend == "both":
        return ["gauss", "closed-form"]
    return [raw_backend]


def resolve_cases(raw_cases: list[str]) -> list[CaseSpec]:
    if not raw_cases:
        return list(CASE_SPECS)

    selected_ids: set[str] = set()
    unknown: list[str] = []
    normalized_case_ids = {
        normalize_token(spec.case_id): spec.case_id for spec in CASE_SPECS
    }

    for raw in raw_cases:
        token = normalize_token(raw)
        if token in normalized_case_ids:
            selected_ids.add(normalized_case_ids[token])
            continue
        if token in ALIAS_GROUPS:
            selected_ids.update(ALIAS_GROUPS[token])
            continue
        unknown.append(raw)

    if unknown:
        examples = ", ".join(
            ["2.1", "2.2.3", "3.1", "helm10_rect", "fem3d0_air", "table13"]
        )
        raise SystemExit(
            "Unknown --case value(s): "
            + ", ".join(unknown)
            + f". Try values like: {examples}."
        )

    return sorted(
        (CASE_BY_ID[case_id] for case_id in selected_ids),
        key=lambda spec: CASE_ORDER[spec.case_id],
    )


def rel_to_root(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except Exception:
        return str(path.resolve())


def format_number_tag(value: float) -> str:
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return text.replace(".", "p")


def default_grid_3d(case_name: str) -> tuple[int, int, int]:
    grids = {
        "air": (6, 3, 3),
        "half": (5, 5, 4),
        "cyl": (7, 7, 4),
        "sphere": (6, 5, 5),
    }
    return grids[case_name]


def validation_grids(profile: str) -> dict[str, list[tuple[int, int, int]]]:
    quick = {
        "air": [(5, 3, 3), (6, 3, 3)],
        "half": [(4, 4, 3), (5, 5, 4)],
        "cyl": [(6, 6, 3), (7, 7, 4)],
        "sphere": [(5, 4, 4), (6, 5, 5)],
    }
    full = {
        "air": quick["air"] + [(7, 4, 4)],
        "half": quick["half"] + [(6, 6, 5)],
        "cyl": quick["cyl"] + [(8, 8, 5)],
        "sphere": quick["sphere"] + [(7, 6, 6)],
    }
    return full if profile == "full" else quick


def snap_even_at_least(value: float, minimum: int = 4) -> int:
    out = max(minimum, int(round(value)))
    if out % 2:
        out += 1
    return out


def scale_rect_grid_2d(nx: int, ny: int, node_scale: float) -> tuple[int, int]:
    factor = math.sqrt(node_scale)
    return (
        max(1, int(round((nx + 1) * factor - 1))),
        max(1, int(round((ny + 1) * factor - 1))),
    )


def scale_rect_grid_3d(nx: int, ny: int, nz: int, node_scale: float) -> tuple[int, int, int]:
    factor = node_scale ** (1.0 / 3.0)
    return (
        max(1, int(round((nx + 1) * factor - 1))),
        max(1, int(round((ny + 1) * factor - 1))),
        max(1, int(round((nz + 1) * factor - 1))),
    )


def scale_circle_grid(nr: int, nt: int, node_scale: float) -> tuple[int, int]:
    factor = math.sqrt(node_scale)
    return (
        max(1, int(round(nr * factor))),
        snap_even_at_least(nt * factor),
    )


def scale_coax_grid(nr: int, nt: int, node_scale: float) -> tuple[int, int]:
    factor = math.sqrt(node_scale)
    return (
        max(1, int(round((nr + 1) * factor - 1))),
        snap_even_at_least(nt * factor),
    )


def config_id_for_case(case_id: str, params: dict[str, object]) -> str:
    if case_id in {"helm10_rect", "edge_rect"}:
        return f"nx{params['nx']}_ny{params['ny']}_modes{params['mode_export']}"
    if case_id in {"helm10_circle", "helm10_coax", "edge_circle", "edge_coax"}:
        return f"nr{params['nr']}_nt{params['nt']}_modes{params['mode_export']}"
    if case_id == "mixed_rect":
        return f"nx{params['nx']}_ny{params['ny']}"
    if case_id in {"mixed_circle", "mixed_coax"}:
        return f"nr{params['nr']}_nt{params['nt']}"
    if case_id == "helmvec2_rect":
        return f"beta{format_number_tag(float(params['beta']))}_nx{params['nx']}_ny{params['ny']}"
    if case_id == "helmvec3_fig12_rect":
        return f"nx{params['nx']}_ny{params['ny']}"
    if case_id == "helmvec3_fig13_rect":
        return f"d_over_a{format_number_tag(float(params['d_over_a']))}_nx{params['nx']}_ny{params['ny']}"
    if case_id.startswith("fem3d"):
        return f"nx{params['nx']}_ny{params['ny']}_nz{params['nz']}"
    raise SystemExit(f"Missing config id formatter for case: {case_id}")


def scale_config_params(case_id: str, params: dict[str, object], node_scale: float) -> dict[str, object]:
    scaled = dict(params)
    if math.isclose(node_scale, 1.0, rel_tol=0.0, abs_tol=1e-12):
        return scaled

    if case_id in {"helm10_rect", "edge_rect", "mixed_rect", "helmvec2_rect", "helmvec3_fig12_rect", "helmvec3_fig13_rect"}:
        scaled["nx"], scaled["ny"] = scale_rect_grid_2d(int(params["nx"]), int(params["ny"]), node_scale)
        return scaled

    if case_id in {"helm10_circle", "edge_circle", "mixed_circle"}:
        scaled["nr"], scaled["nt"] = scale_circle_grid(int(params["nr"]), int(params["nt"]), node_scale)
        return scaled

    if case_id in {"helm10_coax", "edge_coax", "mixed_coax"}:
        scaled["nr"], scaled["nt"] = scale_coax_grid(int(params["nr"]), int(params["nt"]), node_scale)
        return scaled

    if case_id.startswith("fem3d"):
        scaled["nx"], scaled["ny"], scaled["nz"] = scale_rect_grid_3d(
            int(params["nx"]),
            int(params["ny"]),
            int(params["nz"]),
            node_scale,
        )
        return scaled

    raise SystemExit(f"Missing node scaling rule for case: {case_id}")


def finalize_config_specs(case: CaseSpec, configs: list[ConfigSpec], node_scale: float) -> list[ConfigSpec]:
    if math.isclose(node_scale, 1.0, rel_tol=0.0, abs_tol=1e-12):
        return configs

    out: list[ConfigSpec] = []
    for config in configs:
        scaled_params = scale_config_params(case.case_id, config.params, node_scale)
        out.append(
            ConfigSpec(
                config_id=config_id_for_case(case.case_id, scaled_params),
                label=config.label,
                purposes=config.purposes,
                params=scaled_params,
            )
        )
    return out


def config_specs_for_case(
    case: CaseSpec,
    *,
    mode_export: int,
    profile: str,
    article_only: bool,
    node_scale: float,
) -> list[ConfigSpec]:
    if case.case_id in {"helm10_rect", "edge_rect"}:
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id=f"nx14_ny14_modes{mode_export}",
                label="Primary rectangular 2D run",
                purposes=("direct_default",),
                params={"nx": 14, "ny": 14, "mode_export": mode_export},
            )
        ], node_scale)
    if case.case_id in {"helm10_circle", "helm10_coax", "edge_circle", "edge_coax"}:
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id=f"nr10_nt48_modes{mode_export}",
                label="Primary polar 2D run",
                purposes=("direct_default",),
                params={"nr": 10, "nt": 48, "mode_export": mode_export},
            )
        ], node_scale)
    if case.case_id == "mixed_rect":
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id="nx12_ny6",
                label="Mixed rectangular cutoff run",
                purposes=("direct_default",),
                params={"nx": 12, "ny": 6},
            )
        ], node_scale)
    if case.case_id in {"mixed_circle", "mixed_coax"}:
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id="nr10_nt48",
                label="Mixed polar cutoff run",
                purposes=("direct_default",),
                params={"nr": 10, "nt": 48},
            )
        ], node_scale)
    if case.case_id == "helmvec2_rect":
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id="beta10_nx6_ny6",
                label="Figure 11 / Table 8 default run",
                purposes=("direct_default", "validation_default"),
                params={"beta": 10.0, "nx": 6, "ny": 6, "legacy_debug": 0},
            )
        ], node_scale)
    if case.case_id == "helmvec3_fig12_rect":
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id="nx10_ny5",
                label="Figure 12 / Table 9 default run",
                purposes=("direct_default", "validation_default"),
                params={"nx": 10, "ny": 5, "legacy_debug": 0},
            )
        ], node_scale)
    if case.case_id == "helmvec3_fig13_rect":
        return finalize_config_specs(case, [
            ConfigSpec(
                config_id=f"d_over_a{format_number_tag(0.20)}_nx10_ny5",
                label="Figure 13 / Table 10 default run",
                purposes=("direct_default", "validation_default"),
                params={"d_over_a": 0.20, "nx": 10, "ny": 5, "legacy_debug": 0},
            )
        ], node_scale)
    if case.case_id.startswith("fem3d0_") or case.case_id.startswith("fem3d1_"):
        case_name = case.case_id.split("_", 1)[1]
        by_grid: dict[tuple[int, int, int], list[str]] = {}
        default_grid = default_grid_3d(case_name)
        by_grid.setdefault(default_grid, []).append("direct_default")
        if not article_only:
            for idx, grid in enumerate(validation_grids(profile)[case_name], start=1):
                by_grid.setdefault(grid, []).append(f"validation_{profile}_{idx:02d}")
        out: list[ConfigSpec] = []
        for grid in sorted(by_grid.keys()):
            nx, ny, nz = grid
            out.append(
                ConfigSpec(
                    config_id=f"nx{nx}_ny{ny}_nz{nz}",
                    label="3D cavity grid configuration",
                    purposes=tuple(by_grid[grid]),
                    params={"nx": nx, "ny": ny, "nz": nz, "profile": profile},
                )
            )
        return finalize_config_specs(case, out, node_scale)
    raise SystemExit(f"Missing configuration scaffold for case: {case.case_id}")


def build_entries(
    *,
    backends: list[str],
    selected_cases: list[CaseSpec],
    mode_export: int,
    profile: str,
    article_only: bool,
    node_scale: float,
) -> list[ConfigEntry]:
    entries: list[ConfigEntry] = []
    for backend in backends:
        backend_dir = backend_dirname(backend)
        for case in selected_cases:
            for config in config_specs_for_case(
                case,
                mode_export=mode_export,
                profile=profile,
                article_only=article_only,
                node_scale=node_scale,
            ):
                entries.append(
                    ConfigEntry(
                        backend=backend,
                        backend_dir=backend_dir,
                        case=case,
                        config=config,
                    )
                )
    return entries


def build_manifest(
    *,
    args: argparse.Namespace,
    now: dt.datetime,
    run_name: str,
    run_root: Path,
    backends: list[str],
    selected_cases: list[CaseSpec],
    entries: list[ConfigEntry],
    execution_summary: dict[str, object] | None = None,
) -> dict[str, object]:
    sections = sorted(
        {spec.section for spec in selected_cases},
        key=lambda value: SECTION_ORDER.get(value, 999),
    )
    manifest = {
        "schema_version": 3,
        "stage": "task_3_execute_structured_tree",
        "created_at": now.isoformat(timespec="seconds"),
        "run_name": run_name,
        "run_root": str(run_root),
        "run_root_from_repo": rel_to_root(run_root),
        "out_base": str(args.out_base.resolve()),
        "build_dir": str(args.build_dir.resolve()),
        "build_dir_exists": args.build_dir.resolve().exists(),
        "requested_backend": args.backend,
        "backends": [
            {
                "backend": backend,
                "dir_name": backend_dirname(backend),
                "path": str(run_root / backend_dirname(backend)),
            }
            for backend in backends
        ],
        "requested_cases": args.case if args.case else ["all"],
        "sections": sections,
        "profile": args.profile,
        "article_only": args.article_only,
        "node_scale": args.node_scale,
        "mode_export": args.mode_export,
        "skip_images": args.skip_images,
        "skip_validate": args.skip_validate,
        "skip_build": args.skip_build,
        "dry_run": args.dry_run,
        "artifact_layout": {
            "levels": ["backend", "section", "case", "configuration"],
            "leaf_subdirs": list(ARTIFACT_SUBDIRS),
            "tp3485_out_dir_per_config": "meta/tp3485_out",
        },
        "case_count": len(selected_cases),
        "config_count_per_backend": len(entries) // max(1, len(backends)),
        "config_count_total": len(entries),
        "cases": [
            {
                "case_id": spec.case_id,
                "section": spec.section,
                "family": spec.family,
                "program": spec.program,
                "label": spec.label,
                "config_count": len(
                    config_specs_for_case(
                        spec,
                        mode_export=args.mode_export,
                        profile=args.profile,
                        article_only=args.article_only,
                        node_scale=args.node_scale,
                    )
                ),
            }
            for spec in selected_cases
        ],
        "entries": [
            {
                "backend": entry.backend,
                "backend_dir": entry.backend_dir,
                "section": entry.case.section,
                "case_id": entry.case.case_id,
                "program": entry.case.program,
                "config_id": entry.config.config_id,
                "purposes": list(entry.config.purposes),
                "params": entry.config.params,
                "path": str(run_root / entry.rel_dir),
                "path_from_repo": rel_to_root(run_root / entry.rel_dir),
            }
            for entry in entries
        ],
    }
    if execution_summary is not None:
        manifest["execution_summary"] = execution_summary
    return manifest


def render_run_md(
    *,
    now: dt.datetime,
    run_name: str,
    run_root: Path,
    build_dir: Path,
    backends: list[str],
    selected_cases: list[CaseSpec],
    entries: list[ConfigEntry],
    args: argparse.Namespace,
    execution_summary: dict[str, object] | None = None,
    results: list[dict[str, object]] | None = None,
) -> str:
    lines: list[str] = []
    lines.append("# Structured Execution")
    lines.append("")
    lines.append(f"Run name: `{run_name}`")
    lines.append(f"Generated at: `{now.isoformat(timespec='seconds')}`")
    lines.append(f"Run root: `{rel_to_root(run_root)}`")
    lines.append("")
    lines.append("This run creates one root per execution, separates backends at the top level,")
    lines.append("and stores per-configuration artifacts in `vtk/`, `img/`, `csv/`, `logs/` and `meta/`.")
    lines.append("")
    lines.append("## Base Configuration")
    lines.append("")
    lines.append(f"- Build dir: `{rel_to_root(build_dir)}`")
    lines.append(f"- Build dir exists: `{'yes' if build_dir.exists() else 'no'}`")
    lines.append(f"- Backends: `{', '.join(backends)}`")
    lines.append(f"- Requested cases: `{', '.join(args.case) if args.case else 'all'}`")
    lines.append(f"- 3D profile: `{args.profile}`")
    lines.append(f"- Article-only configs: `{'yes' if args.article_only else 'no'}`")
    lines.append(f"- Node scale: `{args.node_scale:g}x`")
    lines.append(f"- 2D mode export target: `{args.mode_export}`")
    lines.append(f"- Skip images: `{'yes' if args.skip_images else 'no'}`")
    lines.append(f"- Skip validations: `{'yes' if args.skip_validate else 'no'}`")
    lines.append(f"- Skip build: `{'yes' if args.skip_build else 'no'}`")
    lines.append(f"- Manifest: [manifest.json](manifest.json)")
    lines.append("")
    lines.append("## Artifact Layout")
    lines.append("")
    lines.append("```text")
    lines.append(f"{run_name}/")
    lines.append("  <backend>/")
    lines.append("    <section>/")
    lines.append("      <case>/")
    lines.append("        <configuration>/")
    for leaf in ARTIFACT_SUBDIRS:
        lines.append(f"          {leaf}/")
    lines.append("          config.json")
    lines.append("          meta/result.json")
    lines.append("```")
    lines.append("")
    lines.append("## Selected Cases")
    lines.append("")
    lines.append("| Case ID | Section | Family | Program | Description | Configs |")
    lines.append("|---|---|---|---|---|---:|")
    for spec in selected_cases:
        cfg_count = len(
            config_specs_for_case(
                spec,
                mode_export=args.mode_export,
                profile=args.profile,
                article_only=args.article_only,
                node_scale=args.node_scale,
            )
        )
        lines.append(
            f"| `{spec.case_id}` | `{spec.section}` | `{spec.family}` | "
            f"`{spec.program}` | {spec.label} | {cfg_count} |"
        )
    lines.append("")
    lines.append("## Configuration Roots")
    lines.append("")
    lines.append("| Backend | Section | Case | Config | Purposes | Path |")
    lines.append("|---|---|---|---|---|---|")
    for entry in entries:
        rel_dir = rel_to_root(run_root / entry.rel_dir)
        lines.append(
            f"| `{entry.backend}` | `{entry.case.section}` | `{entry.case.case_id}` | "
            f"`{entry.config.config_id}` | `{', '.join(entry.config.purposes)}` | "
            f"[{rel_dir}]({entry.rel_dir.as_posix()}/) |"
        )
    if execution_summary is not None:
        lines.append("")
        lines.append("## Execution Summary")
        lines.append("")
        lines.append(f"- Executed configurations: `{execution_summary['executed_configurations']}`")
        lines.append(f"- Total VTK files: `{execution_summary['vtk_files']}`")
        lines.append(f"- Total PNG files: `{execution_summary['png_files']}`")
        lines.append(f"- Total CSV files: `{execution_summary['csv_files']}`")
        lines.append(f"- Total LOG files: `{execution_summary['log_files']}`")
    if results is not None:
        result_by_key = {
            (str(result["backend"]), str(result["case_id"]), str(result["config_id"])): result
            for result in results
        }
        grouped_entries: dict[tuple[str, str], list[ConfigEntry]] = defaultdict(list)
        for entry in entries:
            grouped_entries[(entry.case.section, entry.case.case_id)].append(entry)

        lines.append("")
        lines.append("## Index by Section/Case")
        lines.append("")
        for (section, case_id), case_entries in sorted(
            grouped_entries.items(),
            key=lambda item: (SECTION_ORDER.get(item[0][0], 999), CASE_ORDER[item[0][1]]),
        ):
            sample = case_entries[0]
            backend_links: list[str] = []
            seen_backend_dirs: set[str] = set()
            for entry in case_entries:
                if entry.backend_dir in seen_backend_dirs:
                    continue
                seen_backend_dirs.add(entry.backend_dir)
                backend_links.append(
                    f"`{entry.backend}`: {run_rel_link(run_root, run_root / entry.backend_dir / section / case_id, entry.backend_dir)}"
                )
            lines.append(
                f"- `{section} / {case_id}`: {sample.case.label}. "
                + " | ".join(backend_links)
            )

        lines.append("")
        lines.append("## Results by Section/Case")
        lines.append("")
        for (section, case_id), case_entries in sorted(
            grouped_entries.items(),
            key=lambda item: (SECTION_ORDER.get(item[0][0], 999), CASE_ORDER[item[0][1]]),
        ):
            sample = case_entries[0]
            lines.append(f"### {section} / {case_id}")
            lines.append("")
            lines.append(sample.case.label)
            lines.append("")
            lines.append("| Backend | Config | Purposes | Mesh | Nodes | Elements | DOF | wall_ms | total_ms | CSVs | Images | Logs | Meta |")
            lines.append("|---|---|---|---|---:|---|---:|---:|---:|---|---|---|---|")
            for entry in sorted(case_entries, key=lambda item: (item.backend_dir, item.config.config_id)):
                result = result_by_key[(entry.backend, entry.case.case_id, entry.config.config_id)]
                mesh = result.get("mesh_summary", {})
                links = primary_artifact_links(run_root=run_root, entry=entry, result=result)
                first_step = result["steps"][0] if result["steps"] else None
                timing = first_step.get("timing", {}) if first_step else {}
                wall_ms = first_step.get("wall_ms") if first_step else None
                total_ms = timing.get("total_ms") if timing else None

                elements = "-"
                if mesh.get("tris") is not None:
                    elements = f"tris={mesh['tris']}"
                elif mesh.get("tets") is not None:
                    elements = f"tets={mesh['tets']}"
                if mesh.get("edges") is not None:
                    elements = f"{elements}, edges={mesh['edges']}" if elements != "-" else f"edges={mesh['edges']}"
                if mesh.get("nnz_lower_s") is not None:
                    nnz_bits = f"nnzS={mesh['nnz_lower_s']}, nnzT={mesh.get('nnz_lower_t', '-')}"
                    elements = f"{elements}, {nnz_bits}" if elements != "-" else nnz_bits

                lines.append(
                    f"| `{entry.backend}` | `{entry.config.config_id}` | `{', '.join(entry.config.purposes)}` | "
                    f"`{mesh.get('mesh', '-')}` | {format_opt_int(mesh.get('nodes'))} | `{elements}` | "
                    f"{format_opt_int(mesh.get('dof'))} | {format_opt_float(wall_ms)} | {format_opt_float(total_ms)} | "
                    f"{links['csv']} | {links['img']} | {links['logs']} | {links['meta']} |"
                )
            lines.append("")
    lines.append("")
    return "\n".join(lines) + "\n"


def ensure_execution_root(
    *,
    out_base: Path,
    requested_run_name: str | None,
    now: dt.datetime,
) -> tuple[str, Path]:
    out_base = out_base.resolve()
    if requested_run_name:
        run_name = sanitize_run_name(requested_run_name)
        run_root = out_base / run_name
        if run_root.exists():
            raise SystemExit(f"Execution folder already exists: {run_root}")
        return run_name, run_root

    run_name = unique_auto_run_name(out_base, now)
    return run_name, out_base / run_name


def write_json(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def artifact_paths(config_root: Path) -> dict[str, Path]:
    meta_dir = config_root / "meta"
    return {
        "root": config_root,
        "vtk": config_root / "vtk",
        "img": config_root / "img",
        "csv": config_root / "csv",
        "logs": config_root / "logs",
        "meta": meta_dir,
        "tp3485_out": meta_dir / "tp3485_out",
    }


def create_structure(run_root: Path, entries: list[ConfigEntry]) -> None:
    run_root.mkdir(parents=False, exist_ok=False)

    case_groups: dict[tuple[str, str, str], list[ConfigEntry]] = defaultdict(list)
    for entry in entries:
        backend_root = run_root / entry.backend_dir
        section_root = backend_root / entry.case.section
        case_root = section_root / entry.case.case_id
        config_root = case_root / entry.config.config_id
        paths = artifact_paths(config_root)

        backend_root.mkdir(parents=True, exist_ok=True)
        section_root.mkdir(parents=True, exist_ok=True)
        case_root.mkdir(parents=True, exist_ok=True)
        config_root.mkdir(parents=True, exist_ok=True)
        for leaf in ARTIFACT_SUBDIRS:
            (config_root / leaf).mkdir(parents=True, exist_ok=True)
        paths["tp3485_out"].mkdir(parents=True, exist_ok=True)

        case_groups[(entry.backend_dir, entry.case.section, entry.case.case_id)].append(entry)

        config_meta = {
            "backend": entry.backend,
            "backend_dir": entry.backend_dir,
            "section": entry.case.section,
            "case_id": entry.case.case_id,
            "program": entry.case.program,
            "case_label": entry.case.label,
            "config_id": entry.config.config_id,
            "config_label": entry.config.label,
            "purposes": list(entry.config.purposes),
            "params": entry.config.params,
            "root": str(config_root),
            "root_from_repo": rel_to_root(config_root),
            "artifact_dirs": {
                leaf: rel_to_root(config_root / leaf) for leaf in ARTIFACT_SUBDIRS
            },
            "tp3485_out_dir": rel_to_root(paths["tp3485_out"]),
        }
        write_json(config_root / "config.json", config_meta)

    for (backend_dir, section, case_id), grouped_entries in case_groups.items():
        sample = grouped_entries[0]
        case_root = run_root / backend_dir / section / case_id
        case_meta = {
            "backend": sample.backend,
            "backend_dir": backend_dir,
            "section": section,
            "case_id": case_id,
            "program": sample.case.program,
            "case_label": sample.case.label,
            "family": sample.case.family,
            "configurations": [
                {
                    "config_id": entry.config.config_id,
                    "config_label": entry.config.label,
                    "purposes": list(entry.config.purposes),
                    "params": entry.config.params,
                    "path_from_repo": rel_to_root(run_root / entry.rel_dir),
                }
                for entry in grouped_entries
            ],
        }
        write_json(case_root / "CASE.json", case_meta)


def build_project(build_dir: Path) -> None:
    jobs = os.cpu_count() or 4
    subprocess.run(
        ["cmake", "-S", str(ROOT), "-B", str(build_dir), "-DCMAKE_BUILD_TYPE=Release"],
        check=True,
    )
    subprocess.run(["cmake", "--build", str(build_dir), "-j", str(jobs)], check=True)


def step_file_pair(logs_dir: Path, step_name: str) -> tuple[Path, Path]:
    stdout_path = logs_dir / f"{step_name}.stdout.txt"
    stderr_path = logs_dir / f"{step_name}.stderr.txt"
    return stdout_path, stderr_path


def run_logged_command(
    *,
    cmd: list[str],
    cwd: Path,
    env: dict[str, str],
    logs_dir: Path,
    step_name: str,
) -> CommandResult:
    stdout_path, stderr_path = step_file_pair(logs_dir, step_name)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    wall_ms = (time.perf_counter() - t0) * 1000.0
    stdout_path.write_text(proc.stdout, encoding="utf-8")
    stderr_path.write_text(proc.stderr, encoding="utf-8")
    timing = parse_timing(proc.stdout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Falha ao executar: {' '.join(cmd)}\n"
            f"stdout: {stdout_path}\n"
            f"stderr: {stderr_path}"
        )
    return CommandResult(
        step_name=step_name,
        command=cmd,
        wall_ms=wall_ms,
        returncode=proc.returncode,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        stdout_text=proc.stdout,
        stderr_text=proc.stderr,
        timing=timing,
    )


def run_logged_callable(
    *,
    logs_dir: Path,
    step_name: str,
    fn,
) -> dict[str, object]:
    stdout_path, stderr_path = step_file_pair(logs_dir, step_name)
    returncode = 0
    t0 = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8") as stdout_f, stderr_path.open("w", encoding="utf-8") as stderr_f:
        with contextlib.redirect_stdout(stdout_f), contextlib.redirect_stderr(stderr_f):
            try:
                fn()
            except Exception:
                traceback.print_exc(file=stderr_f)
                returncode = 1
    wall_ms = (time.perf_counter() - t0) * 1000.0
    if returncode != 0:
        raise RuntimeError(
            f"Falha ao executar passo Python: {step_name}\n"
            f"stdout: {stdout_path}\n"
            f"stderr: {stderr_path}"
        )
    return {
        "step_name": step_name,
        "command": [],
        "wall_ms": wall_ms,
        "returncode": 0,
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
    }


def write_rows_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def list_files_with_suffix(base: Path, suffix: str) -> list[str]:
    if not base.exists():
        return []
    return [rel_to_root(path) for path in sorted(base.rglob(f"*{suffix}"))]


def case_geometry_spec(case_id: str) -> tuple[str, str, str]:
    mapping = {
        "helm10_rect": ("rect", "scalar", "n"),
        "helm10_circle": ("circle", "scalar", "p"),
        "helm10_coax": ("coax", "scalar", "p"),
        "edge_rect": ("rect", "edge", "n"),
        "edge_circle": ("circle", "edge", "p"),
        "edge_coax": ("coax", "edge", "p"),
    }
    return mapping[case_id]


def format_scalar(value: object) -> str:
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:g}"
    return str(value)


def extract_int(pattern: str, text: str) -> int | None:
    match = re.search(pattern, text)
    return int(match.group(1)) if match else None


def mesh_label_from_entry(entry: ConfigEntry) -> str:
    params = entry.config.params
    parts: list[str] = []
    if "nr" in params and "nt" in params:
        parts.append(f"nr={format_scalar(params['nr'])}")
        parts.append(f"nt={format_scalar(params['nt'])}")
    if "nx" in params and "ny" in params:
        parts.append(f"nx={format_scalar(params['nx'])}")
        parts.append(f"ny={format_scalar(params['ny'])}")
    if "nz" in params:
        parts.append(f"nz={format_scalar(params['nz'])}")
    if "beta" in params:
        parts.append(f"beta={format_scalar(params['beta'])}")
    if "d_over_a" in params:
        parts.append(f"d/a={format_scalar(params['d_over_a'])}")
    return ", ".join(parts) if parts else "-"


def parse_mesh_summary(entry: ConfigEntry, stdout_text: str) -> dict[str, object]:
    params = entry.config.params
    nodes = extract_int(r"\bnodes=(\d+)", stdout_text)
    tris = extract_int(r"\btris=(\d+)", stdout_text)
    tets = extract_int(r"\btets=(\d+)", stdout_text)
    edges = extract_int(r"\bedges=(\d+)", stdout_text)
    edge_dofs = extract_int(r"\bedge_dofs=(\d+)", stdout_text)
    dof = extract_int(r"\bdof=(\d+)", stdout_text)
    nnz_lower_s = extract_int(r"\bnnz_lower\(S\)=(\d+)", stdout_text)
    nnz_lower_t = extract_int(r"\bnnz_lower\(T\)=(\d+)", stdout_text)

    if nodes is None and "nx" in params and "ny" in params and "nz" not in params:
        nodes = (int(params["nx"]) + 1) * (int(params["ny"]) + 1)
    if tris is None and "nx" in params and "ny" in params and "nz" not in params:
        tris = 2 * int(params["nx"]) * int(params["ny"])
    if dof is None and edge_dofs is not None:
        dof = edge_dofs

    return {
        "mesh": mesh_label_from_entry(entry),
        "nodes": nodes,
        "tris": tris,
        "tets": tets,
        "edges": edges,
        "dof": dof,
        "edge_dofs": edge_dofs,
        "nnz_lower_s": nnz_lower_s,
        "nnz_lower_t": nnz_lower_t,
    }


def repo_rel_to_abs(path_str: str) -> Path:
    return (ROOT / path_str).resolve()


def run_rel_link(run_root: Path, path: Path, label: str | None = None) -> str:
    rel = path.resolve().relative_to(run_root.resolve())
    return f"[{label or rel.name}]({rel.as_posix()})"


def format_opt_int(value: object | None) -> str:
    if value in {None, "", "None"}:
        return "-"
    return str(value)


def format_opt_float(value: object | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, float) and math.isnan(value):
        return "-"
    return f"{float(value):.3f}"


def primary_artifact_links(
    *,
    run_root: Path,
    entry: ConfigEntry,
    result: dict[str, object],
) -> dict[str, str]:
    csv_dir = run_root / entry.rel_dir / "csv"
    img_dir = run_root / entry.rel_dir / "img"
    logs_dir = run_root / entry.rel_dir / "logs"
    meta_dir = run_root / entry.rel_dir / "meta"

    csv_paths = [repo_rel_to_abs(p) for p in result["artifacts"]["csv"]]
    img_paths = [repo_rel_to_abs(p) for p in result["artifacts"]["img"]]
    steps = result["steps"]

    csv_bits = [run_rel_link(run_root, csv_dir, "csv/")]
    csv_bits.extend(run_rel_link(run_root, path) for path in csv_paths)

    img_bits = [run_rel_link(run_root, img_dir, f"img/ ({len(img_paths)})")]
    if img_paths:
        img_bits.append(run_rel_link(run_root, img_paths[0], "sample"))

    log_bits = [run_rel_link(run_root, logs_dir, "logs/")]
    if steps:
        first_step = steps[0]
        log_bits.append(run_rel_link(run_root, repo_rel_to_abs(first_step["stdout"]), "stdout"))
        log_bits.append(run_rel_link(run_root, repo_rel_to_abs(first_step["stderr"]), "stderr"))

    meta_bits = [
        run_rel_link(run_root, meta_dir, "meta/"),
        run_rel_link(run_root, meta_dir / "result.json", "result.json"),
    ]

    return {
        "csv": " ".join(csv_bits),
        "img": " ".join(img_bits),
        "logs": " ".join(log_bits),
        "meta": " ".join(meta_bits),
    }


def build_command(entry: ConfigEntry, build_dir: Path) -> list[str]:
    exe = build_dir / entry.case.program
    if not exe.exists():
        raise FileNotFoundError(f"Executavel nao encontrado: {exe}")

    params = entry.config.params
    case_id = entry.case.case_id

    if case_id in {"helm10_rect", "edge_rect"}:
        return [
            str(exe),
            str(params["nx"]),
            str(params["ny"]),
            str(params["mode_export"]),
            "--backend",
            entry.backend,
        ]
    if case_id in {"helm10_circle", "helm10_coax", "edge_circle", "edge_coax"}:
        return [
            str(exe),
            str(params["nr"]),
            str(params["nt"]),
            str(params["mode_export"]),
            "--backend",
            entry.backend,
        ]
    if case_id == "mixed_rect":
        return [str(exe), str(params["nx"]), str(params["ny"]), "--backend", entry.backend]
    if case_id in {"mixed_circle", "mixed_coax"}:
        return [str(exe), str(params["nr"]), str(params["nt"]), "--backend", entry.backend]
    if case_id == "helmvec2_rect":
        return [
            str(exe),
            str(params["beta"]),
            str(params["nx"]),
            str(params["ny"]),
            str(params["legacy_debug"]),
            "--backend",
            entry.backend,
        ]
    if case_id == "helmvec3_fig12_rect":
        return [
            str(exe),
            str(params["nx"]),
            str(params["ny"]),
            str(params["legacy_debug"]),
            "--backend",
            entry.backend,
        ]
    if case_id == "helmvec3_fig13_rect":
        return [
            str(exe),
            str(params["d_over_a"]),
            str(params["nx"]),
            str(params["ny"]),
            str(params["legacy_debug"]),
            "--backend",
            entry.backend,
        ]
    if case_id.startswith("fem3d0_") or case_id.startswith("fem3d1_"):
        return [
            str(exe),
            "--nx",
            str(params["nx"]),
            "--ny",
            str(params["ny"]),
            "--nz",
            str(params["nz"]),
            "--backend",
            entry.backend,
        ]
    raise SystemExit(f"Missing command builder for case: {case_id}")


def normalize_vtk_outputs(staging_root: Path, vtk_dir: Path) -> list[Path]:
    vtk_dir.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    for src in sorted(staging_root.rglob("*.vtk")):
        dest = vtk_dir / src.name
        if dest.exists():
            prefix = "__".join(src.relative_to(staging_root).parts[:-1])
            dest = vtk_dir / f"{prefix}__{src.name}" if prefix else vtk_dir / src.name
        shutil.copy2(src, dest)
        copied.append(dest)
    return copied


def build_validation_2d_rows(entry: ConfigEntry, stdout_text: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    case_id = entry.case.case_id

    if case_id == "mixed_rect":
        rows.extend(
            parse_mixed_rect_table(
                stdout_text,
                "[E-formulation] TE cutoffs (edge block)",
                "mixed_rect_E_TE_table",
            )
        )
        rows.extend(
            parse_mixed_rect_table(
                stdout_text,
                "[E-formulation] TM cutoffs (scalar block)",
                "mixed_rect_E_TM_table",
            )
        )
    elif case_id == "mixed_circle":
        snapshots = [
            ("2.2.2", "mixed_circle_TE_edge", parse_first_kc_block(stdout_text, "TE (edge block)")),
            ("2.2.2", "mixed_circle_TM_scalar", parse_first_kc_block(stdout_text, "TM (scalar block)")),
        ]
        for section, case_name, values in snapshots:
            for idx, value in enumerate(values, start=1):
                rows.append(
                    {
                        "section": section,
                        "case": case_name,
                        "mode": idx,
                        "fem": value,
                        "ref_primary": "",
                        "ref_secondary": "",
                        "err_primary_pct": "",
                        "err_secondary_pct": "",
                    }
                )
    elif case_id == "mixed_coax":
        snapshots = [
            ("2.2.2", "mixed_coax_TE_edge", parse_first_kc_block(stdout_text, "TE (edge block)")),
            ("2.2.2", "mixed_coax_TM_scalar", parse_first_kc_block(stdout_text, "TM (scalar block)")),
        ]
        for section, case_name, values in snapshots:
            for idx, value in enumerate(values, start=1):
                rows.append(
                    {
                        "section": section,
                        "case": case_name,
                        "mode": idx,
                        "fem": value,
                        "ref_primary": "",
                        "ref_secondary": "",
                        "err_primary_pct": "",
                        "err_secondary_pct": "",
                    }
                )
    elif case_id == "helmvec2_rect":
        rows.extend(parse_helmvec2_table(stdout_text))
    elif case_id == "helmvec3_fig12_rect":
        rows.extend(parse_helmvec3_table9(stdout_text))
    elif case_id == "helmvec3_fig13_rect":
        rows.extend(parse_helmvec3_table10(stdout_text))

    for row in rows:
        row["backend"] = entry.backend
    return rows


def write_validation_2d_csv(csv_path: Path, rows: list[dict[str, object]]) -> None:
    write_rows_csv(
        csv_path,
        rows,
        [
            "backend",
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


def build_validation_3d_rows(entry: ConfigEntry, stdout_text: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    parsed_rows = parse_mode_table_3d(stdout_text)
    solver = entry.case.case_id.split("_", 1)[0]
    case_name = entry.case.case_id.split("_", 1)[1]
    nx = int(entry.config.params["nx"])
    ny = int(entry.config.params["ny"])
    nz = int(entry.config.params["nz"])

    mode_rows: list[dict[str, object]] = []
    for row in parsed_rows:
        mode_rows.append(
            {
                "section": "3.1",
                "solver": solver,
                "backend": entry.backend,
                "case": case_name,
                "nx": nx,
                "ny": ny,
                "nz": nz,
                **row,
            }
        )

    if not mode_rows:
        return [], []

    err_ana = [abs(float(r["err_ana_pct"])) for r in mode_rows]
    err_ref = [abs(float(r["err_ref_pct"])) for r in mode_rows]
    summary_rows = [
        {
            "section": "3.1",
            "solver": solver,
            "backend": entry.backend,
            "case": case_name,
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "n_modes": len(mode_rows),
            "max_err_ana_pct": max(err_ana),
            "mean_err_ana_pct": statistics.fmean(err_ana),
            "max_err_ref_pct": max(err_ref),
            "mean_err_ref_pct": statistics.fmean(err_ref),
        }
    ]
    return mode_rows, summary_rows


def write_validation_3d_csvs(
    *,
    modes_csv: Path,
    summary_csv: Path,
    mode_rows: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
) -> None:
    write_rows_csv(
        modes_csv,
        mode_rows,
        [
            "section",
            "solver",
            "backend",
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
    write_rows_csv(
        summary_csv,
        summary_rows,
        [
            "section",
            "solver",
            "backend",
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


def execute_entry(entry: ConfigEntry, *, build_dir: Path, args: argparse.Namespace, run_root: Path) -> dict[str, object]:
    config_root = run_root / entry.rel_dir
    paths = artifact_paths(config_root)
    env = os.environ.copy()
    env["TP3485_OUT_DIR"] = str(paths["tp3485_out"])
    env.setdefault("MPLBACKEND", "Agg")

    step_index = 1
    logs_dir = paths["logs"]
    direct_step = f"{step_index:02d}_run_{entry.case.program}"
    command = build_command(entry, build_dir)
    direct_result = run_logged_command(
        cmd=command,
        cwd=build_dir,
        env=env,
        logs_dir=logs_dir,
        step_name=direct_step,
    )
    mesh_summary = parse_mesh_summary(entry, direct_result.stdout_text)

    steps: list[dict[str, object]] = [
        {
            "step_name": direct_result.step_name,
            "kind": "command",
            "command": direct_result.command,
            "wall_ms": round(direct_result.wall_ms, 3),
            "returncode": direct_result.returncode,
            "stdout": rel_to_root(direct_result.stdout_path),
            "stderr": rel_to_root(direct_result.stderr_path),
            "timing": {
                key: (None if math.isnan(value) else round(value, 3))
                for key, value in direct_result.timing.items()
            },
        }
    ]

    artifact_rel: dict[str, list[str]] = {"vtk": [], "img": [], "csv": [], "logs": [], "meta": []}
    artifact_rel["logs"] = list_files_with_suffix(paths["logs"], ".txt")

    copied_vtk = normalize_vtk_outputs(paths["tp3485_out"], paths["vtk"])
    artifact_rel["vtk"] = [rel_to_root(path) for path in copied_vtk]

    if entry.case.case_id in FIELD_CASE_IDS:
        geometry, formulation, index2_name = case_geometry_spec(entry.case.case_id)
        mode_rows = parse_mode_table_2d(
            direct_result.stdout_text,
            geometry,
            formulation,
            index2_name,
            parse_dims_2d(direct_result.stdout_text),
        )
        mode_csv = paths["csv"] / "mode_summary.csv"

        def post_fields() -> None:
            write_mode_summary_csv(mode_rows, mode_csv)
            if not args.skip_images and copied_vtk:
                plot_all_field_images(
                    vtk_root=paths["vtk"],
                    out_dir=paths["img"],
                    rows=mode_rows,
                    stride=2,
                    scale=22.0,
                    show_mesh=True,
                    dpi=210,
                    max_rank=int(entry.config.params.get("mode_export", args.mode_export)),
                )

        step_index += 1
        post_result = run_logged_callable(
            logs_dir=logs_dir,
            step_name=f"{step_index:02d}_fields_post",
            fn=post_fields,
        )
        steps.append(
            {
                "step_name": post_result["step_name"],
                "kind": "python",
                "command": [],
                "wall_ms": round(float(post_result["wall_ms"]), 3),
                "returncode": 0,
                "stdout": rel_to_root(Path(post_result["stdout_path"])),
                "stderr": rel_to_root(Path(post_result["stderr_path"])),
            }
        )
    elif entry.case.case_id in VALIDATION_2D_CASE_IDS and not args.skip_validate:
        validation_csv = paths["csv"] / "validation_2d_22.csv"
        validation_rows = build_validation_2d_rows(entry, direct_result.stdout_text)
        write_validation_2d_csv(validation_csv, validation_rows)

        if validation_rows and not args.skip_images:
            step_index += 1
            plot_cmd = [
                "python3",
                str(ROOT / "scripts" / "plot_validation_2d_22.py"),
                "--in-csv",
                str(validation_csv),
                "--out-dir",
                str(paths["img"]),
                "--backend",
                entry.backend,
            ]
            plot_result = run_logged_command(
                cmd=plot_cmd,
                cwd=ROOT,
                env=env,
                logs_dir=logs_dir,
                step_name=f"{step_index:02d}_plot_validation_2d_22",
            )
            steps.append(
                {
                    "step_name": plot_result.step_name,
                    "kind": "command",
                    "command": plot_result.command,
                    "wall_ms": round(plot_result.wall_ms, 3),
                    "returncode": plot_result.returncode,
                    "stdout": rel_to_root(plot_result.stdout_path),
                    "stderr": rel_to_root(plot_result.stderr_path),
                    "timing": {
                        key: (None if math.isnan(value) else round(value, 3))
                        for key, value in plot_result.timing.items()
                    },
                }
            )
    elif entry.case.case_id.startswith("fem3d") and not args.skip_validate:
        mode_rows, summary_rows = build_validation_3d_rows(entry, direct_result.stdout_text)
        if mode_rows:
            write_validation_3d_csvs(
                modes_csv=paths["csv"] / "validation_3d_31_modes.csv",
                summary_csv=paths["csv"] / "validation_3d_31_summary.csv",
                mode_rows=mode_rows,
                summary_rows=summary_rows,
            )

    artifact_rel["img"] = list_files_with_suffix(paths["img"], ".png")
    artifact_rel["csv"] = list_files_with_suffix(paths["csv"], ".csv")
    artifact_rel["logs"] = list_files_with_suffix(paths["logs"], ".txt")
    artifact_rel["meta"] = [
        rel_to_root(path)
        for path in sorted(paths["meta"].rglob("*"))
        if path.is_file() and path.name != "result.json"
    ]

    result_payload = {
        "backend": entry.backend,
        "section": entry.case.section,
        "case_id": entry.case.case_id,
        "program": entry.case.program,
        "config_id": entry.config.config_id,
        "config_label": entry.config.label,
        "purposes": list(entry.config.purposes),
        "params": entry.config.params,
        "tp3485_out_dir": rel_to_root(paths["tp3485_out"]),
        "status": "completed",
        "mesh_summary": mesh_summary,
        "steps": steps,
        "artifacts": artifact_rel,
    }
    write_json(paths["meta"] / "result.json", result_payload)

    config_json_path = config_root / "config.json"
    config_payload = read_json(config_json_path)
    config_payload["execution_result"] = {
        "status": "completed",
        "result_json": rel_to_root(paths["meta"] / "result.json"),
        "mesh_summary": mesh_summary,
        "generated_counts": {
            "vtk": len(artifact_rel["vtk"]),
            "img": len(artifact_rel["img"]),
            "csv": len(artifact_rel["csv"]),
            "logs": len(artifact_rel["logs"]),
        },
    }
    write_json(config_json_path, config_payload)
    return result_payload


def execution_summary_from_results(results: list[dict[str, object]]) -> dict[str, object]:
    return {
        "status": "completed",
        "executed_configurations": len(results),
        "vtk_files": sum(len(result["artifacts"]["vtk"]) for result in results),
        "png_files": sum(len(result["artifacts"]["img"]) for result in results),
        "csv_files": sum(len(result["artifacts"]["csv"]) for result in results),
        "log_files": sum(len(result["artifacts"]["logs"]) for result in results),
    }


def main() -> int:
    args = parse_args()
    now = dt.datetime.now()

    if args.mode_export < 1:
        raise SystemExit("Invalid --mode-export: expected integer >= 1.")
    if args.node_scale <= 0:
        raise SystemExit("Invalid --node-scale: expected a positive number.")

    out_base = args.out_base.resolve()
    build_dir = args.build_dir.resolve()
    backends = resolve_backends(args.backend)
    selected_cases = resolve_cases(args.case)
    entries = build_entries(
        backends=backends,
        selected_cases=selected_cases,
        mode_export=args.mode_export,
        profile=args.profile,
        article_only=args.article_only,
        node_scale=args.node_scale,
    )
    run_name, run_root = ensure_execution_root(
        out_base=out_base,
        requested_run_name=args.run_name,
        now=now,
    )

    print(f"[structured-run] run_name={run_name}")
    print(f"[structured-run] run_root={run_root}")
    print(f"[structured-run] backends={','.join(backends)}")
    print(f"[structured-run] cases={len(selected_cases)}")
    print(f"[structured-run] config_dirs={len(entries)}")

    if args.dry_run:
        for entry in entries[:8]:
            print(f"[structured-run] plan={entry.rel_dir.as_posix()}")
        if len(entries) > 8:
            print(f"[structured-run] ... {len(entries) - 8} more configuration roots omitted")
        print("[structured-run] dry-run active: no files or directories were created.")
        return 0

    out_base.mkdir(parents=True, exist_ok=True)
    create_structure(run_root, entries)

    initial_manifest = build_manifest(
        args=args,
        now=now,
        run_name=run_name,
        run_root=run_root,
        backends=backends,
        selected_cases=selected_cases,
        entries=entries,
    )
    initial_run_md = render_run_md(
        now=now,
        run_name=run_name,
        run_root=run_root,
        build_dir=build_dir,
        backends=backends,
        selected_cases=selected_cases,
        entries=entries,
        args=args,
    )
    write_json(run_root / "manifest.json", initial_manifest)
    (run_root / "RUN.md").write_text(initial_run_md, encoding="utf-8")

    if not args.skip_build:
        print("[structured-run] configuring and building project...")
        build_project(build_dir)
    elif not build_dir.exists():
        raise SystemExit(f"Build directory not found: {build_dir}")

    results: list[dict[str, object]] = []
    for entry in entries:
        print(f"[structured-run] executing={entry.rel_dir.as_posix()}")
        results.append(execute_entry(entry, build_dir=build_dir, args=args, run_root=run_root))

    summary = execution_summary_from_results(results)
    final_manifest = build_manifest(
        args=args,
        now=now,
        run_name=run_name,
        run_root=run_root,
        backends=backends,
        selected_cases=selected_cases,
        entries=entries,
        execution_summary=summary,
    )
    final_run_md = render_run_md(
        now=now,
        run_name=run_name,
        run_root=run_root,
        build_dir=build_dir,
        backends=backends,
        selected_cases=selected_cases,
        entries=entries,
        args=args,
        execution_summary=summary,
        results=results,
    )
    write_json(run_root / "manifest.json", final_manifest)
    (run_root / "RUN.md").write_text(final_run_md, encoding="utf-8")

    print(f"[structured-run] manifest={run_root / 'manifest.json'}")
    print(f"[structured-run] report={run_root / 'RUN.md'}")
    print(
        "[structured-run] generated="
        f"vtk:{summary['vtk_files']} png:{summary['png_files']} csv:{summary['csv_files']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
