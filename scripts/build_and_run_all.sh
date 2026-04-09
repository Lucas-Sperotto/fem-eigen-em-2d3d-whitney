#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"
OUT_DIR="$ROOT_DIR/out"
BUILD_TYPE="Release"
PROFILE="quick"
MODE_EXPORT=8
BACKEND="closed-form"
SKIP_VALIDATE=0
SKIP_IMAGES=0
SKIP_3D=0
VERBOSE=0
NO_LOG=0
FORCE_VALIDATE=0
FORCE_IMAGES=0
RUN_BACKEND_BENCH=0
BACKEND_BENCH_SUITE="core"
BACKEND_BENCH_REPEATS=3
DEBUG_LOCAL_BLOCKS=0
DEBUG_CANDIDATES=0
SHOW_VALIDATE_OUTPUT=0
JOBS="$(nproc 2>/dev/null || echo 4)"
LOG_FILE=""
declare -a CASES=()

# Run flags (default: full pipeline).
RUN_21=1
RUN_221=1
RUN_222=1
RUN_223=1
RUN_224=1
RUN_3D_AIR=1
RUN_3D_HALF=1
RUN_3D_CYL=1
RUN_3D_SPHERE=1

usage() {
  cat <<'EOF'
Usage: scripts/build_and_run_all.sh [options]

Compile and run the full project pipeline:
- build CMake targets
- run 2D executables (helm10, helmvec, helmvec1, helmvec2, helmvec3)
- run 3D executables (fem3d0, fem3d1)
- run validation scripts (2D + 3D)
- generate post-processing images and mode summary CSV

Options:
  --build-dir <dir>      Build directory (default: build)
  --out-dir <dir>        Root output directory (default: out)
  --build-type <type>    CMake build type (default: Release)
  --jobs <N>             Parallel build jobs (default: nproc)
  --profile <quick|full> Profile for validate_3d_31.py (default: quick)
  --mode-export <N>      Number of exported 2D modes per TE/TM block (default: 8)
  --backend <name>       Backend: closed-form|gauss (default: closed-form)
  --case <id>            Run only selected section/case (repeatable)
  --log-file <path>      Write console log to file (default: <out-dir>/run_all.log)
  --no-log               Disable log file output
  --skip-3d              Skip fem3d0/fem3d1 runs and 3D validation
  --skip-validate        Skip validation scripts (2D CSV aggregators, validate_2d_22.py and validate_3d_31.py)
  --with-validate        Force validations even in --case mode
  --skip-images          Skip post-processing image generation (plot_vtk_quiver.py / plot_validation_2d_22.py)
  --with-images          Force post-processing image generation even in --case mode
  --benchmark-backends   Run gauss vs closed-form benchmark after build
  --benchmark-suite <s>  Benchmark suite: core|all (default: core)
  --benchmark-repeats N  Benchmark repeats per backend (default: 3)
  --debug-local-blocks   Propagate local closed-form block debug to supported executables
  --debug-candidates     Propagate candidate/root debug to supported executables
  --show-validation-output  Echo raw stdout captured by validation scripts
  --verbose              Print executed commands
  --help                 Show this help

Case aliases examples:
  --case 2.1             Scalar 2D block (helm10: tables 1-3)
  --case 2.2.1           Edge-only transverse vector block
  --case 2.2.2           Coupled vector+scalar cutoff block
  --case 2.2.3           Figure 11 / Table 8 (helmvec2)
  --case 2.2.4           Figure 12-13 / Tables 9-10 (helmvec3)
  --case 3.1             All 3D tables (12-15) in fem3d0/fem3d1
  --case table13         Only Table 13 / Figure 16 (half-filled cavity)
EOF
}

log() {
  printf '[run_all] %s\n' "$*"
}

run() {
  if [[ "$VERBOSE" -eq 1 ]]; then
    log "CMD: $*"
  fi
  "$@"
}

disable_all_cases() {
  RUN_21=0
  RUN_221=0
  RUN_222=0
  RUN_223=0
  RUN_224=0
  RUN_3D_AIR=0
  RUN_3D_HALF=0
  RUN_3D_CYL=0
  RUN_3D_SPHERE=0
}

enable_2d_all() {
  RUN_21=1
  RUN_221=1
  RUN_222=1
  RUN_223=1
  RUN_224=1
}

enable_3d_all() {
  RUN_3D_AIR=1
  RUN_3D_HALF=1
  RUN_3D_CYL=1
  RUN_3D_SPHERE=1
}

normalize_case_id() {
  local v="$1"
  v="${v,,}"
  v="${v//./}"
  v="${v//_/}"
  v="${v//-/}"
  v="${v//\//}"
  v="${v//:/}"
  echo "$v"
}

select_case_id() {
  local raw="$1"
  local c
  c="$(normalize_case_id "$raw")"
  case "$c" in
    all)
      enable_2d_all
      enable_3d_all
      ;;

    2|2d|sec2|section2)
      enable_2d_all
      ;;
    21|sec21|section21|scalar|helm10|table1|table2|table3|tabela1|tabela2|tabela3)
      RUN_21=1
      ;;
    221|sec221|section221|edge)
      RUN_221=1
      ;;
    222|sec222|section222|mixed|coupled|acoplado)
      RUN_222=1
      ;;
    223|sec223|section223|fig11|figure11|table8|tabela8|helmvec2)
      RUN_223=1
      ;;
    224|sec224|section224|fig12|figure12|fig13|figure13|table9|table10|tabela9|tabela10|helmvec3)
      RUN_224=1
      ;;

    31|sec31|section31|3d|fem3d|all3d)
      enable_3d_all
      ;;
    table12|tabela12|fig15|figure15|air)
      RUN_3D_AIR=1
      ;;
    table13|tabela13|fig16|figure16|half|halffilled)
      RUN_3D_HALF=1
      ;;
    table14|tabela14|fig17|figure17|cyl|cylinder|cil|cilindro)
      RUN_3D_CYL=1
      ;;
    table15|tabela15|sphere|spherical|esfera)
      RUN_3D_SPHERE=1
      ;;
    *)
      echo "Unknown --case value: $raw" >&2
      echo "Try: 2.1, 2.2.3, 2.2.4, 3.1, table12, table13, table14, table15" >&2
      exit 2
      ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      BUILD_DIR="$2"
      shift 2
      ;;
    --build-type)
      BUILD_TYPE="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --mode-export)
      MODE_EXPORT="$2"
      shift 2
      ;;
    --backend)
      BACKEND="$2"
      shift 2
      ;;
    --case)
      CASES+=("$2")
      shift 2
      ;;
    --log-file)
      LOG_FILE="$2"
      shift 2
      ;;
    --no-log)
      NO_LOG=1
      shift
      ;;
    --skip-3d)
      SKIP_3D=1
      shift
      ;;
    --skip-validate)
      SKIP_VALIDATE=1
      FORCE_VALIDATE=0
      shift
      ;;
    --with-validate)
      SKIP_VALIDATE=0
      FORCE_VALIDATE=1
      shift
      ;;
    --skip-images)
      SKIP_IMAGES=1
      FORCE_IMAGES=0
      shift
      ;;
    --with-images)
      SKIP_IMAGES=0
      FORCE_IMAGES=1
      shift
      ;;
    --benchmark-backends)
      RUN_BACKEND_BENCH=1
      shift
      ;;
    --benchmark-suite)
      BACKEND_BENCH_SUITE="$2"
      shift 2
      ;;
    --benchmark-repeats)
      BACKEND_BENCH_REPEATS="$2"
      shift 2
      ;;
    --debug-local-blocks)
      DEBUG_LOCAL_BLOCKS=1
      shift
      ;;
    --debug-candidates)
      DEBUG_CANDIDATES=1
      shift
      ;;
    --show-validation-output)
      SHOW_VALIDATE_OUTPUT=1
      shift
      ;;
    --verbose)
      VERBOSE=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ "$PROFILE" != "quick" && "$PROFILE" != "full" ]]; then
  echo "Invalid --profile: $PROFILE (expected quick|full)" >&2
  exit 2
fi

if [[ "$BACKEND" != "gauss" && "$BACKEND" != "closed-form" ]]; then
  echo "Invalid --backend: $BACKEND (expected gauss|closed-form)" >&2
  exit 2
fi

if [[ "$BACKEND_BENCH_SUITE" != "core" && "$BACKEND_BENCH_SUITE" != "all" ]]; then
  echo "Invalid --benchmark-suite: $BACKEND_BENCH_SUITE (expected core|all)" >&2
  exit 2
fi

if ! [[ "$BACKEND_BENCH_REPEATS" =~ ^[0-9]+$ ]] || [[ "$BACKEND_BENCH_REPEATS" -lt 1 ]]; then
  echo "Invalid --benchmark-repeats: $BACKEND_BENCH_REPEATS (expected integer >= 1)" >&2
  exit 2
fi

case "$BUILD_DIR" in
  /*) ;;
  *) BUILD_DIR="$ROOT_DIR/$BUILD_DIR" ;;
esac
case "$OUT_DIR" in
  /*) ;;
  *) OUT_DIR="$ROOT_DIR/$OUT_DIR" ;;
esac

if ! [[ "$MODE_EXPORT" =~ ^[0-9]+$ ]] || [[ "$MODE_EXPORT" -lt 1 ]]; then
  echo "Invalid --mode-export: $MODE_EXPORT (expected integer >= 1)" >&2
  exit 2
fi

if [[ -n "$LOG_FILE" ]]; then
  case "$LOG_FILE" in
    /*) ;;
    *) LOG_FILE="$ROOT_DIR/$LOG_FILE" ;;
  esac
else
  LOG_FILE="$OUT_DIR/run_all.log"
fi

if [[ "$NO_LOG" -eq 0 ]]; then
  mkdir -p "$(dirname "$LOG_FILE")"
  exec > >(tee "$LOG_FILE") 2>&1
fi

mkdir -p "$OUT_DIR"
export TP3485_OUT_DIR="$OUT_DIR"

CASE_MODE=0
if [[ "${#CASES[@]}" -gt 0 ]]; then
  CASE_MODE=1
  disable_all_cases
  for c in "${CASES[@]}"; do
    select_case_id "$c"
  done
fi

if [[ "$CASE_MODE" -eq 1 ]]; then
  # In selective mode we avoid failing because unrelated outputs are missing.
  if [[ "$SKIP_VALIDATE" -eq 0 && "$FORCE_VALIDATE" -eq 0 ]]; then
    SKIP_VALIDATE=1
    log "Case mode: validations disabled by default. Use --with-validate to force."
  fi
  if [[ "$SKIP_IMAGES" -eq 0 && "$FORCE_IMAGES" -eq 0 ]]; then
    SKIP_IMAGES=1
    log "Case mode: image generation disabled by default. Use --with-images to force."
  fi
fi

log "ROOT_DIR=$ROOT_DIR"
log "BUILD_DIR=$BUILD_DIR"
log "OUT_DIR=$OUT_DIR"
log "BUILD_TYPE=$BUILD_TYPE"
log "JOBS=$JOBS"
log "PROFILE=$PROFILE"
log "MODE_EXPORT=$MODE_EXPORT"
log "BACKEND=$BACKEND"
if [[ "$NO_LOG" -eq 0 ]]; then
  log "LOG_FILE=$LOG_FILE"
fi
if [[ "$RUN_BACKEND_BENCH" -eq 1 ]]; then
  log "BACKEND_BENCH=1 suite=$BACKEND_BENCH_SUITE repeats=$BACKEND_BENCH_REPEATS"
fi
if [[ "$DEBUG_LOCAL_BLOCKS" -eq 1 ]]; then
  log "DEBUG_LOCAL_BLOCKS=1"
fi
if [[ "$DEBUG_CANDIDATES" -eq 1 ]]; then
  log "DEBUG_CANDIDATES=1"
fi
if [[ "$SHOW_VALIDATE_OUTPUT" -eq 1 ]]; then
  log "SHOW_VALIDATE_OUTPUT=1"
fi
if [[ "$CASE_MODE" -eq 1 ]]; then
  log "CASE_MODE=1 selected: ${CASES[*]}"
fi

declare -a DEBUG_ARGS=()
if [[ "$DEBUG_LOCAL_BLOCKS" -eq 1 ]]; then
  DEBUG_ARGS+=("--debug-local-blocks")
fi
if [[ "$DEBUG_CANDIDATES" -eq 1 ]]; then
  DEBUG_ARGS+=("--debug-candidates")
fi

log "Configuring CMake..."
run cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"

log "Building targets..."
run cmake --build "$BUILD_DIR" -j "$JOBS"

if [[ "$RUN_21" -eq 1 || "$RUN_221" -eq 1 || "$RUN_222" -eq 1 || "$RUN_223" -eq 1 || "$RUN_224" -eq 1 ]]; then
  log "Running 2D executables..."
  (
    cd "$BUILD_DIR"
    if [[ "$RUN_21" -eq 1 ]]; then
      HELM10_RECT_MODES="$MODE_EXPORT"
      HELM10_CIRCLE_MODES="$MODE_EXPORT"
      HELM10_COAX_MODES="$MODE_EXPORT"
      if [[ "$HELM10_RECT_MODES" -lt 10 ]]; then
        HELM10_RECT_MODES=10
      fi
      if [[ "$HELM10_CIRCLE_MODES" -lt 14 ]]; then
        HELM10_CIRCLE_MODES=14
      fi
      if [[ "$HELM10_COAX_MODES" -lt 8 ]]; then
        HELM10_COAX_MODES=8
      fi
      run ./helm10_rect --ar-m 1.0 --nx 14 --ny 14 --nmodos "$HELM10_RECT_MODES" --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./helm10_circle --nr 10 --nt 48 --nmodos "$HELM10_CIRCLE_MODES" --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./helm10_coax --nr 10 --nt 48 --nmodos "$HELM10_COAX_MODES" --backend "$BACKEND" "${DEBUG_ARGS[@]}"
    fi

    if [[ "$RUN_221" -eq 1 ]]; then
      run ./edge_rect --nx 14 --ny 14 --nmodos 20 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./edge_circle --nr 10 --nt 48 --nmodos 20 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./edge_coax 10 48 "$MODE_EXPORT" --backend "$BACKEND" "${DEBUG_ARGS[@]}"
    fi

    if [[ "$RUN_222" -eq 1 ]]; then
      run ./mixed_rect 12 6 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./mixed_circle 10 48 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./mixed_coax 10 48 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
    fi

    if [[ "$RUN_223" -eq 1 ]]; then
      run ./helmvec2_rect 10 6 6 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
    fi

    if [[ "$RUN_224" -eq 1 ]]; then
      run ./helmvec3_fig12_rect 10 5 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      run ./helmvec3_fig13_rect 0.20 10 5 --backend "$BACKEND" "${DEBUG_ARGS[@]}"
    fi
  )
else
  log "No 2D case selected."
fi

if [[ "$SKIP_3D" -eq 0 ]]; then
  if [[ "$RUN_3D_AIR" -eq 1 || "$RUN_3D_HALF" -eq 1 || "$RUN_3D_CYL" -eq 1 || "$RUN_3D_SPHERE" -eq 1 ]]; then
    log "Running 3D executables..."
    (
      cd "$BUILD_DIR"
      if [[ "$RUN_3D_AIR" -eq 1 ]]; then
        run ./fem3d0_air --backend "$BACKEND" "${DEBUG_ARGS[@]}"
        run ./fem3d1_air --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      fi
      if [[ "$RUN_3D_HALF" -eq 1 ]]; then
        run ./fem3d0_half --backend "$BACKEND" "${DEBUG_ARGS[@]}"
        run ./fem3d1_half --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      fi
      if [[ "$RUN_3D_CYL" -eq 1 ]]; then
        run ./fem3d0_cyl --backend "$BACKEND" "${DEBUG_ARGS[@]}"
        run ./fem3d1_cyl --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      fi
      if [[ "$RUN_3D_SPHERE" -eq 1 ]]; then
        run ./fem3d0_sphere --backend "$BACKEND" "${DEBUG_ARGS[@]}"
        run ./fem3d1_sphere --backend "$BACKEND" "${DEBUG_ARGS[@]}"
      fi
    )
  else
    log "No 3D case selected."
  fi
else
  log "Skipping 3D executables (--skip-3d)."
fi

if [[ "$SKIP_VALIDATE" -eq 0 ]]; then
  declare -a VALIDATE_OUTPUT_ARGS=()
  if [[ "$SHOW_VALIDATE_OUTPUT" -eq 1 ]]; then
    VALIDATE_OUTPUT_ARGS+=("--show-output")
  fi

  log "Running 2D validation..."
  run python3 "$ROOT_DIR/scripts/validate_2d_22.py" \
    --build-dir "$BUILD_DIR" \
    --backend "$BACKEND" \
    --out-csv "$OUT_DIR/validation/validation_2d_22.csv" \
    "${VALIDATE_OUTPUT_ARGS[@]}" \
    "${DEBUG_ARGS[@]}"

  if [[ "$RUN_21" -eq 1 ]]; then
    log "Running 2.1 CSV-based validation..."
    run python3 "$ROOT_DIR/scripts/validate_2d_21_csv.py" \
      --out-root "$OUT_DIR" \
      --backend "$BACKEND" \
      --out-csv "$OUT_DIR/validation/validation_2d_21.csv"
  else
    log "Skipping 2.1 CSV-based validation (no scalar 2D case selected)."
  fi

  if [[ "$RUN_221" -eq 1 ]]; then
    log "Running 2.2.1 CSV-based validation..."
    run python3 "$ROOT_DIR/scripts/validate_2d_221_csv.py" \
      --out-root "$OUT_DIR" \
      --backend "$BACKEND" \
      --out-csv "$OUT_DIR/validation/validation_2d_221.csv"
  else
    log "Skipping 2.2.1 CSV-based validation (no edge 2D case selected)."
  fi

  if [[ "$RUN_222" -eq 1 ]]; then
    log "Running 2.2.2 Table 6 CSV-based validation..."
    run python3 "$ROOT_DIR/scripts/validate_2d_222_table6_csv.py" \
      --out-root "$OUT_DIR" \
      --backend "$BACKEND" \
      --out-csv "$OUT_DIR/validation/validation_2d_222_table6.csv"

    log "Running 2.2.2 Table 7 CSV-based validation..."
    run python3 "$ROOT_DIR/scripts/validate_2d_222_table7_csv.py" \
      --out-root "$OUT_DIR" \
      --backend "$BACKEND" \
      --out-csv "$OUT_DIR/validation/validation_2d_222_table7.csv"
  else
    log "Skipping 2.2.2 Table 6 CSV-based validation (no mixed 2D case selected)."
    log "Skipping 2.2.2 Table 7 CSV-based validation (no mixed 2D case selected)."
  fi

  if [[ "$RUN_224" -eq 1 ]]; then
    log "Running 2.2.4 Figure 13 / Table 10 CSV-based validation..."
    run python3 "$ROOT_DIR/scripts/validate_2d_224_table10_csv.py" \
      --out-root "$OUT_DIR" \
      --backend "$BACKEND" \
      --out-csv "$OUT_DIR/validation/validation_2d_224_table10.csv"
  else
    log "Skipping 2.2.4 Figure 13 / Table 10 CSV-based validation (no HELMVEC3 Figure 13 case selected)."
  fi

  if [[ "$SKIP_3D" -eq 0 ]]; then
    log "Running 3D validation..."
    run python3 "$ROOT_DIR/scripts/validate_3d_31.py" \
      --profile "$PROFILE" \
      --solver both \
      --build-dir "$BUILD_DIR" \
      --backend "$BACKEND" \
      --out-modes "$OUT_DIR/validation/validation_3d_31_modes.csv" \
      --out-summary "$OUT_DIR/validation/validation_3d_31_summary.csv" \
      "${VALIDATE_OUTPUT_ARGS[@]}" \
      "${DEBUG_ARGS[@]}"
  else
    log "Skipping 3D validation because --skip-3d is active."
  fi

  log "Generating consolidated 2D validation index..."
  run python3 "$ROOT_DIR/scripts/generate_validation_2d_index.py" \
    --validation-dir "$OUT_DIR/validation" \
    --out-csv "$OUT_DIR/validation/validation_2d_index.csv" \
    --out-md "$OUT_DIR/validation/VALIDATION_2D_INDEX.md"

  log "Generating consolidated 3D validation index..."
  run python3 "$ROOT_DIR/scripts/generate_validation_3d_index.py" \
    --validation-dir "$OUT_DIR/validation" \
    --out-csv "$OUT_DIR/validation/validation_3d_index.csv" \
    --out-md "$OUT_DIR/validation/VALIDATION_3D_INDEX.md"

  log "Generating master scientific validation index..."
  run python3 "$ROOT_DIR/scripts/generate_validation_master_index.py" \
    --validation-dir "$OUT_DIR/validation" \
    --out-csv "$OUT_DIR/validation/validation_master_index.csv" \
    --out-md "$OUT_DIR/validation/VALIDATION_MASTER_INDEX.md"

  log "Generating scientific validation verdict..."
  run python3 "$ROOT_DIR/scripts/generate_validation_verdict.py" \
    --master-index "$OUT_DIR/validation/validation_master_index.csv" \
    --policy-csv "$ROOT_DIR/docs/validation_thresholds.csv" \
    --out-csv "$OUT_DIR/validation/validation_verdict.csv" \
    --out-md "$OUT_DIR/validation/VALIDATION_VERDICT.md"
else
  log "Skipping validations (--skip-validate)."
fi

if [[ "$SKIP_IMAGES" -eq 0 ]]; then
  if [[ "$RUN_21" -eq 1 || "$RUN_221" -eq 1 ]]; then
    log "Generating all images and mode summary CSV..."
    run python3 "$ROOT_DIR/scripts/plot_vtk_quiver.py" \
      --all-img \
      --build-dir "$BUILD_DIR" \
      --vtk-root "$OUT_DIR" \
      --out-dir "$OUT_DIR/img_all" \
      --csv "$OUT_DIR/img_all/mode_summary.csv" \
      --mode-export "$MODE_EXPORT" \
      --max-rank "$MODE_EXPORT"
  else
    log "Skipping 2.1/2.2.1 field images (no scalar/edge 2D case selected)."
  fi

  if [[ -f "$OUT_DIR/validation/validation_2d_22.csv" ]]; then
    log "Generating 2.2.2/2.2.3/2.2.4 validation figures..."
    run python3 "$ROOT_DIR/scripts/plot_validation_2d_22.py" \
      --in-csv "$OUT_DIR/validation/validation_2d_22.csv" \
      --out-dir "$OUT_DIR/img_all/validation_2d_22" \
      --backend "$BACKEND"
  else
    log "Skipping 2.2.x validation figures (missing $OUT_DIR/validation/validation_2d_22.csv)."
  fi
else
  log "Skipping image generation (--skip-images)."
fi

if [[ "$RUN_BACKEND_BENCH" -eq 1 ]]; then
  log "Running gauss vs closed-form benchmark..."
  run python3 "$ROOT_DIR/scripts/benchmark_backends.py" \
    --build-dir "$BUILD_DIR" \
    --out-dir "$OUT_DIR/benchmark" \
    --suite "$BACKEND_BENCH_SUITE" \
    --repeats "$BACKEND_BENCH_REPEATS"
else
  log "Skipping backend benchmark."
fi

log "Generating consolidated Markdown report..."
run python3 "$ROOT_DIR/scripts/generate_results_md.py" \
  --out-dir "$OUT_DIR" \
  --report "$OUT_DIR/RESULTS_REPORT.md"

log "Pipeline completed."
log "Main outputs:"
if [[ "$SKIP_VALIDATE" -eq 0 ]]; then
  printf '  - %s\n' "$OUT_DIR/validation/validation_2d_22.csv"
  if [[ -f "$OUT_DIR/validation/validation_2d_21.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_21.csv"
  fi
  if [[ -f "$OUT_DIR/validation/validation_2d_221.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_221.csv"
  fi
  if [[ -f "$OUT_DIR/validation/validation_2d_222_table6.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_222_table6.csv"
  fi
  if [[ -f "$OUT_DIR/validation/validation_2d_222_table7.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_222_table7.csv"
  fi
  if [[ -f "$OUT_DIR/validation/validation_2d_224_table10.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_224_table10.csv"
  fi
  if [[ -f "$OUT_DIR/validation/validation_2d_index.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_2d_index.csv"
  fi
  if [[ -f "$OUT_DIR/validation/VALIDATION_2D_INDEX.md" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/VALIDATION_2D_INDEX.md"
  fi
  if [[ -f "$OUT_DIR/validation/validation_3d_index.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_3d_index.csv"
  fi
  if [[ -f "$OUT_DIR/validation/VALIDATION_3D_INDEX.md" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/VALIDATION_3D_INDEX.md"
  fi
  if [[ -f "$OUT_DIR/validation/validation_master_index.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_master_index.csv"
  fi
  if [[ -f "$OUT_DIR/validation/VALIDATION_MASTER_INDEX.md" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/VALIDATION_MASTER_INDEX.md"
  fi
  if [[ -f "$OUT_DIR/validation/validation_verdict.csv" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_verdict.csv"
  fi
  if [[ -f "$OUT_DIR/validation/VALIDATION_VERDICT.md" ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/VALIDATION_VERDICT.md"
  fi
  if [[ "$SKIP_3D" -eq 0 ]]; then
    printf '  - %s\n' "$OUT_DIR/validation/validation_3d_31_modes.csv"
    printf '  - %s\n' "$OUT_DIR/validation/validation_3d_31_summary.csv"
  fi
fi
if [[ "$SKIP_IMAGES" -eq 0 ]]; then
  printf '  - %s\n' "$OUT_DIR/img_all/mode_summary.csv"
  printf '  - %s\n' "$OUT_DIR/img_all/"
  if [[ -d "$OUT_DIR/img_all/validation_2d_22" ]]; then
    printf '  - %s\n' "$OUT_DIR/img_all/validation_2d_22/"
  fi
fi
if [[ "$NO_LOG" -eq 0 ]]; then
  printf '  - %s\n' "$LOG_FILE"
fi
if [[ "$RUN_BACKEND_BENCH" -eq 1 ]]; then
  printf '  - %s\n' "$OUT_DIR/benchmark/backend_benchmark_detail.csv"
  printf '  - %s\n' "$OUT_DIR/benchmark/backend_benchmark_summary.csv"
  printf '  - %s\n' "$OUT_DIR/benchmark/BACKEND_BENCHMARK.md"
fi
printf '  - %s\n' "$OUT_DIR/RESULTS_REPORT.md"
