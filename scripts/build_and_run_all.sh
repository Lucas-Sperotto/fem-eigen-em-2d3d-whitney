#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"
OUT_DIR="$ROOT_DIR/out"
BUILD_TYPE="Release"
PROFILE="quick"
MODE_EXPORT=8
SKIP_VALIDATE=0
SKIP_IMAGES=0
SKIP_3D=0
VERBOSE=0
NO_LOG=0
FORCE_VALIDATE=0
FORCE_IMAGES=0
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
- generate all quiver images and mode summary CSV

Options:
  --build-dir <dir>      Build directory (default: build)
  --out-dir <dir>        Root output directory (default: out)
  --build-type <type>    CMake build type (default: Release)
  --jobs <N>             Parallel build jobs (default: nproc)
  --profile <quick|full> Profile for validate_3d_31.py (default: quick)
  --mode-export <N>      Number of exported 2D modes per TE/TM block (default: 8)
  --case <id>            Run only selected section/case (repeatable)
  --log-file <path>      Write console log to file (default: <out-dir>/run_all.log)
  --no-log               Disable log file output
  --skip-3d              Skip fem3d0/fem3d1 runs and 3D validation
  --skip-validate        Skip validate_2d_22.py and validate_3d_31.py
  --with-validate        Force validations even in --case mode
  --skip-images          Skip plot_vtk_quiver.py --all-img
  --with-images          Force image generation even in --case mode
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
if [[ "$NO_LOG" -eq 0 ]]; then
  log "LOG_FILE=$LOG_FILE"
fi
if [[ "$CASE_MODE" -eq 1 ]]; then
  log "CASE_MODE=1 selected: ${CASES[*]}"
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
      run ./helm10_rect 14 14 "$MODE_EXPORT"
      run ./helm10_circle 10 48 "$MODE_EXPORT"
      run ./helm10_coax 10 48 "$MODE_EXPORT"
    fi

    if [[ "$RUN_221" -eq 1 ]]; then
      run ./edge_rect 14 14 "$MODE_EXPORT"
      run ./edge_circle 10 48 "$MODE_EXPORT"
      run ./edge_coax 10 48 "$MODE_EXPORT"
    fi

    if [[ "$RUN_222" -eq 1 ]]; then
      run ./mixed_rect 12 6
      run ./mixed_circle 10 48
      run ./mixed_coax 10 48
    fi

    if [[ "$RUN_223" -eq 1 ]]; then
      run ./helmvec2_rect 10 6 6
    fi

    if [[ "$RUN_224" -eq 1 ]]; then
      run ./helmvec3_rect 0.20 10 5
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
        run ./fem3d0_rect --air
        run ./fem3d1_rect --air
      fi
      if [[ "$RUN_3D_HALF" -eq 1 ]]; then
        run ./fem3d0_rect --half
        run ./fem3d1_rect --half
      fi
      if [[ "$RUN_3D_CYL" -eq 1 ]]; then
        run ./fem3d0_rect --cyl
        run ./fem3d1_rect --cyl
      fi
      if [[ "$RUN_3D_SPHERE" -eq 1 ]]; then
        run ./fem3d0_rect --sphere
        run ./fem3d1_rect --sphere
      fi
    )
  else
    log "No 3D case selected."
  fi
else
  log "Skipping 3D executables (--skip-3d)."
fi

if [[ "$SKIP_VALIDATE" -eq 0 ]]; then
  log "Running 2D validation..."
  run python3 "$ROOT_DIR/scripts/validate_2d_22.py" \
    --build-dir "$BUILD_DIR" \
    --out-csv "$OUT_DIR/validation/validation_2d_22.csv"

  if [[ "$SKIP_3D" -eq 0 ]]; then
    log "Running 3D validation..."
    run python3 "$ROOT_DIR/scripts/validate_3d_31.py" \
      --profile "$PROFILE" \
      --solver both \
      --build-dir "$BUILD_DIR" \
      --out-modes "$OUT_DIR/validation/validation_3d_31_modes.csv" \
      --out-summary "$OUT_DIR/validation/validation_3d_31_summary.csv"
  else
    log "Skipping 3D validation because --skip-3d is active."
  fi
else
  log "Skipping validations (--skip-validate)."
fi

if [[ "$SKIP_IMAGES" -eq 0 ]]; then
  log "Generating all images and mode summary CSV..."
  run python3 "$ROOT_DIR/scripts/plot_vtk_quiver.py" \
    --all-img \
    --build-dir "$BUILD_DIR" \
    --vtk-root "$OUT_DIR/2d" \
    --out-dir "$OUT_DIR/img_all" \
    --csv "$OUT_DIR/img_all/mode_summary.csv" \
    --mode-export "$MODE_EXPORT" \
    --max-rank "$MODE_EXPORT"

  if [[ -f "$OUT_DIR/validation/validation_2d_22.csv" ]]; then
    log "Generating 2.2.2/2.2.3/2.2.4 validation figures..."
    run python3 "$ROOT_DIR/scripts/plot_validation_2d_22.py" \
      --in-csv "$OUT_DIR/validation/validation_2d_22.csv" \
      --out-dir "$OUT_DIR/img_all/validation_2d_22"
  else
    log "Skipping 2.2.x validation figures (missing $OUT_DIR/validation/validation_2d_22.csv)."
  fi
else
  log "Skipping image generation (--skip-images)."
fi

log "Generating consolidated Markdown report..."
run python3 "$ROOT_DIR/scripts/generate_results_md.py" \
  --out-dir "$OUT_DIR" \
  --report "$OUT_DIR/RESULTS_REPORT.md"

log "Pipeline completed."
log "Main outputs:"
if [[ "$SKIP_VALIDATE" -eq 0 ]]; then
  printf '  - %s\n' "$OUT_DIR/validation/validation_2d_22.csv"
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
printf '  - %s\n' "$OUT_DIR/RESULTS_REPORT.md"
