#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_SH="$ROOT_DIR/scripts/build_and_run_all.sh"

BACKEND_MODE="both"
OUT_ROOT="$ROOT_DIR/out/backend_compare"
INTERACTIVE=0
declare -a FORWARD_ARGS=()

usage() {
  cat <<'EOF'
Usage: scripts/run_backend_compare.sh [options] [-- <args for build_and_run_all.sh>]

Wrapper para rodar o pipeline em um ou dois backends, sempre em pastas
separadas para nao sobrescrever VTK/PNG/CSV.

Options:
  --backend-mode <m>   gauss|closed-form|both (default: both)
  --out-root <dir>     Diretorio raiz de saida (default: out/backend_compare)
  --interactive        Faz perguntas simples antes de disparar o pipeline
  --help               Mostra esta ajuda

Examples:
  scripts/run_backend_compare.sh --backend-mode both -- --case 2.2.4 --with-images
  scripts/run_backend_compare.sh --interactive
EOF
}

prompt_default() {
  local label="$1"
  local default="$2"
  local value
  read -r -p "$label [$default]: " value || true
  if [[ -z "${value:-}" ]]; then
    echo "$default"
  else
    echo "$value"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --backend-mode)
      BACKEND_MODE="$2"
      shift 2
      ;;
    --out-root)
      OUT_ROOT="$2"
      shift 2
      ;;
    --interactive)
      INTERACTIVE=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    --)
      shift
      FORWARD_ARGS+=("$@")
      break
      ;;
    *)
      FORWARD_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ "$INTERACTIVE" -eq 1 ]]; then
  BACKEND_MODE="$(prompt_default "Backend mode" "$BACKEND_MODE")"
  OUT_ROOT="$(prompt_default "Diretorio de saida" "$OUT_ROOT")"

  run_validate="$(prompt_default "Rodar validacoes? (yes/no)" "yes")"
  run_images="$(prompt_default "Gerar imagens? (yes/no)" "yes")"
  run_bench="$(prompt_default "Rodar benchmark de backends? (yes/no)" "no")"
  run_debug_blocks="$(prompt_default "Debug blocos locais? (yes/no)" "no")"
  run_debug_cands="$(prompt_default "Debug candidatos/raizes? (yes/no)" "no")"
  show_validate_output="$(prompt_default "Mostrar saida bruta das validacoes? (yes/no)" "no")"

  if [[ "${run_validate,,}" == "no" ]]; then
    FORWARD_ARGS+=("--skip-validate")
  fi
  if [[ "${run_images,,}" == "no" ]]; then
    FORWARD_ARGS+=("--skip-images")
  fi
  if [[ "${run_bench,,}" == "yes" ]]; then
    bench_suite="$(prompt_default "Suite do benchmark" "core")"
    bench_repeats="$(prompt_default "Repeticoes do benchmark" "3")"
    FORWARD_ARGS+=("--benchmark-backends" "--benchmark-suite" "$bench_suite" "--benchmark-repeats" "$bench_repeats")
  fi
  if [[ "${run_debug_blocks,,}" == "yes" ]]; then
    FORWARD_ARGS+=("--debug-local-blocks")
  fi
  if [[ "${run_debug_cands,,}" == "yes" ]]; then
    FORWARD_ARGS+=("--debug-candidates")
  fi
  if [[ "${show_validate_output,,}" == "yes" ]]; then
    FORWARD_ARGS+=("--show-validation-output")
  fi
fi

case "$BACKEND_MODE" in
  gauss|closed-form|both) ;;
  *)
    echo "Backend mode invalido: $BACKEND_MODE (use gauss|closed-form|both)" >&2
    exit 2
    ;;
esac

run_backend() {
  local backend="$1"
  local safe_backend="${backend//-/_}"
  local out_dir="$OUT_ROOT/$safe_backend"
  echo "[backend-compare] backend=$backend out_dir=$out_dir"
  "$PIPELINE_SH" --out-dir "$out_dir" --backend "$backend" "${FORWARD_ARGS[@]}"
}

mkdir -p "$OUT_ROOT"

if [[ "$BACKEND_MODE" == "both" ]]; then
  run_backend "gauss"
  run_backend "closed-form"
else
  run_backend "$BACKEND_MODE"
fi

echo "[backend-compare] finalizado."
