#!/bin/bash

set -euo pipefail

# Helper legado/local para validacao e reproducao local das saidas do artigo.
# Nao e o fluxo publico/canonico do repositorio.
# Pipeline recomendado: scripts/build_and_run_all.sh

# Ir para a raiz do projeto (script está em scripts/)
cd "$(dirname "$0")/.."

printf '%s\n' \
    '[legacy-helper] scripts/build_and_run.sh e um helper legado/local para validacao do artigo.' \
    '[legacy-helper] Ele nao e a interface canonica; use scripts/build_and_run_all.sh para o pipeline publico.' \
    '[legacy-helper] Aviso: este script remove build/ e out/ antes de recompilar.' >&2

echo "🧹 Limpando diretórios antigos..."

# Segurança: só remove se as pastas existirem
if [ -d "build" ]; then
    rm -rf build
fi

if [ -d "out" ]; then
    rm -rf out
fi

# (opcional) recriar out
mkdir -p out

echo "🔧 Configurando projeto..."
cmake -S . -B build

JOBS="$(nproc)"

echo "⚙️ Compilando a suíte didática 2D + 3D..."
cmake --build build \
    --target \
    helm10_rect helm10_circle helm10_coax \
    edge_rect edge_circle edge_coax \
    mixed_rect mixed_circle mixed_coax \
    helmvec2_rect \
    helmvec3_fig12_rect helmvec3_fig13_rect \
    fem3d0_air fem3d0_half fem3d0_cyl fem3d0_sphere \
    fem3d1_air fem3d1_half fem3d1_cyl fem3d1_sphere \
    -j"${JOBS}"

echo "🚀 Executando HELM10 com a discretização de referência do artigo..."
./build/helm10_rect 1.0 10 20 10 --backend closed-form #400 elementos
./build/helm10_circle 8 15 10 --backend closed-form #200 elementos
./build/helm10_coax 10 17 10 --backend closed-form #340 elementos

echo "📊 Gerando imagens do HELM10..."
python3 scripts/helm10.py

echo "🚀 Executando HELMVEC com a discretização de referência do artigo..."
./build/edge_rect 10 20 10 --backend closed-form #400 elementos
./build/edge_circle 8 15 10 --backend closed-form #200 elementos
./build/edge_coax 10 17 10 --backend closed-form #340 elementos

echo "📊 Gerando imagens do HELMVEC..."
python3 scripts/helmvec.py

echo "🚀 Executando HELMVEC1 com a discretização de referência do artigo..."
./build/mixed_rect --nx 10 --ny 20 --backend closed-form #400 elementos
./build/mixed_circle --nr 8 --nt 15 --backend closed-form #200 elementos
./build/mixed_coax --nr 10 --nt 17 --backend closed-form #340 elementos

echo "📊 Gerando imagens do HELMVEC1..."
python3 scripts/helmvec1.py

echo "🚀 Executando HELMVEC2 com a discretização de referência do artigo..."
./build/helmvec2_rect --beta 10 --nx 20 --ny 20 --backend closed-form #beta=10 com malha retangular 20x20

echo "📊 Gerando imagens do HELMVEC2..."
python3 scripts/helmvec2.py

echo "🚀 Executando HELMVEC3 com a discretização de referência do artigo..."
./build/helmvec3_fig12_rect 10 5 --backend closed-form #100 elementos
./build/helmvec3_fig13_rect 0.20 10 5 --backend closed-form #100 elementos

echo "📊 Gerando imagens do HELMVEC3..."
python3 scripts/helmvec3.py

echo "🚀 Executando FEM3D0/FEM3D1 com as discretizações de referência do artigo..."
./build/fem3d0_air --backend closed-form #343 elementos
./build/fem3d1_air --backend closed-form #343 elementos
./build/fem3d0_half --backend closed-form #615 elementos
./build/fem3d1_half --backend closed-form #615 elementos
./build/fem3d0_cyl --backend closed-form #633 elementos
./build/fem3d1_cyl --backend closed-form #633 elementos
./build/fem3d0_sphere --backend closed-form #473 elementos
./build/fem3d1_sphere --backend closed-form #473 elementos

echo "✅ Tudo finalizado!"
