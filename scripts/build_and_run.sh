#!/bin/bash

set -e  # para parar se der erro

# Ir para a raiz do projeto (script está em scripts/)
cd "$(dirname "$0")/.."

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

echo "⚙️ Compilando helm..."
cmake --build build --target helm10_rect helm10_circle helm10_coax -j$(nproc)

echo "🚀 Executando helm..."
./build/helm10_rect 1.0 13 28 20 --backend closed-form
./build/helm10_circle 10 40 20 --backend closed-form
./build/helm10_coax 10 40 20 --backend closed-form

echo "📊 Plotando helm..."
python3 scripts/helm10.py

echo "⚙️ Compilando edge..."
cmake --build build --target edge_rect edge_circle edge_coax -j$(nproc)

echo "🚀 Executando edge..."
./build/edge_rect 13 28 20 --backend closed-form
./build/edge_circle 10 40 20 --backend closed-form
./build/edge_coax 10 40 20 --backend closed-form

echo "📊 Plotando edge..."
python3 scripts/helmvec.py

echo "⚙️ Compilando mixed..."
cmake --build build --target mixed_rect mixed_circle mixed_coax -j$(nproc)

echo "🚀 Executando mixed..."
./build/mixed_rect 13 28 --backend closed-form
./build/mixed_circle 10 40 --backend closed-form
./build/mixed_coax 10 40 --backend closed-form

echo "📊 Plotando mixed..."
python3 scripts/helmvec1.py

echo "✅ Tudo finalizado!"
