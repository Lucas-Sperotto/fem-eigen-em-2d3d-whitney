# NASA TP-3485 reproduction (FEM eigenvalue EM)

Reproducao numerica da NASA Technical Paper 3485:
*Finite Element Method for Eigenvalue Problems in Electromagnetics* (1994).

Este repositorio implementa, valida e organiza os blocos 2D e 3D do artigo,
com foco em:
- formulacao didatica,
- comparacao com tabelas de referencia,
- fluxo reprodutivel via executaveis C++ e scripts Python.

## 1) Mapa do artigo para o codigo

- Sec. 2.1 (escalar 2D, cutoff `kc`): `src/helm10`
- Sec. 2.2.1 (vetorial transversal, edge 2D, cutoff `kc`): `src/helmvec`
- Sec. 2.2.2 (misto vetorial+escalar, cutoff `kc`): `src/helmvec1`
- Sec. 2.2.3 (`k0` dado `beta`): `src/helmvec2`
- Sec. 2.2.4 (`beta` dado `k0`): `src/helmvec3`
- Sec. 3.1 (cavidades 3D, edge tetra): `src/fem3d0` e `src/fem3d1`

## 2) Estrutura principal (src)

- `src/core`: malhas, estruturas de matriz, solver LAPACK, utilitarios de I/O VTK
- `src/edge`: base/DOFs/montagem edge 2D
- `src/edge3d`: base/DOFs/montagem edge 3D
- `src/helm10`: executaveis escalares 2D
- `src/helmvec`: executaveis vetoriais transversais 2D
- `src/helmvec1`: executaveis mistos para `kc`
- `src/helmvec2`: sistema acoplado para `k0` com `beta` dado
- `src/helmvec3`: sistema acoplado para `beta` com `k0` dado
- `src/fem3d`: utilitarios compartilhados de casos/tabelas 3D
- `src/fem3d0`: solver 3D denso
- `src/fem3d1`: solver 3D com montagem esparsa (solve denso fallback)

## 3) Dependencias e build

Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libopenblas-dev liblapacke-dev python3 python3-pip
```

Build:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
```

## 4) Executaveis 2D

### 4.1) Secao 2.1 (`helm10`)

```bash
./build/helm10_rect 14 14 8
./build/helm10_circle 10 48 8
./build/helm10_coax 10 48 8
```

Saidas tipicas:
- lista de `kc`
- tabela de comparacao FEM x analitico com correlacao modal (`rho`)
- VTK em `out/2d/2.1_scalar/{rect,circle,coax}` (inclui varios modos por rank)

### 4.2) Secao 2.2.1 (`helmvec`)

```bash
./build/edge_rect 14 14 8
./build/edge_circle 10 48 8
./build/edge_coax 10 48 8
```

Saidas tipicas:
- lista de `kc`
- tabela FEM x analitico (matching por correlacao em massa)
- VTK em `out/2d/2.2.1_edge/{rect,circle,coax}` (inclui varios modos por rank)

### 4.3) Secao 2.2.2 (`helmvec1`)

```bash
./build/mixed_rect 12 6
./build/mixed_circle 10 48
./build/mixed_coax 10 48
```

Saidas tipicas:
- espectros separados por energia de bloco (edge vs escalar)
- comparacao analitica detalhada no caso retangular

### 4.4) Secao 2.2.3 (`helmvec2`)

```bash
./build/helmvec2_rect 10 6 6
# args: beta nx ny [debug]
```

Saidas tipicas:
- Tabela 8 (Figura 11): `k0L(FEM matched)` vs HELMVEC2/Hayata

### 4.5) Secao 2.2.4 (`helmvec3`)

```bash
./build/helmvec3_rect 0.20 10 5
./build/helmvec3_rect 0.20 10 5 1
# args: d_over_a nx ny [debug]
```

Saidas tipicas:
- Tabela 9 (Figura 12)
- preview de ramo para Figura 13
- validacao da Tabela 10 (Figura 13)

## 5) Executaveis 3D (Secao 3.1)

### 5.1) FEM3D0 (denso)

```bash
./build/fem3d0_rect
./build/fem3d0_rect --half
./build/fem3d0_rect --cyl
./build/fem3d0_rect --sphere
./build/fem3d0_rect --all --nx 5 --ny 4 --nz 4
```

### 5.2) FEM3D1 (montagem esparsa)

```bash
./build/fem3d1_rect --air
./build/fem3d1_rect --half
./build/fem3d1_rect --cyl
./build/fem3d1_rect --sphere
./build/fem3d1_rect --all --nx 6 --ny 4 --nz 4
```

## 6) Scripts de validacao

### 6.1) Validacao 2D (Secao 2.2)

```bash
python3 scripts/validate_2d_22.py
```

CLI util:

```bash
python3 scripts/validate_2d_22.py \
  --build-dir build \
  --out-csv out/validation/validation_2d_22.csv \
  --rect-nx 12 --rect-ny 6 \
  --circle-nr 10 --circle-nt 48 \
  --coax-nr 10 --coax-nt 48 \
  --beta 10 --hv2-nx 6 --hv2-ny 6 \
  --d-over-a 0.20 --hv3-nx 10 --hv3-ny 5 \
  --verbose
```

Saida:
- `out/validation/validation_2d_22.csv`

### 6.2) Validacao 3D (Secao 3.1)

```bash
python3 scripts/validate_3d_31.py --profile quick --solver both
python3 scripts/validate_3d_31.py --profile full --solver fem3d1
```

CLI util:

```bash
python3 scripts/validate_3d_31.py \
  --profile quick \
  --solver both \
  --cases air,half,cyl,sphere \
  --build-dir build \
  --out-modes out/validation/validation_3d_31_modes.csv \
  --out-summary out/validation/validation_3d_31_summary.csv \
  --verbose
```

Saidas:
- `out/validation/validation_3d_31_modes.csv`
- `out/validation/validation_3d_31_summary.csv`

## 7) Plot de campos VTK (quiver)

Modo arquivo unico:

```bash
python3 scripts/plot_vtk_quiver.py out/2d/2.1_scalar/rect/te10_rect_sv.vtk --out out/img/te10_rect.png --stride 2 --scale 22 --dpi 210
python3 scripts/plot_vtk_quiver.py out/2d/2.1_scalar/rect/tm11_rect_sv.vtk --out out/img/tm11_rect.png --stride 2 --scale 22 --dpi 210
python3 scripts/plot_vtk_quiver.py out/2d/2.2.1_edge/rect/edge_rect_Et.vtk --out out/img/edge_rect_Et.png --stride 2 --scale 25 --dpi 210
```

Modo lote (gera imagens e CSV resumo):

```bash
python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out/2d --out-dir out/img_all --csv out/img_all/mode_summary.csv --mode-export 8 --max-rank 8
```

Imagens de validacao para os casos 2.2.2, 2.2.3 e 2.2.4:

```bash
python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22
```

Observacao:
- `--all_img` continua aceito como alias de compatibilidade.

Saidas do lote:
- imagens em `out/img_all/` preservando a arvore de `out/2d/`
- `out/img_all/mode_summary.csv`
- `out/img_all/validation_2d_22/` (graficos de 2.2.2/2.2.3/2.2.4)

## 8) Fluxo recomendado de reproducao

1. Compilar (`cmake --build . -j`).
2. Rodar casos 2D-base (`helm10_*`, `edge_*`, `mixed_*`).
3. Rodar acoplados (`helmvec2_rect`, `helmvec3_rect`).
4. Rodar validacao 2D (`validate_2d_22.py`).
5. Rodar validacao 3D (`validate_3d_31.py`).
6. Gerar figuras (`plot_vtk_quiver.py --all-img ...`).

## 9) Notas numericas

- Matching modal usa correlacao e tratamento de degenerescencia para reduzir troca artificial de ordem.
- Em problemas generalizados nao simetricos (`helmvec2`, `helmvec3`), o pipeline filtra raizes nao fisicas (parte imaginaria, sinal e faixa fisica).
- `fem3d1` usa montagem esparsa simetrica; o solve atual ainda converte para denso antes de `dsygv`.

## 10) Guia de reproducao do paper (checklist)

Use esta sequencia para reproduzir os blocos numericos em ordem de tabelas/figuras.

1. Tabela 1 (retangular escalar, Sec. 2.1):
```bash
./build/helm10_rect 14 14 8
```

2. Tabela 2 (circular escalar, Sec. 2.1):
```bash
./build/helm10_circle 10 48 8
```

3. Tabela 3 (coax escalar, Sec. 2.1):
```bash
./build/helm10_coax 10 48 8
```

4. Figuras de campo vetorial 2D (edge, Sec. 2.2.1):
```bash
./build/edge_rect 14 14 8
./build/edge_circle 10 48 8
./build/edge_coax 10 48 8
```

5. Sistema misto no cutoff (Sec. 2.2.2):
```bash
./build/mixed_rect 12 6
./build/mixed_circle 10 48
./build/mixed_coax 10 48
```

6. Figura 11 / Tabela 8 (`k0` dado `beta`, Sec. 2.2.3):
```bash
./build/helmvec2_rect 10 6 6
```

7. Figura 12 / Tabela 9 e Figura 13 / Tabela 10 (`beta` dado `k0`, Sec. 2.2.4):
```bash
./build/helmvec3_rect 0.20 10 5
```

8. Secao 3.1 em cavidades 3D (Tabelas 12-15):
```bash
./build/fem3d0_rect --all
./build/fem3d1_rect --all
```

9. Validacao automatica consolidada:
```bash
python3 scripts/validate_2d_22.py --build-dir build --out-csv out/validation/validation_2d_22.csv
python3 scripts/validate_3d_31.py --profile quick --solver both --build-dir build --out-modes out/validation/validation_3d_31_modes.csv --out-summary out/validation/validation_3d_31_summary.csv
```

10. Geracao de imagens e CSV de modos:
```bash
python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out/2d --out-dir out/img_all --csv out/img_all/mode_summary.csv --mode-export 8 --max-rank 8
```

## 11) Script unico (compila e roda tudo)

Foi adicionado:
- `scripts/build_and_run_all.sh`

Ele executa:
- configuracao e build CMake,
- todos os executaveis 2D,
- todos os executaveis 3D (`fem3d0` e `fem3d1`),
- validacoes 2D e 3D,
- geracao de imagens e `mode_summary.csv`,
- log automatico em arquivo (`out/run_all.log` por padrao).

Uso padrao:

```bash
./scripts/build_and_run_all.sh
```

Uso com perfil completo de validacao 3D:

```bash
./scripts/build_and_run_all.sh --profile full
```

Controlar quantos modos 2D sao exportados e plotados:

```bash
./scripts/build_and_run_all.sh --mode-export 10
```

Log customizado:

```bash
./scripts/build_and_run_all.sh --log-file out/logs/full_run.log
./scripts/build_and_run_all.sh --no-log
```

Execucao seletiva por secao/tabela (pode repetir `--case`):

```bash
# So secao 2.2.3 (Figura 11 / Tabela 8)
./scripts/build_and_run_all.sh --case 2.2.3

# So Tabela 13 (Figura 16, cavidade half-filled)
./scripts/build_and_run_all.sh --case table13

# Combinar blocos
./scripts/build_and_run_all.sh --case 2.2.4 --case table15
```

Observacao:
- em modo `--case`, validacao e geracao de imagens ficam desativadas por padrao
  (para evitar falha por faltarem casos nao executados);
- para forcar:

```bash
./scripts/build_and_run_all.sh --case 2.2.3 --with-validate
./scripts/build_and_run_all.sh --case 2.2.1 --with-images
```

Opcoes uteis:

```bash
./scripts/build_and_run_all.sh --help
./scripts/build_and_run_all.sh --build-dir build --jobs 8 --build-type Release
./scripts/build_and_run_all.sh --out-dir out
./scripts/build_and_run_all.sh --case 2.1
./scripts/build_and_run_all.sh --case 3.1
./scripts/build_and_run_all.sh --skip-3d
./scripts/build_and_run_all.sh --skip-validate
./scripts/build_and_run_all.sh --skip-images
./scripts/build_and_run_all.sh --verbose
```
