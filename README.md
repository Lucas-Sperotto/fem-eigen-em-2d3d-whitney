# NASA TP-3485 reproduction (FEM eigenvalue EM)

Reproducao numerica da NASA Technical Paper 3485:
*Finite Element Method for Eigenvalue Problems in Electromagnetics* (1994).

Este repositorio implementa, valida e organiza os blocos 2D e 3D do artigo,
com foco em:

- formulacao didatica,
- comparacao com tabelas de referencia,
- fluxo reprodutivel via executaveis C++ e scripts Python.

## 0) Trilha teorica e documentacao

Antes de entrar no codigo, a trilha teorica revisada do artigo esta em:

- [docs/README.md](docs/README.md): indice geral da documentacao, em ordem de estudo.
- [docs/19950011772.pdf](docs/19950011772.pdf): PDF original do paper.
- [docs/Rastreabilidade_Equacoes_Artigo_Codigo.md](docs/Rastreabilidade_Equacoes_Artigo_Codigo.md): trilha central de rastreabilidade entre equacoes do artigo e funcoes/arquivos C++.
- [docs/Tabela_Executaveis_Entradas_Saidas.md](docs/Tabela_Executaveis_Entradas_Saidas.md): tabela unificada com todos os executaveis gerados, suas entradas e suas saidas.
- [docs/results/README.md](docs/results/README.md): resultados curados, figuras e validacoes preservadas no repositorio.

Os arquivos principais em `docs/` preservam a traducao-base do artigo e receberam comentarios complementares, explicacoes intermediarias entre equacoes e notas de consistencia para estudo.

## 1) Mapa do artigo para o codigo

- Sec. 2.1 (escalar 2D, cutoff `kc`): `src/helm10`
- Sec. 2.2.1 (vetorial transversal, edge 2D, cutoff `kc`): `src/helmvec`
- Sec. 2.2.2 (misto vetorial+escalar, cutoff `kc`): `src/helmvec1`
- Sec. 2.2.3 (`k0` dado `beta`): `src/helmvec2`
- Sec. 2.2.4 (`beta` dado `k0`): `src/helmvec3`
- Sec. 3.1 (cavidades 3D, edge tetra): `src/fem3d0` e `src/fem3d1`
- Fachada didatica publica por equacao global: `src/article`

## 1.1) Mapa do Apendice FORTRAN para o repositorio

| Programa do artigo | Secao | Equacao global | Modulo atual | Arquivo principal de montagem |
|---|---|---:|---|---|
| `HELM10` | 2.1 | Eq. (43) | `src/helm10` | `src/core/helm10_scalar_system.cpp` |
| `HELMVEC` | 2.2.1 | Eq. (65) | `src/helmvec` | `src/edge/edge_assembly.cpp` |
| `HELMVEC1` | 2.2.2 | Eq. (92) | `src/helmvec1` | `src/helmvec1/helmvec1_mixed_system.cpp` |
| `HELMVEC2` | 2.2.3 | Eq. (119) | `src/helmvec2` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| `HELMVEC3` | 2.2.4 | Eq. (136) | `src/helmvec3` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| `FEM3D0` | 3.1 | Eq. (178) | `src/fem3d0` | `src/edge3d/edge3d_assembly.cpp` |
| `FEM3D1` | 3.1 | Eq. (178) | `src/fem3d1` | `src/edge3d/edge3d_assembly.cpp` |

Observacao:

- os `main_*` fazem o papel dos drivers dos programas do apendice;
- as montagens globais efetivas estao concentradas nos arquivos listados na
  ultima coluna;
- as formas fechadas locais ficam organizadas em `src/explicit`.
- as entradas publicas mais didaticas, nomeadas pelas equacoes globais do
  artigo, ficam em `src/article/tp3485_systems.hpp`.

## 1.2) Resumo matematico do repositorio

O repositorio cobre seis problemas generalizados principais do artigo:

### Notacao comum

Ao longo dos modulos, a notacao base e sempre esta:

```text
Gamma    = secao transversal 2D do guia
T        = triangulo 2D de area A
Tet      = tetraedro 3D de volume V
N_i      = funcao de forma nodal linear em 2D
lambda_i = coordenada simplex linear em 2D/3D
W_m      = base de Whitney associada a aresta local m
L_m      = comprimento da aresta local m
eps_r    = permissividade relativa por elemento
mu_r     = permeabilidade relativa por elemento
```

Nos triangulos 2D:

```text
N_i(x,y)      = (a_i + b_i x + c_i y) / (2A)
grad N_i      = [b_i, c_i] / (2A)
x_tri, y_tri  = coordenadas do baricentro do triangulo
```

Nos tetraedros 3D:

```text
lambda_i(x,y,z) = (a_i + b_i x + c_i y + d_i z) / (6V)
grad lambda_i   = [b_i, c_i, d_i] / (6V)
x_tet, y_tet, z_tet = coordenadas do baricentro do tetraedro
```

Nas formulacoes vetoriais 2D e 3D:

```text
W_m = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
```

Por isso, quando um README fala em `St`, `Tt`, `Sz`, `Tz`, `C`, `Mt^(1/mu)` ou
`Gz`, ele esta apenas reorganizando integrais construidas a partir desses
blocos geometricos comuns.

### Secao 2.1 - escalar 2D

```text
S phi = kc^2 T phi
```

com, por elemento triangular linear:

```text
N_i(x,y) = (a_i + b_i x + c_i y) / (2A)
grad N_i = [b_i, c_i] / (2A)

S_e(i,j) = (1/mu_r) * (b_i b_j + c_i c_j) / (4A)
T_e(i,j) = eps_r * (A/12) * [2 1 1; 1 2 1; 1 1 2]_(i,j)
```

### Secao 2.2.1 - edge 2D transversal

```text
S e = kc^2 T e
```

com base de Whitney:

```text
W_m = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
lambda_i = (a_i + b_i x + c_i y) / (2A)
```

e, na forma fechada:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(4 A^3) * D_m D_n
T_e(m,n) = eps_r * (L_m L_n)/(16 A^3) * sum_{k=1}^5 It_k
```

### Secao 2.2.2 - sistema misto para `kc`

```text
[ St   0 ] [et] = kc^2 [ Tt   0 ] [et]
[  0  Sz ] [ez]        [  0  Tz ] [ez]
```

onde:

```text
St <- bloco edge 2D
Tt <- massa edge 2D
Sz <- grad-grad escalar
Tz <- massa escalar
```

### Secao 2.2.3 - `k0` dado `beta`

```text
A x = k0^2 B x
x = [Et; Ez]
```

com:

```text
A_tt = St + beta^2 Mt^(1/mu)
A_tz = beta^2 C
A_zt = beta^2 C^T
A_zz = beta^2 Sz

B_tt = Tt
B_zz = beta^2 Tz
```

### Secao 2.2.4 - `beta` dado `k0`

```text
P x = beta^2 Q x
x = [Et; Ez]
```

com:

```text
P_tt = St - k0^2 Tt
P_zz = k0^2 Tz

Q_tt = -Mt^(1/mu)
Q_tz = C
Q_zt = C^T
Q_zz = Sz
```

### Secao 3.1 - edge 3D tetraedrico

```text
S e = k0^2 T e
```

com coordenadas simplex tetraedricas:

```text
lambda_i(x,y,z) = (a_i + b_i x + c_i y + d_i z) / (6V)
W_m = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
```

e, na forma fechada usada no repositorio:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(324 V^3) * (...)
T_e(m,n) = eps_r * (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

Os detalhes de cada bloco, incluindo coeficientes geometricos e relacao com as
equacoes numeradas do artigo, estao nos READMEs de cada modulo.

## 2) Estrutura principal (src)

- `src/core`: malhas, estruturas de matriz, solver LAPACK, utilitarios de I/O VTK
- `src/article`: fachada publica nomeada pelas equacoes globais do artigo
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

## 2.1) Documentacao por modulo (READMEs)

- [src/helm10/README.md](src/helm10/README.md): formulacao escalar 2D (Secao 2.1).
- [src/helmvec/README.md](src/helmvec/README.md): formulacao vetorial transversal com elementos de aresta (Secao 2.2.1).
- [src/helmvec1/README.md](src/helmvec1/README.md): sistema misto vetorial + escalar para `kc` (Secao 2.2.2).
- [src/helmvec2/README.md](src/helmvec2/README.md): sistema acoplado para obter `k0` com `beta` dado (Secao 2.2.3).
- [src/helmvec3/README.md](src/helmvec3/README.md): sistema acoplado para obter `beta` com `k0` dado (Secao 2.2.4).
- [docs/Artefatos_Espectrais_CSV_Referencia.md](docs/Artefatos_Espectrais_CSV_Referencia.md): formato dos CSVs em `linop/`, com matrizes em CRS, autovalores ordenados e autovetores ordenados das familias `HELM10` a `FEM3D1`.
- [docs/FEM3D_CSV_Referencia.md](docs/FEM3D_CSV_Referencia.md): referencia operacional dos artefatos 3D (`run.log`, `run_timing.csv`, `modes.csv`, `fields.csv`, `vtk/` e `linop/`) por caso de `FEM3D0` e `FEM3D1`.
- [docs/FEM3D_Imagens_Referencia.md](docs/FEM3D_Imagens_Referencia.md): guia do script `scripts/fem3d.py`, com resumos modais e projecoes ortogonais `XY/XZ/YZ` dos campos 3D exportados por modo.
- [src/explicit/README.md](src/explicit/README.md): backend `closed-form`, mapeamento de equacoes locais e rearranjos.
- [src/fem3d/README.md](src/fem3d/README.md): infraestrutura comum de validacao 3D (Secao 3.1).
- [src/fem3d0/README.md](src/fem3d0/README.md): solver 3D denso (`FEM3D0`).
- [src/fem3d1/README.md](src/fem3d1/README.md): solver 3D com montagem esparsa (`FEM3D1`).

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

## 3.1) Backends e depuracao

Os executaveis principais aceitam, quando aplicavel:

- `--backend closed-form`
- `--backend gauss`
- `--debug-local-blocks`
- `--debug-candidates`
- `--debug` ou `--debug-all`

Em termos práticos:

- `closed-form` usa as formas fechadas ligadas diretamente as equacoes do artigo;
- `closed-form` e o fluxo principal do repositorio;
- `gauss` preserva a montagem por quadratura/cubatura para verificacao auxiliar;
- `--debug-local-blocks` imprime o primeiro elemento local com rastreabilidade
  matematica;
- `--debug-candidates` imprime as primeiras raizes/candidatos antes do matching.

Exemplos:

```bash
./build/helm10_rect 1.0 14 14 8 --backend closed-form --debug-local-blocks
./build/edge_rect 14 14 8 --debug-candidates
./build/mixed_rect 12 6 --backend closed-form --debug-local-blocks --debug-candidates
./build/helmvec2_rect 10 6 6 --backend closed-form --debug-local-blocks
./build/helmvec3_fig12_rect 10 5 --backend closed-form --debug-candidates
./build/helmvec3_fig13_rect 0.20 10 5 --backend closed-form --debug-candidates
./build/fem3d0_air --backend closed-form --debug-local-blocks
```

## 4) Executaveis 2D

### 4.1) Secao 2.1 (`helm10`)

```bash
./build/helm10_rect 1.0 14 14 8
./build/helm10_circle 10 48 8
./build/helm10_coax 10 48 8
```

Saidas tipicas:

- lista de `kc`
- tabela de comparacao FEM x analitico com correlacao modal (`rho`)
- VTK em `out/helm10/{rect,circle,coax}/vtk`
- CSVs didaticos em `out/helm10/{rect,circle,coax}/csv`, com `modes.csv` e
  `fields_<modo>.csv`

### 4.2) Secao 2.2.1 (`helmvec`)

```bash
./build/edge_rect 14 14 8
./build/edge_circle 10 48 8
./build/edge_coax 10 48 8
```

Saidas tipicas:

- lista de `kc`
- tabela FEM x analitico (matching por correlacao em massa)
- `run.log` e `run_timing.csv` em `out/helmvec/{rect,circle,coax}`
- CSVs didaticos em `out/helmvec/{rect,circle,coax}/csv`, com `modes.csv` e
  `fields_<modo>.csv`
- VTK em `out/helmvec/{rect,circle,coax}/vtk`
- imagens geradas por `python3 scripts/helmvec.py` em
  `out/helmvec/{rect,circle,coax}/img`

### 4.3) Secao 2.2.2 (`helmvec1`)

```bash
./build/mixed_rect 12 6
./build/mixed_circle 10 48
./build/mixed_coax 10 48
```

Saidas tipicas:

- espectros separados por energia de bloco (edge vs escalar)
- matching analitico por correlacao de massa no bloco dominante
- `run.log` e `run_timing.csv` em `out/helmvec1/{rect,circle,coax}`
- CSVs didaticos em `out/helmvec1/{rect,circle,coax}/csv`, com
  `mixed_<caso>_modes.csv`, incluindo `rho_abs`, `match_space` e
  `match_method`
- imagens-resumo geradas por `python3 scripts/helmvec1.py` em
  `out/helmvec1/{rect,circle,coax}/img`, incluindo cutoff, `rho`, energia
  dominante e energias de bloco

### 4.4) Secao 2.2.3 (`helmvec2`)

```bash
./build/helmvec2_rect 10 6 6
# args: beta nx ny [debug]
```

Saidas tipicas:

- Tabela 8 (Figura 11): `k0L(FEM matched)` vs HELMVEC2/Hayata
- `run.log` e `run_timing.csv` em `out/helmvec2/rect`
- CSVs didaticos em `out/helmvec2/rect/csv`, com `helmvec2_rect_modes.csv`
  e `helmvec2_rect_candidates.csv`
- CSVs espaciais por modo casado, com `Et` por celula e `Ez` por no
- VTKs por modo casado em `out/helmvec2/rect/vtk`
- imagens geradas por `python3 scripts/helmvec2.py` em
  `out/helmvec2/rect/img`, incluindo resumo da Tabela 8, espectro de
  candidatos, magnitude/quiver de `Et` e mapa escalar de `Ez`

Observacao:

- a Eq. `(120)` impressa no artigo esta incoerente com as Eq. `(66)`, `(67)` e
  `(113)`: falta o fator `beta^2` no termo de massa vetorial e, na forma
  fatorada por `1/(16 A^3)`, tambem falta um fator `4` no termo `D_m D_n`.
- o codigo deste repositorio permanece correto porque nao copia a Eq. `(120)`
  isoladamente; o bloco `A_tt` e reconstruido por reaproveitamento dos blocos
  validados de `curl-curl` e massa vetorial.
- ver a nota detalhada em [`src/helmvec2/README.md`](src/helmvec2/README.md).

### 4.5) Secao 2.2.4 (`helmvec3`)

```bash
./build/helmvec3_fig12_rect 10 5
./build/helmvec3_fig12_rect 10 5 1
./build/helmvec3_fig12_rect 10 5 --backend closed-form
# args: nx ny [debug]

./build/helmvec3_fig13_rect 0.20 10 5
./build/helmvec3_fig13_rect 0.20 10 5 1
./build/helmvec3_fig13_rect 0.20 10 5 --backend closed-form
# args: d_over_a_preview nx ny [debug]
```

Saidas tipicas:

- `helmvec3_fig12_rect`
  - Tabela 9 (Figura 12)
  - `out/helmvec3/fig12_rect/run.log`
  - `out/helmvec3/fig12_rect/run_timing.csv`
  - `out/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_table9.csv`
  - CSVs espaciais por ponto exportado da Figura 12, com `Et` por celula e `Ez` por no
  - VTKs por ponto exportado em `out/helmvec3/fig12_rect/vtk`
- `helmvec3_fig13_rect`
  - preview de ramo para Figura 13
  - validacao da Tabela 10 (Figura 13)
  - `out/helmvec3/fig13_rect/run.log`
  - `out/helmvec3/fig13_rect/run_timing.csv`
  - `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_preview.csv`
  - `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv`
  - CSVs espaciais por ponto exportado, com `Et` por celula e `Ez` por no
  - VTKs por ponto exportado em `out/helmvec3/fig13_rect/vtk`

## 5) Executaveis 3D (Secao 3.1)

### 5.1) FEM3D0 (denso)

```bash
./build/fem3d0_air
./build/fem3d0_half
./build/fem3d0_cyl
./build/fem3d0_sphere
./build/fem3d0_air --nx 5 --ny 4 --nz 4
```

Saidas tipicas por caso (`air`, `half`, `cyl`, `sphere`):

- `out/fem3d0/run.log` com a trilha textual completa da execucao
- `out/fem3d0/<caso>/run.log`
- `out/fem3d0/<caso>/run_timing.csv`
- `out/fem3d0/<caso>/csv/fem3d0_<caso>_modes.csv`
- `out/fem3d0/<caso>/csv/fem3d0_<caso>_modeXX_<modo>_E_fields.csv`
- `out/fem3d0/<caso>/vtk/fem3d0_<caso>_modeXX_<modo>_E.vtk`
- `out/fem3d0/<caso>/img/` com resumos e projecoes `XY/XZ/YZ`
- `out/fem3d0/<caso>/linop/` com `S`, `T`, autovalores e autovetores em CSV

### 5.2) FEM3D1 (montagem esparsa)

```bash
./build/fem3d1_air
./build/fem3d1_half
./build/fem3d1_cyl
./build/fem3d1_sphere
./build/fem3d1_half --nx 6 --ny 4 --nz 4
```

Saidas tipicas por caso (`air`, `half`, `cyl`, `sphere`):

- `out/fem3d1/run.log` com a trilha textual completa da execucao
- `out/fem3d1/<caso>/run.log`
- `out/fem3d1/<caso>/run_timing.csv`
- `out/fem3d1/<caso>/csv/fem3d1_<caso>_modes.csv`
- `out/fem3d1/<caso>/csv/fem3d1_<caso>_modeXX_<modo>_E_fields.csv`
- `out/fem3d1/<caso>/vtk/fem3d1_<caso>_modeXX_<modo>_E.vtk`
- `out/fem3d1/<caso>/img/` com resumos e projecoes `XY/XZ/YZ`
- `out/fem3d1/<caso>/linop/` com `S`, `T`, autovalores e autovetores em CSV

## 6) Scripts de validacao

### 6.1) Validacao 2D (Secao 2.2)

```bash
python3 scripts/validate_2d_22.py
```

CLI util:

```bash
python3 scripts/validate_2d_22.py \
  --build-dir build \
  --backend closed-form \
  --out-csv out/validation/validation_2d_22.csv \
  --rect-nx 12 --rect-ny 6 \
  --circle-nr 10 --circle-nt 48 \
  --coax-nr 10 --coax-nt 48 \
  --beta 10 --hv2-nx 6 --hv2-ny 6 \
  --d-over-a 0.20 --hv3-nx 10 --hv3-ny 5 \
  --show-output \
  --debug-local-blocks \
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
  --backend closed-form \
  --out-modes out/validation/validation_3d_31_modes.csv \
  --out-summary out/validation/validation_3d_31_summary.csv \
  --show-output \
  --debug-candidates \
  --verbose
```

Saidas:

- `out/validation/validation_3d_31_modes.csv`
- `out/validation/validation_3d_31_summary.csv`

## 7) Plot de campos VTK (quiver)

Modo arquivo unico:

```bash
python3 scripts/plot_vtk_quiver.py out/helm10/rect/vtk/te10_rect_sv.vtk --out out/img/te10_rect.png --stride 2 --scale 22 --dpi 210
python3 scripts/plot_vtk_quiver.py out/helm10/rect/vtk/tm11_rect_sv.vtk --out out/img/tm11_rect.png --stride 2 --scale 22 --dpi 210
python3 scripts/plot_vtk_quiver.py out/helmvec/rect/vtk/edge_rect_Et.vtk --out out/img/edge_rect_Et.png --stride 2 --scale 25 --dpi 210
```

Modo lote para campos 2D-base (Sec. 2.1 e 2.2.1; gera imagens e CSV resumo):

```bash
python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out --out-dir out/img_all --csv out/img_all/mode_summary.csv --mode-export 8 --max-rank 8
```

O modo lote descobre VTKs no layout atual por familia:

- `out/helm10/*/vtk` para `2.1_scalar`
- `out/helmvec/*/vtk` para `2.2.1_edge`

Compatibilidade barata:

- se uma arvore legada `out/2d/` existir, ela continua aceita como entrada.

Imagens de validacao para os casos 2.2.2, 2.2.3 e 2.2.4:

```bash
python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22
```

Observacao:

- `--all_img` continua aceito como alias de compatibilidade.

Saidas de `plot_vtk_quiver.py --all-img`:

- imagens em `out/img_all/2.1_scalar/` e `out/img_all/2.2.1_edge/`
- `out/img_all/mode_summary.csv`

Saidas complementares de `plot_validation_2d_22.py`:

- `out/img_all/validation_2d_22/` (graficos de 2.2.2/2.2.3/2.2.4)

## 8) Fluxo recomendado de reproducao

1. Compilar (`cmake --build . -j`).
2. Rodar casos 2D-base (`helm10_*`, `edge_*`, `mixed_*`).
3. Rodar acoplados (`helmvec2_rect`, `helmvec3_fig12_rect`, `helmvec3_fig13_rect`).
4. Rodar validacao 2D (`validate_2d_22.py`).
5. Rodar validacao 3D (`validate_3d_31.py`).
6. Gerar figuras de campos 2D-base (`plot_vtk_quiver.py --all-img ...`) e, separadamente, as figuras de validacao (`plot_validation_2d_22.py`).

## 9) Script unico e comparacao entre backends

Pipeline principal:

```bash
./scripts/build_and_run_all.sh
./scripts/build_and_run_all.sh --backend closed-form
./scripts/build_and_run_all.sh --backend closed-form --debug-local-blocks
./scripts/build_and_run_all.sh --backend closed-form --show-validation-output
./scripts/build_and_run_all.sh --case 2.2.2 --backend closed-form --debug-candidates
```

Wrapper para comparar `gauss` e `closed-form` em diretorios separados:

```bash
./scripts/run_backend_compare.sh --backend-mode both -- --with-validate --with-images
./scripts/run_backend_compare.sh --interactive
```

## 10) Notas numericas

- Matching modal usa correlacao e tratamento de degenerescencia para reduzir troca artificial de ordem.
- Em problemas generalizados nao simetricos (`helmvec2`, `helmvec3`), o pipeline filtra raizes nao fisicas (parte imaginaria, sinal e faixa fisica).
- Na Sec. `2.2.3`, a Eq. `(120)` do artigo tem inconsistencias de impressao; a
  implementacao usa a decomposicao em blocos elementares, que preserva a forma
  correta da montagem local.
- `fem3d1` usa montagem esparsa simetrica; o solve atual ainda converte para denso antes de `dsygv`.

## 11) Guia de reproducao do paper (checklist)

Use esta sequencia para reproduzir os blocos numericos em ordem de tabelas/figuras.

1. Tabela 1 (retangular escalar, Sec. 2.1):

```bash
./build/helm10_rect 1.0 14 14 8
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
./build/helmvec3_fig12_rect 10 5
./build/helmvec3_fig13_rect 0.20 10 5
```

8. Secao 3.1 em cavidades 3D (Tabelas 12-15):

```bash
./build/fem3d0_air
./build/fem3d0_half
./build/fem3d0_cyl
./build/fem3d0_sphere
./build/fem3d1_air
./build/fem3d1_half
./build/fem3d1_cyl
./build/fem3d1_sphere
```

9. Validacao automatica consolidada:

```bash
python3 scripts/validate_2d_22.py --build-dir build --out-csv out/validation/validation_2d_22.csv
python3 scripts/validate_3d_31.py --profile quick --solver both --build-dir build --out-modes out/validation/validation_3d_31_modes.csv --out-summary out/validation/validation_3d_31_summary.csv
```

10. Geracao de imagens de campos 2D-base e CSV de modos:

```bash
python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out --out-dir out/img_all --csv out/img_all/mode_summary.csv --mode-export 8 --max-rank 8
```

11. Graficos de validacao 2.2.2 / 2.2.3 / 2.2.4:

```bash
python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22
```

## 12) Script unico (compila e roda tudo)

Foi adicionado:

- `scripts/build_and_run_all.sh`

Ele executa:

- configuracao e build CMake,
- todos os executaveis 2D,
- todos os executaveis 3D (`fem3d0` e `fem3d1`),
- validacoes 2D e 3D,
- geracao das imagens de campos 2D-base (`2.1` e `2.2.1`) e de `mode_summary.csv`,
- geracao dos graficos de validacao `2.2.2` / `2.2.3` / `2.2.4` quando `out/validation/validation_2d_22.csv` estiver disponivel,
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
