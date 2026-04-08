# `fem3d` - Infra compartilhada para casos 3D (Sec. 3.1)

Versao didatica: `2.0`

Este diretorio centraliza componentes comuns usados por:
- `src/fem3d0` (montagem densa),
- `src/fem3d1` (montagem esparsa simetrica).

## 1) Objetivo do modulo

Evitar duplicacao de:
- definicao de casos e geometrias,
- tabelas de referencia,
- parsing/controle de CLI para selecao de casos,
- comparacao modal FEM x analitico x referencia do paper.

## 2) Arquivos e responsabilidades

- `fem3d_reference_tables.hpp`
  - geometrias de referencia:
    - Figura 15 (`1 x 0.5 x 0.75 cm`),
    - Figura 16 (`1 x 1 x 1 cm`),
    - Figura 17 (cilindro `diam=1`, `altura=0.5`),
    - Tabela 15 (esfera `raio=1`).
  - malhas padrao (`Grid3D`) por caso,
  - tabelas de referencia (12 a 15),
  - regra de `scan_limit_for_table`.

- `fem3d_compare.hpp`
  - extracao de primeiras raizes positivas (`k0`),
  - matching com agrupamento de degenerescencia analitica,
  - impressao de tabela comparativa estavel para parser.

- `fem3d_case_driver.hpp`
  - parser de CLI comum (`--air`, `--half`, `--cyl`, `--sphere`, `--all`, `--nx --ny --nz`),
  - parser dedicado por caso (`parse_single_case_cli(...)`),
  - selecao de casos,
  - construcao de malha/material por caso (`PreparedCase`),
  - fluxo `for_each_selected_case(...)`.
- `fem3d_case_output.hpp`
  - padronizacao de pastas `out/fem3d{0,1}/<caso>/`,
  - escrita de `modes.csv`,
  - apoio a `run.log` e `linop/`.
- `fem3d_field_output.hpp`
  - reconstrucao vetorial 3D por tetraedro,
  - escrita de `fields.csv` e VTK por modo casado.

## 3) Estruturas principais

- `Grid3D {nx, ny, nz}`
- `RefRow {mode, analytical, ref_paper}`
- `PreparedCase`
  - `id`, `case_name`, `header`, `grid`, `mesh`,
  - `eps_r_tet`, `mu_r_tet`,
  - `rows` (referencias modais).

## 3.1) Forma matematica compartilhada em 3D

A formulacao vetorial 3D do repositorio parte de:

```text
curl((1/mu_r) curl(E)) = k0^2 eps_r E
```

e, apos discretizacao por elementos de aresta tetraedricos:

```text
S e = k0^2 T e
```

No tetraedro:

```text
lambda_i(x,y,z) = (a_i + b_i x + c_i y + d_i z) / (6V)
grad lambda_i   = [b_i, c_i, d_i] / (6V)
W_m = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
```

O backend `closed-form` usa:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(324 V^3) * K_mn
T_e(m,n) = eps_r * (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

onde `K_mn` e os `I_k` sao construidos a partir dos coeficientes geometricos
das Eq. `(162)` a `(172)`.

Na implementacao atual:

```text
K_mn = C_zm C_zn + C_xm C_xn + B_ym B_yn
```

e os termos de massa seguem:

```text
I1  = A_xm A_xn + A_ym A_yn + A_zm A_zn
I2  = (A_ym B_yn + A_yn B_ym + A_zm B_zn + A_zn B_zm) x_tet
I3  = (A_xm B_xn + A_xn B_xm + A_zm C_zn + A_zn C_zm) y_tet
I4  = (A_xm C_xn + A_xn C_xm + A_ym C_yn + A_yn C_ym) z_tet
I5  = (B_zm C_zn + B_zn C_zm) * (1/V) int_T xy dV
I6  = (B_xm C_xn + B_xn C_xm) * (1/V) int_T yz dV
I7  = (B_ym C_yn + B_yn C_ym) * (1/V) int_T xz dV
I8  = (B_ym B_yn + B_zm B_zn) * (1/V) int_T x^2 dV
I9  = (B_xm B_xn + C_zm C_zn) * (1/V) int_T y^2 dV
I10 = (C_xm C_xn + C_ym C_yn) * (1/V) int_T z^2 dV
```

Os momentos geometricos usados nesses termos sao os equivalentes tetraedricos
dos momentos do triangulo 2D:

```text
x_tet = (x1 + x2 + x3 + x4) / 4
y_tet = (y1 + y2 + y3 + y4) / 4
z_tet = (z1 + z2 + z3 + z4) / 4

(1/V) int_T x^2 dV = (sum x_i^2 + 16 x_tet^2) / 20
(1/V) int_T y^2 dV = (sum y_i^2 + 16 y_tet^2) / 20
(1/V) int_T z^2 dV = (sum z_i^2 + 16 z_tet^2) / 20
(1/V) int_T xy dV  = (sum x_i y_i + 16 x_tet y_tet) / 20
(1/V) int_T xz dV  = (sum x_i z_i + 16 x_tet z_tet) / 20
(1/V) int_T yz dV  = (sum y_i z_i + 16 y_tet z_tet) / 20
```

Essas expressoes sao as que aparecem materializadas no backend
`closed-form` em `src/explicit/tet3d_edge_explicit.hpp`.

## 4) Beneficios para o projeto

- `main_fem3d0_rect.cpp` e `main_fem3d1_rect.cpp` concentram o nucleo compartilhado do solver.
- os wrappers `main_fem3d{0,1}_{air,half,cyl,sphere}.cpp` deixam explicito
  qual executavel reproduz cada caso do artigo.
- Mudancas de caso/parametro de referencia ficam em um unico lugar.
- Scripts de validacao (`validate_3d_31.py`) continuam com formato de saida
  estavel mesmo apos refatoracoes internas.
- Cada caso 3D agora gera um pacote didatico rastreavel com:
  - `run.log`,
  - `run_timing.csv`,
  - `csv/<solver>_<caso>_modes.csv`,
  - `csv/<solver>_<caso>_modeXX_<modo>_E_fields.csv`,
  - `vtk/<solver>_<caso>_modeXX_<modo>_E.vtk`,
  - `img/` com resumos modais e projecoes `XY/XZ/YZ`,
  - `linop/` com `S`, `T`, autovalores e autovetores.

## 4.1) Flags de depuracao compartilhadas

Os modulos internos compartilhados `main_fem3d0_rect.cpp` e
`main_fem3d1_rect.cpp` aceitam:

- `--debug-local-blocks`
  - imprime os blocos locais `Sel` e `Tel` do primeiro tetraedro,
    correspondentes as Eq. `(181)` e `(182)`
- `--debug-candidates`
  - imprime as primeiras raizes positivas `k0` antes do matching com a tabela
- `--debug` ou `--debug-all`
  - ativa os dois comportamentos

Essas flags tambem sao repassadas pelos scripts principais:

- `scripts/build_and_run_all.sh`
- `scripts/run_backend_compare.sh`
- `scripts/validate_3d_31.py`

Para acompanhar a saida completa durante a validacao automatica:

```bash
python3 scripts/validate_3d_31.py --backend closed-form --show-output --debug-candidates
```

Na casca publica atual, a recomendacao e usar os executaveis por caso:

- `fem3d0_air`, `fem3d0_half`, `fem3d0_cyl`, `fem3d0_sphere`
- `fem3d1_air`, `fem3d1_half`, `fem3d1_cyl`, `fem3d1_sphere`

## 5) Relacao com a secao 3.1 do artigo

Este modulo nao monta matrizes nem resolve EVP diretamente.
Ele organiza a camada de reproducao dos experimentos numericos da Sec. 3.1,
onde as montagens reais estao em `src/edge3d` e os solves em `fem3d0/fem3d1`.

## 5.1) Trilha de rastreabilidade

Para conectar o artigo ao codigo, a trilha principal dos casos 3D e:

1. Eq. `(181)` e Eq. `(182)` -> formas locais closed-form em
   `src/explicit/tet3d_edge_explicit.hpp`
2. Eq. `(178)` no nivel do tetraedro ->
   `build_eq178_local_tet_blocks(...)` em
   `src/edge3d/edge3d_assembly.cpp`
3. Eq. `(178)` na montagem global compartilhada ->
   `assemble_eq178_global_generic(...)` em
   `src/edge3d/edge3d_assembly.cpp`
4. Eq. `(178)` no fechamento por politica de armazenamento ->
   `assemble_eq178_global_dense(...)` e
   `assemble_eq178_global_sparse(...)` em
   `src/edge3d/edge3d_assembly.cpp`
5. entradas publicas didaticas da Eq. `(178)` ->
   `tp3485::build_eq178_fem3d_system_dense(...)` e
   `tp3485::build_eq178_fem3d_system_sparse(...)` em
   `src/article/tp3485_systems.hpp`
6. Eq. `(178)` -> sistema global vetorial 3D montado em
   `src/edge3d/edge3d_assembly.cpp`
7. funcoes de montagem principais subjacentes ->
   `build_helm3d_edge_system(...)` e
   `build_helm3d_edge_system_sparse(...)`

Os coeficientes locais usados pelo backend `closed-form` aparecem em:

```text
A_x, B_x, C_x, A_y, B_y, C_y, A_z, B_z, C_z
```

e alimentam os termos auxiliares:

```text
I1, I2, ..., I10
```

descritos em `src/explicit/tet3d_edge_explicit.hpp`.
