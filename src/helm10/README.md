# `helm10` - Formulacao escalar 2D para `kc` (Sec. 2.1)

Este modulo implementa a etapa escalar do artigo para guias homogeneos,
usando elementos triangulares nodais `P1` e solve generalizado simetrico.

## 1) Objetivo do modulo

Resolver:

- problema TE escalar (Neumann) para cutoff `kc`,
- problema TM escalar (Dirichlet) para cutoff `kc`,
- comparacao modal com referencias analiticas de geometria padrao.

Executaveis:
- `main_helm10_rect.cpp`
- `main_helm10_circle.cpp`
- `main_helm10_coax.cpp`

Utilitario compartilhado:
- `scalar_mode_post.hpp`

## 2) Modelo continuo e forma fraca

No corte transversal `Gamma`:

`nabla_t^2 phi + kc^2 phi = 0`

Forma fraca:

`int_Gamma (nabla_t T_s . nabla_t phi) dA = kc^2 int_Gamma (T_s phi) dA + int_dGamma T_s (dphi/dn) dl`

No codigo, isso vira o EVP:

`S u = lambda T u`, com `lambda = kc^2`.

## 3) Condicoes de contorno

- `ScalarBC::TE_Neumann`
  - condicao natural `dphi/dn = 0`,
  - mantem nos de contorno,
  - possui modo constante (`lambda ~ 0`) removido no pos-processamento.

- `ScalarBC::TM_Dirichlet`
  - condicao essencial `phi = 0` em PEC,
  - elimina nos de contorno do espaco discreto.

## 4) Discretizacao FEM (`P1`)

Por triangulo:

`phi_h = sum_{i=1..3} phi_i N_i`

com gradiente constante:

`grad N_i = [b_i, c_i] / (2A)`

Matrizes elementares:

`S_e(i,j) = int_e (1/mu_r) grad N_i . grad N_j dA`

`T_e(i,j) = int_e eps_r N_i N_j dA`

## 4.1) Trilha de rastreabilidade

Para conectar o artigo ao codigo, a trilha principal deste modulo e:

1. Eq. `(30)` e Eq. `(33)` -> formas locais closed-form em
   `src/explicit/tri2d_scalar_explicit.hpp`
2. Eq. `(43)` -> sistema global escalar montado em
   `src/core/helm10_scalar_system.cpp`
3. funcao de montagem principal ->
   `build_helm10_scalar_system(...)`

## 5) Solver e pos-processamento

Solver:
- `generalized_eigs_sym_vec` (`LAPACKE_dsygv`).

Pos-processamento (`scalar_mode_post.hpp`):
- filtro de autovalores fisicos,
- extracao de autovetor nodal,
- reconstrucao de campo transversal para visualizacao:

`Ft = z_hat x grad(phi) = (-dphi/dy, dphi/dx)`

- suavizacao nodal por media ponderada de area,
- normalizacao para quiver.

## 6) Arquivos e responsabilidades

- `src/core/helm10_scalar_system.{hpp,cpp}`: montagem escalar global.
- `src/core/fem_scalar_helm10.cpp`: blocos elementares e assembleia.
- `src/core/mode_match_rect.hpp`: matching retangular por correlacao.
- `src/core/mode_match_circle.hpp`: matching circular por correlacao.
- `src/core/mode_match_coax.hpp`: matching coaxial por correlacao.
- `src/core/io_vtk_sv.hpp`: export VTK escalar+vetorial nodal.

## 7) Uso

Retangular:

```bash
./build/helm10_rect 14 14 --backend gauss
```

Circular:

```bash
./build/helm10_circle 10 48 --backend gauss
```

Coaxial:

```bash
./build/helm10_coax 10 48 --backend gauss
```

## 7.1) Backends e depuracao

Flags disponiveis:

- `--backend gauss`
- `--backend closed-form`
- `--debug-local-blocks`
- `--debug-candidates`
- `--debug` ou `--debug-all`

Exemplos:

```bash
./build/helm10_rect 14 14 8 --backend closed-form --debug-local-blocks
./build/helm10_circle 10 48 8 --backend gauss --debug-candidates
./build/helm10_coax 10 48 8 --backend closed-form --debug-local-blocks --debug-candidates
```

Interpretacao:

- `--debug-local-blocks`: imprime o primeiro triangulo com os blocos locais
  ligados as Eq. `(30)` e `(33)` e sua contribuicao para a Eq. `(43)`;
- `--debug-candidates`: imprime os primeiros `kc` positivos antes do matching.

## 8) Saidas tipicas

Console:
- primeiros `kc`,
- tabela FEM x analitico (`kc_ana`, `kc_fem`, erro %, `rho`).

Arquivos VTK gerados pelos `main`:
- retangular: `te10_rect_sv.vtk`, `tm11_rect_sv.vtk`
- circular: `te_circle_sv.vtk`, `tm_circle_sv.vtk`
- coaxial: `te_coax_sv.vtk`, `tm_coax_sv.vtk`

## 9) Relacao com o artigo

Corresponde ao fluxo da Sec. 2.1:
- formulacao escalar,
- discretizacao triangular nodal,
- validacao em guias retangular/circular/coaxial.

Extensoes de engenharia no repositorio:
- matching modal automatico por correlacao em massa,
- exportacao VTK pronta para pipeline de figuras.
