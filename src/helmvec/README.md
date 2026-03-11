# `helmvec` - Formulacao vetorial transversal 2D para `kc` (Sec. 2.2.1)

Este modulo implementa o problema vetorial transversal com elementos de aresta
(Whitney 1-form) para reduzir modos espurios classicos de discretizacoes nodais
vetoriais.

## 1) Objetivo do modulo

Resolver o cutoff `kc` em guias homogeneos usando apenas campo transversal:

- ramo TE (campo eletrico transversal),
- ramo dual TM (campo magnetico transversal).

Executaveis:
- `main_edge_rect.cpp`
- `main_edge_circle.cpp`
- `main_edge_coax.cpp`

Utilitario compartilhado:
- `edge_mode_post.hpp`

## 2) Espaco discreto de aresta

Em cada triangulo:

`Et_h(x,y) = sum_{m=1..3} e_m W_m(x,y)`

com Whitney local:

`W_m = L_m (lambda_i grad(lambda_j) - lambda_j grad(lambda_i))`

Propriedades usadas no codigo:
- DOF associado a circulacao tangencial na aresta,
- continuidade tangencial entre elementos,
- orientacao global de arestas com fator de sinal (`sgn`).

## 3) Sistema matricial

Problema generalizado simetrico:

`S e = lambda T e`, `lambda = kc^2`

com:

`S(m,n) = int (1/mu_r) curl_t(W_m) curl_t(W_n) dA`

`T(m,n) = int eps_r W_m . W_n dA`

Solver:
- `generalized_eigs_sym_vec` (`LAPACKE_dsygv`).

## 4) Condicoes de contorno

- `EdgeBC::TE_PEC_TangentialZero`
  - PEC para `Et`: zera componente tangencial no contorno,
  - elimina arestas de contorno do espaco de DOFs.

- `EdgeBC::TM_PEC_NormalZero`
  - formulacao dual para ramo TM,
  - mantem arestas de contorno (condicao natural no fraco).

## 5) Pos-processamento e VTK

`edge_mode_post.hpp` centraliza:
- impressao de primeiros `kc` positivos,
- escolha do primeiro modo fisico,
- reconstrucao do campo no centro de cada celula:

`Ft(xc) = sum_m e_m W_m(xc)`

- normalizacao para visualizacao.

Exportacao:
- `write_vtk_unstructured_tri_cell_vector(...)` em `src/core/io_vtk_sv.hpp`
- campo vetorial em `CELL_DATA`.

## 6) Arquivos e responsabilidades

- `src/edge/edge_dofs.{hpp,cpp}`: numeracao global de arestas e BC.
- `src/edge/edge_basis.{hpp,cpp}`: funcoes Whitney locais.
- `src/edge/edge_assembly.{hpp,cpp}`: montagem `S` e `T` vetoriais.
- `src/edge/mode_match_*_edge.hpp`: matching analitico/FEM por correlacao.

## 7) Uso

Retangular:

```bash
./build/edge_rect 14 14
```

Circular:

```bash
./build/edge_circle 10 48
```

Coaxial:

```bash
./build/edge_coax 10 48
```

## 8) Saidas tipicas

Console:
- primeiros `kc`,
- tabela FEM x analitico (`kc_ana`, `kc_fem`, erro %, `rho`).

VTK:
- retangular: `edge_rect_Et.vtk`, `edge_rect_Ht.vtk`
- circular: `edge_circle_Et.vtk`, `edge_circle_Ht.vtk`
- coaxial: `edge_coax_Et.vtk`, `edge_coax_Ht.vtk`

## 9) Relacao com o artigo

Corresponde a Sec. 2.2.1 (formulacao vetorial transversal).

Extensoes de engenharia no repositorio:
- matching modal automatico,
- reconstrucao de campo por celula pronta para quiver e figuras comparativas.
