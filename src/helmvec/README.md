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

## 3.1) Coeficientes locais e formas fechadas

Para a aresta local `m = (i,j)`, o codigo usa os coeficientes:

```text
A_m = a_i b_j - a_j b_i
B_m = c_i b_j - c_j b_i
C_m = a_i c_j - a_j c_i
D_m = b_i c_j - b_j c_i
```

e os momentos geometricos:

```text
x_tri = (x1 + x2 + x3) / 3
y_tri = (y1 + y2 + y3) / 3

(1/A) int_T x^2 dA = (sum x_i^2 + 9 x_tri^2) / 12
(1/A) int_T y^2 dA = (sum y_i^2 + 9 y_tri^2) / 12
```

Assim, a base de Whitney pode ser lida como:

```text
W_m(x,y) = (L_m / (4 A^2)) * [ (A_m + B_m y) x_hat + (C_m + D_m x) y_hat ]
```

e as matrizes elementares usadas pelo backend `closed-form` ficam:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(4 A^3) * D_m D_n
```

```text
T_e(m,n) = eps_r * (L_m L_n)/(16 A^3) * (It1 + It2 + It3 + It4 + It5)
```

com:

```text
It1 = A_m A_n + C_m C_n
It2 = (C_m D_n + C_n D_m) x_tri
It3 = (A_m B_n + A_n B_m) y_tri
It4 = B_m B_n * (1/A) int_T y^2 dA
It5 = D_m D_n * (1/A) int_T x^2 dA
```

Na forma global:

```text
S e = kc^2 T e
```

## 3.2) Trilha de rastreabilidade

Para conectar o artigo ao codigo, a trilha principal deste modulo e:

1. Eq. `(66)` e Eq. `(67)` -> formas locais closed-form em
   `src/explicit/tri2d_edge_explicit.hpp`
2. entrada publica didatica da Eq. `(65)` ->
   `tp3485::build_eq65_helmvec_system(...)` em `src/article/tp3485_systems.hpp`
3. Eq. `(65)` -> sistema global vetorial montado em
   `src/edge/edge_assembly.cpp`
4. funcao de montagem principal subjacente ->
   `build_helm10_edge_system(...)`

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
./build/edge_rect 14 14 20
```

Circular:

```bash
./build/edge_circle 10 48
./build/edge_circle 10 48 20
```

Coaxial:

```bash
./build/edge_coax 10 17
./build/edge_coax 10 17 10
```

Convencoes publicas atuais:

- `nmodos` indica quantos modos positivos de cada ramo serao exportados;
- se `nmodos` nao for informado, o padrao publico agora e `10`;
- sem `--backend`, o fluxo principal do executavel e `closed-form`.

## 7.1) Backends e depuracao

Flags disponiveis:

- `--backend closed-form`
- `--backend gauss`
- `--debug-local-blocks`
- `--debug-candidates`
- `--debug` ou `--debug-all`

Sem `--backend`, o padrao publico do executavel e `closed-form`.

Exemplos:

```bash
./build/edge_rect 14 14 8 --backend closed-form --debug-local-blocks
./build/edge_circle 10 48 8 --debug-candidates
./build/edge_coax 10 48 8 --backend closed-form --debug-local-blocks --debug-candidates
```

Interpretacao:

- `--debug-local-blocks`: imprime o primeiro triangulo com os blocos locais das
  Eq. `(66)` e `(67)` e sua relacao com a Eq. `(65)`;
- `--debug-candidates`: imprime os primeiros `kc` positivos antes do matching
  modal.

## 8) Saidas tipicas

Console:
- primeiros `kc`,
- tabela FEM x analitico (`kc_ana`, `kc_fem`, erro %, `rho`).

Arvore de saida:

- `out/helmvec/rect/`
- `out/helmvec/circle/`
- `out/helmvec/coax/`

Cada caso agora grava:

- `run.log`
  - copia completa do que foi impresso no terminal durante a execucao.
- `run_timing.csv`
  - configuracao da rodada e tempos de montagem, solve e pos-processamento.
- `csv/edge_<caso>_modes.csv`
  - uma linha por modo exportado.
- `csv/edge_<caso>_fields_<familia>_m<...>_<n/p><...>_rank<rr>.csv`
  - campo transversal por celula.
- `vtk/`
  - VTK por modo e aliases legados do primeiro modo.
- `img/`
  - imagens geradas por `python3 scripts/helmvec.py`.
- `linop/`
  - matrizes globais `S` e `T` em CSV no formato CRS.
  - `*_eigenvalues.csv` com os autovalores ordenados.
  - `*_eigenvectors.csv` com os autovetores na mesma ordem.
  - no `HELMVEC`, isso aparece separado por ramo `TE` e `TM`.

Cabecalho atual do `modes.csv`:

- retangular:
  - `family,transverse_label,mode_label,positive_rank,eig_index,m,n,ar_m,b_m,kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,error_percent,rho_abs,field_status,fields_csv_file,vtk_file`
- circular:
  - `family,transverse_label,mode_label,positive_rank,eig_index,m,p,r_m,kc_fem,kc_ana,kc_r_fem,kc_r_ana,error_percent,rho_abs,field_status,fields_csv_file,vtk_file`
- coaxial:
  - `family,transverse_label,mode_label,positive_rank,eig_index,m,p,r1_m,r2_m,kc_fem,kc_ana,kc_r1_fem,kc_r1_ana,error_percent,rho_abs,field_status,fields_csv_file,vtk_file`

Cabecalho atual do `fields_<modo>.csv`:

- no ramo `TE`:
  - `cell_id,xc_m,yc_m,Ex,Ey,Emag`
- no ramo `TM`:
  - `cell_id,xc_m,yc_m,Hx,Hy,Hmag`

Observacao didatica importante:

- o campo salvo em `fields_<modo>.csv` e reconstruido no centroide de cada
  triangulo;
- no ramo `TE`, esse campo e `E_t`, por isso o CSV usa `Ex,Ey,Emag`;
- no ramo `TM`, esse campo e `H_t`, por isso o CSV usa `Hx,Hy,Hmag`;
- no estado atual do projeto, esse campo ja sai normalizado por pico unitario
  para visualizacao;
- isso aparece no `modes.csv` pela coluna `field_status`.

## 8.1) Script de imagens

O script:

```bash
python3 scripts/helmvec.py
```

le os CSVs e gera, para cada modo:

- mapa de magnitude do campo transversal em `img/magnitude/`
- diagrama de setas em `img/quiver/`

e, por caso:

- `img/helmvec_<caso>_error_by_mode.png`
- `img/helmvec_<caso>_rho_by_mode.png`

Referencias didaticas:

- `docs/HELMVEC_CSV_Modos_Referencia.md`
- `docs/HELMVEC_CSV_Campos_Referencia.md`
- `docs/HELMVEC_Imagens_Referencia.md`
- `docs/Artefatos_Espectrais_CSV_Referencia.md`

## 9) Relacao com o artigo

Corresponde a Sec. 2.2.1 (formulacao vetorial transversal).

Extensoes de engenharia no repositorio:
- matching modal automatico,
- reconstrucao de campo por celula pronta para quiver e figuras comparativas.
