# `fem3d0` - Solver vetorial 3D denso para cavidades (Sec. 3.1)

Versao didatica: `2.0`

`fem3d0` e a implementacao baseline 3D:
- elementos de aresta tetraedricos,
- montagem global densa,
- solve generalizado simetrico com LAPACK.

## 1) Objetivo do modulo

Resolver autovalores de cavidades PEC em 3D:

`curl((1/mu_r)curl(E)) = k0^2 eps_r E`

na forma discreta:

`S e = k0^2 T e`

com comparacao contra:
- analitico,
- referencia publicada (FEM3D1 / ref. 17 do artigo).

Arquivo principal:
- `main_fem3d0_rect.cpp`

## 2) Pipeline numerico

1. Monta sistema edge 3D denso:
   - `tp3485::build_eq178_fem3d_system_dense(...)`
   - chamada interna: `build_helm3d_edge_system(...)`
   - inicializacao: `initialize_eq178_dense_global_system(...)`
   - fechamento: `assemble_eq178_global_dense(...)`
2. Resolve EVP simetrico:
   - `generalized_eigs_sym_vec(S, T)`
3. Extrai primeiras raizes positivas `k0`.
4. Faz matching por degenerescencia analitica.
5. Imprime tabela comparativa.

## 2.1) Expressoes matematicas usadas

Problema continuo:

```text
curl((1/mu_r) curl(E)) = k0^2 eps_r E
```

Problema discreto:

```text
S e = k0^2 T e
```

No backend `closed-form`, os blocos locais do tetraedro seguem:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(324 V^3) * K_mn
T_e(m,n) = eps_r * (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

com:

```text
lambda_i(x,y,z) = (a_i + b_i x + c_i y + d_i z) / (6V)
grad lambda_i   = [b_i, c_i, d_i] / (6V)
W_m             = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
```

Na rastreabilidade do artigo:

- Eq. `(181)` -> bloco local `Sel`
- Eq. `(182)` -> bloco local `Tel`
- Eq. `(178)` no tetraedro -> `build_eq178_local_tet_blocks(...)`
- Eq. `(178)` na montagem compartilhada -> `assemble_eq178_global_generic(...)`
- Eq. `(178)` no fechamento denso -> `assemble_eq178_global_dense(...)`
- Eq. `(178)` -> sistema global `S e = k0^2 T e`
- entrada pública didática -> `tp3485::build_eq178_fem3d_system_dense(...)`

Os detalhes de `K_mn` e `I1..I10` ficam documentados em:

- `src/fem3d/README.md`
- `src/explicit/tet3d_edge_explicit.hpp`

## 3) Casos suportados

- `--air` (Figura 15 / Tabela 12)
  - cavidade retangular com ar.
- `--half` (Figura 16 / Tabela 13)
  - cavidade retangular meio preenchida em `z`.
- `--cyl` (Figura 17 / Tabela 14)
  - cavidade cilindrica com ar.
- `--sphere` (Tabela 15)
  - cavidade esferica com ar.
- `--all`
  - executa os quatro casos.

Defaults do executavel:
- sem flags: roda apenas `--air`.
- com `--nx --ny --nz`: sobrescreve malha padrao do caso selecionado.
- sem `--backend`: usa `closed-form` como fluxo publico padrao.

## 4) Uso

```bash
./build/fem3d0_rect
./build/fem3d0_rect --half
./build/fem3d0_rect --cyl
./build/fem3d0_rect --sphere
./build/fem3d0_rect --all
./build/fem3d0_rect --all --nx 5 --ny 4 --nz 4
./build/fem3d0_rect --air --backend closed-form --debug-local-blocks
./build/fem3d0_rect --air --debug-candidates
./build/fem3d0_rect --help
```

## 5) Saida tipica

Para cada caso:
- `nodes`, `tets`, `edges`, `dof`,
- tabela por modo:
  - `k0_ana`,
  - `k0_fem`,
  - `err_ana(%)`,
  - `ref_paper`,
  - `err_ref(%)`.

A tabela usa matching agrupado por degenerescencia para evitar trocas
artificiais de ordem modal.

Flags de depuracao:
- `--debug-local-blocks`: imprime `Sel` e `Tel` do primeiro tetraedro
- `--debug-candidates`: imprime as primeiras raizes positivas `k0`
- `--debug` / `--debug-all`: ativa os dois

## 6) Dependencias internas

- `src/edge3d/edge3d_assembly.{hpp,cpp}`
- `src/edge3d/edge3d_basis.{hpp,cpp}`
- `src/edge3d/edge3d_dofs.{hpp,cpp}`
- `src/fem3d/fem3d_case_driver.hpp`
- `src/fem3d/fem3d_compare.hpp`
- `src/fem3d/fem3d_reference_tables.hpp`

## 7) Custo computacional

Por ser denso:
- memoria cresce aproximadamente com `O(n^2)`,
- solve cresce aproximadamente com `O(n^3)`.

Esse modulo e ideal para:
- baseline de corretude,
- malhas moderadas,
- comparacao rapida com tabelas.

Para montagem esparsa, usar `fem3d1`.
