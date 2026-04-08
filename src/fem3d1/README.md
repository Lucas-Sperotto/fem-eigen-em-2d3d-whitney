# `fem3d1` - Solver vetorial 3D com montagem esparsa simetrica (Sec. 3.1)

Versao didatica: `2.0`

`fem3d1` usa a mesma formulacao de `fem3d0`, mas com acumulacao em
armazenamento esparso simetrico (`SparseSymMat`).

## 1) Objetivo do modulo

Reduzir custo de montagem/memoria mantendo compatibilidade com a validacao
3D do artigo.

Sistema discreto:

`S e = k0^2 T e`

Diferenca principal:

- `S` e `T` sao montadas no triangulo inferior esparso,
- solve atual ainda converte para denso antes de `dsygv`.

Entradas publicas:

- `main_fem3d1_air.cpp`
- `main_fem3d1_half.cpp`
- `main_fem3d1_cyl.cpp`
- `main_fem3d1_sphere.cpp`

Implementacao compartilhada:

- `main_fem3d1_rect.cpp`

## 2) Pipeline numerico

1. Montagem esparsa:
   - `tp3485::build_eq178_fem3d_system_sparse(...)`
   - chamada interna: `build_helm3d_edge_system_sparse(...)`
   - inicializacao: `initialize_eq178_sparse_global_system(...)`
   - fechamento: `assemble_eq178_global_sparse(...)`
2. Conversao para denso:
   - `S = sparse.S.to_dense()`,
   - `T = sparse.T.to_dense()`.
3. Solve generalizado simetrico:
   - `generalized_eigs_sym_vec(S, T)`.
4. Pos-processamento identico ao `fem3d0`:
   - primeiras raizes positivas,
   - matching com degenerescencia,
   - tabela comparativa,
   - exportacao vetorial 3D por modo casado.

## 2.1) Expressoes matematicas usadas

Problema continuo:

```text
curl((1/mu_r) curl(E)) = k0^2 eps_r E
```

Problema discreto:

```text
S e = k0^2 T e
```

Mesmo com montagem esparsa, o elemento local do backend `closed-form` e o
mesmo do `fem3d0`:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(324 V^3) * K_mn
T_e(m,n) = eps_r * (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

com a mesma base simplex:

```text
lambda_i(x,y,z) = (a_i + b_i x + c_i y + d_i z) / (6V)
grad lambda_i   = [b_i, c_i, d_i] / (6V)
W_m             = L_m (lambda_i grad lambda_j - lambda_j grad lambda_i)
```

Na linguagem das equacoes do artigo:

- Eq. `(181)` fornece o bloco local de rigidez `Sel`
- Eq. `(182)` fornece o bloco local de massa `Tel`
- Eq. `(178)` no tetraedro -> `build_eq178_local_tet_blocks(...)`
- Eq. `(178)` na montagem compartilhada -> `assemble_eq178_global_generic(...)`
- Eq. `(178)` no fechamento esparso -> `assemble_eq178_global_sparse(...)`
- Eq. `(178)` e o problema global montado em formato esparso
- entrada pública didática -> `tp3485::build_eq178_fem3d_system_sparse(...)`

A diferenca deste modulo esta apenas na estrutura de armazenamento global
(`SparseSymMat`), nao na formulacao matematica local.

## 3) Casos suportados

- `fem3d1_air` (Figura 15 / Tabela 12)
- `fem3d1_half` (Figura 16 / Tabela 13)
- `fem3d1_cyl` (Figura 17 / Tabela 14)
- `fem3d1_sphere` (Tabela 15)

Cada executavel ja nasce fixo em um caso do artigo:

- sem flags de caso;
- com `--nx --ny --nz`: sobrescreve a malha padrao daquele caso;
- sem `--backend`: usa `closed-form` como fluxo publico padrao.

## 4) Uso

```bash
./build/fem3d1_air
./build/fem3d1_half
./build/fem3d1_cyl
./build/fem3d1_sphere
./build/fem3d1_half --nx 6 --ny 4 --nz 4
./build/fem3d1_half --backend closed-form --debug-local-blocks
./build/fem3d1_half --debug-candidates
./build/fem3d1_half --help
```

## 5) Saida tipica

Para cada caso:

- `nodes`, `tets`, `edges`, `dof`,
- `nnz_lower(S)`, `nnz_lower(T)`,
- tabela modal com:
  - `k0_ana`, `k0_fem`, `err_ana(%)`,
  - `ref_paper`, `err_ref(%)`.

Flags de depuracao:

- `--debug-local-blocks`: imprime `Sel` e `Tel` do primeiro tetraedro
- `--debug-candidates`: imprime as primeiras raizes positivas `k0`
- `--debug` / `--debug-all`: ativa os dois

## 5.1) Artefatos salvos por caso

Cada caso agora grava:

- `out/fem3d1/run.log`
  - trilha textual completa da execucao, cobrindo todos os casos selecionados.
- `out/fem3d1/<caso>/run.log`
  - trilha textual especifica do caso.
- `out/fem3d1/<caso>/run_timing.csv`
  - configuracao da malha, estrutura esparsa e tempos de montagem, solve e pos-processamento.
- `out/fem3d1/<caso>/csv/fem3d1_<caso>_modes.csv`
  - resumo modal alinhado com as Tabelas 12-15.
- `out/fem3d1/<caso>/csv/fem3d1_<caso>_modeXX_<modo>_E_fields.csv`
  - campo vetorial reconstruido no centroide de cada tetraedro.
- `out/fem3d1/<caso>/vtk/*.vtk`
  - malha tetraedrica com `CELL_DATA` para `E` e `Emag`.
- `out/fem3d1/<caso>/img/`
  - resumos modais, projecoes ortogonais e vistas 3D geradas por `python3 scripts/fem3d.py`.
- `out/fem3d1/<caso>/linop/`
  - `S`, `T`, autovalores e autovetores em CSV.

Exemplo para `fem3d1_half`:

```text
out/fem3d1/half/run.log
out/fem3d1/half/run_timing.csv
out/fem3d1/half/csv/fem3d1_half_modes.csv
out/fem3d1/half/csv/fem3d1_half_mode01_TE101_E_fields.csv
out/fem3d1/half/vtk/fem3d1_half_mode01_TE101_E.vtk
out/fem3d1/half/img/fem3d1_half_k0_by_mode.png
out/fem3d1/half/img/quiver/fem3d1_half_mode01_TE101_E_quiver_proj.png
out/fem3d1/half/img/3d_quiver/fem3d1_half_mode01_TE101_quiver3d.png
out/fem3d1/half/linop/fem3d1_half_S_crs.csv
out/fem3d1/half/linop/fem3d1_half_T_crs.csv
out/fem3d1/half/linop/fem3d1_half_eigenvalues.csv
out/fem3d1/half/linop/fem3d1_half_eigenvectors.csv
```

## 6) Dependencias internas

- `src/core/sparse_sym.hpp`
- `src/edge3d/edge3d_assembly.{hpp,cpp}`
- `src/edge3d/edge3d_basis.{hpp,cpp}`
- `src/edge3d/edge3d_dofs.{hpp,cpp}`
- `src/fem3d/fem3d_case_driver.hpp`
- `src/fem3d/fem3d_compare.hpp`
- `src/fem3d/fem3d_reference_tables.hpp`

## 7) Observacoes de engenharia

- Ja prepara o caminho para eigensolver iterativo esparso futuro
  (ARPACK/SLEPc/Lanczos), mantendo a camada de casos e comparacao estavel.
- Hoje, a etapa de solve ainda compartilha o caminho LAPACK denso para
  preservar reproducibilidade dos resultados atuais.

## 8) Relacao com `fem3d0`

- Formulacao fisica e matching: iguais.
- Diferenca: estrutura de dados da montagem (`dense` vs `sparse`).
- Ambos sao validados pelo mesmo script:
  - `scripts/validate_3d_31.py`.
- Ambos tambem compartilham a mesma camada grafica 3D:
  - `scripts/fem3d.py`.
- Ambos tambem compartilham a mesma camada grafica 3D:
  - `scripts/fem3d.py`.
