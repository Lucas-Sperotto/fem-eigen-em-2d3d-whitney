# `fem3d1` - Solver vetorial 3D com montagem esparsa simetrica (Sec. 3.1)

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

Arquivo principal:
- `main_fem3d1_rect.cpp`

## 2) Pipeline numerico

1. Montagem esparsa:
   - `build_helm3d_edge_system_sparse(...)`
2. Conversao para denso:
   - `S = sparse.S.to_dense()`,
   - `T = sparse.T.to_dense()`.
3. Solve generalizado simetrico:
   - `generalized_eigs_sym_vec(S, T)`.
4. Pos-processamento identico ao `fem3d0`:
   - primeiras raizes positivas,
   - matching com degenerescencia,
   - tabela comparativa.

## 3) Casos suportados

- `--air` (Figura 15 / Tabela 12)
- `--half` (Figura 16 / Tabela 13)
- `--cyl` (Figura 17 / Tabela 14)
- `--sphere` (Tabela 15)
- `--all`

Defaults do executavel:
- sem flags: roda `--air` e `--half`.
- com `--nx --ny --nz`: sobrescreve malha padrao do caso selecionado.

## 4) Uso

```bash
./build/fem3d1_rect --air
./build/fem3d1_rect --half
./build/fem3d1_rect --cyl
./build/fem3d1_rect --sphere
./build/fem3d1_rect --all
./build/fem3d1_rect --all --nx 6 --ny 4 --nz 4
./build/fem3d1_rect --help
```

## 5) Saida tipica

Para cada caso:
- `nodes`, `tets`, `edges`, `dof`,
- `nnz_lower(S)`, `nnz_lower(T)`,
- tabela modal com:
  - `k0_ana`, `k0_fem`, `err_ana(%)`,
  - `ref_paper`, `err_ref(%)`.

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
