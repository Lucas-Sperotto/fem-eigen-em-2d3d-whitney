# `fem3d0` - Solver vetorial 3D denso para cavidades (Sec. 3.1)

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
   - `build_helm3d_edge_system(...)`
2. Resolve EVP simetrico:
   - `generalized_eigs_sym_vec(S, T)`
3. Extrai primeiras raizes positivas `k0`.
4. Faz matching por degenerescencia analitica.
5. Imprime tabela comparativa.

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

## 4) Uso

```bash
./build/fem3d0_rect
./build/fem3d0_rect --half
./build/fem3d0_rect --cyl
./build/fem3d0_rect --sphere
./build/fem3d0_rect --all
./build/fem3d0_rect --all --nx 5 --ny 4 --nz 4
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
