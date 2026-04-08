# FEM3D - Referência dos CSVs e Artefatos de Saída

Este arquivo documenta a camada didática de saída dos executáveis:

- `fem3d0_air`, `fem3d0_half`, `fem3d0_cyl`, `fem3d0_sphere`
- `fem3d1_air`, `fem3d1_half`, `fem3d1_cyl`, `fem3d1_sphere`

Ela cobre a Sec. 3.1 do artigo e organiza, por caso 3D, os seguintes artefatos:

- `run.log`
- `run_timing.csv`
- `csv/<solver>_<caso>_modes.csv`
- `csv/<solver>_<caso>_modeXX_<modo>_E_fields.csv`
- `vtk/<solver>_<caso>_modeXX_<modo>_E.vtk`
- `linop/*.csv`

## 1) Organização em disco

Cada solver mantém uma pasta-raiz:

```text
out/fem3d0/
out/fem3d1/
```

Quando um caso é executado, os arquivos ficam em:

```text
out/<solver>/<caso>/
```

com:

```text
run.log
run_timing.csv
csv/
vtk/
linop/
```

Exemplos:

- `out/fem3d0/air/run.log`
- `out/fem3d0/air/csv/fem3d0_air_modes.csv`
- `out/fem3d0/air/csv/fem3d0_air_mode01_TE101_E_fields.csv`
- `out/fem3d0/air/vtk/fem3d0_air_mode01_TE101_E.vtk`
- `out/fem3d0/air/linop/fem3d0_air_S_crs.csv`
- `out/fem3d1/sphere/run_timing.csv`

Observação importante:

- `out/fem3d0/run.log` e `out/fem3d1/run.log` registram a execução completa do solver;
- `out/fem3d0/<caso>/run.log` e `out/fem3d1/<caso>/run.log` registram apenas o caso correspondente;
- cada executável público já nasce fixo em um único caso, então não há mais
  necessidade de flags `--air`, `--half`, `--cyl`, `--sphere` ou `--all`.

## 2) `run_timing.csv`

Cabecalho:

```text
solver,case_name,storage,backend,debug_local_blocks,debug_candidates,
nx,ny,nz,mesh_nodes,mesh_tets,mesh_edges,ndof,nnz_lower_s,nnz_lower_t,
assembly_ms,solve_ms,post_ms,total_ms
```

Significado:

- `solver`
  - `fem3d0` ou `fem3d1`.
- `case_name`
  - `air`, `half`, `cyl` ou `sphere`.
- `storage`
  - `dense` no `FEM3D0`;
  - `sparse-symmetric` no `FEM3D1`.
- `backend`
  - `closed-form` ou `gauss`.
- `debug_local_blocks`, `debug_candidates`
  - flags de depuração efetivamente usadas na rodada.
- `nx`, `ny`, `nz`
  - discretização estruturada usada para gerar a malha do caso.
- `mesh_nodes`, `mesh_tets`, `mesh_edges`
  - tamanhos geométricos da malha final.
- `ndof`
  - número de graus de liberdade de aresta após aplicar a condição PEC.
- `nnz_lower_s`, `nnz_lower_t`
  - apenas para `FEM3D1`: número de não nulos armazenados no triângulo inferior.
  - em `FEM3D0`, esses campos ficam `0`.
- `assembly_ms`
  - tempo da montagem global da Eq. `(178)`.
- `solve_ms`
  - tempo do solve generalizado de autovalores.
- `post_ms`
  - tempo de pós-processamento: extração de `k0`, matching, escrita de CSVs e `linop/`.
- `total_ms`
  - tempo total do caso.

## 3) `modes.csv`

Cada caso gera:

- `fem3d0_air_modes.csv`
- `fem3d0_half_modes.csv`
- `fem3d0_cyl_modes.csv`
- `fem3d0_sphere_modes.csv`
- `fem3d1_air_modes.csv`
- etc.

Cabecalho:

```text
reference_index,mode_label,k0_analytic,k0_fem,ref_paper,
error_percent_analytic,error_percent_ref_paper,matched_eig_index,
match_status,field_status,fields_csv_file,vtk_file
```

Significado:

- `reference_index`
  - posição da linha na tabela de referência do artigo.
- `mode_label`
  - rótulo modal da tabela (`TE101`, `TM110`, `TE2201`, etc.).
- `k0_analytic`
  - valor analítico em forma fechada.
- `k0_fem`
  - valor FEM casado com essa linha de referência.
- `ref_paper`
  - valor publicado na implementação de referência citada pelo artigo.
- `error_percent_analytic`
  - erro relativo absoluto entre `k0_fem` e `k0_analytic`.
- `error_percent_ref_paper`
  - erro relativo absoluto entre `k0_fem` e `ref_paper`.
- `matched_eig_index`
  - indice do autovalor FEM efetivamente escolhido no espectro ordenado.
- `match_status`
  - `matched` quando o casamento foi realizado;
  - `missing` quando não houve raiz positiva disponível para aquela linha.
- `field_status`
  - resume a política usada para exportar o campo espacial do modo.
  - no estado atual, os modos 3D saem como `E_cell_centroid_unit_peak_normalized`.
- `fields_csv_file`
  - nome do CSV com o campo vetorial reconstruído no centroide de cada tetraedro.
- `vtk_file`
  - nome do arquivo VTK tetraédrico correspondente ao mesmo modo.

## 4) Por que o matching funciona

Os casos 3D usam a mesma lógica didática que já aparecia nas famílias 2D:

- primeiro se extraem as primeiras raízes positivas `k0 = sqrt(lambda)`;
- depois se faz o casamento com as linhas da tabela analítica;
- quando há degenerescência analítica, o casamento preserva o grupo degenerado.

No código, isso acontece em:

- `first_positive_k0(...)`
- `first_positive_k0_candidates(...)`
- `match_by_reference_with_degeneracy(...)`
- `match_indices_by_reference_with_degeneracy(...)`

ambas em `src/fem3d/fem3d_compare.hpp`.

A ideia é simples: comparar apenas por ordem numérica pode trocar modos
degenerados de lugar. O agrupamento evita isso e torna a tabela reproduzível.

## 5) `fields.csv` e `vtk/`

Cada linha casada de `modes.csv` aponta para dois artefatos espaciais:

- `csv/<solver>_<caso>_modeXX_<modo>_E_fields.csv`
- `vtk/<solver>_<caso>_modeXX_<modo>_E.vtk`

O CSV de campo usa:

```text
tet_id,xc_m,yc_m,zc_m,Ex,Ey,Ez,Emag
```

Significado:

- `tet_id`
  - índice do tetraedro na malha.
- `xc_m`, `yc_m`, `zc_m`
  - coordenadas do centroide do tetraedro, em metros.
- `Ex`, `Ey`, `Ez`
  - componentes cartesianas do campo vetorial reconstruído no centroide.
- `Emag`
  - magnitude do vetor reconstruído.

Observação didática importante:

- a reconstrução espacial atual usa o autovetor de arestas do modo casado;
- em cada tetraedro, o campo é avaliado no centroide com as funções de base
  de Whitney 1-form;
- o modo inteiro é normalizado pelo pico de `Emag`, para destacar a forma
  espacial sem fixar amplitude física absoluta.

O VTK correspondente grava a malha tetraédrica completa com:

- `CELL_DATA`
- vetor `E`
- escalar `Emag`

Isso permite inspeção posterior em ParaView ou leitura por scripts próprios.

## 6) `linop/`

Os arquivos em `linop/` registram o problema matricial efetivamente resolvido:

- `..._S_crs.csv`
- `..._T_crs.csv`
- `..._eigenvalues.csv`
- `..._eigenvectors.csv`

Eles seguem a convenção descrita em:

- [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md)

Para o `FEM3D0`, `S` e `T` saem diretamente da montagem densa.

Para o `FEM3D1`, os arquivos em `linop/` representam a forma expandida
simétrica do problema armazenado esparsamente. Isso permite comparar os dois
solvers com o mesmo formato externo de auditoria.

## 7) Leitura recomendada

Para entender o papel desses arquivos na Sec. 3.1:

- [src/fem3d/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d/README.md)
- [src/fem3d0/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d0/README.md)
- [src/fem3d1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d1/README.md)
- [traducao/3_Problemas_Tridimensionais.md](traducao/3_Problemas_Tridimensionais.md)
