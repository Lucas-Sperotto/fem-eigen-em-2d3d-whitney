# Artefatos Espectrais em CSV

Este arquivo documenta os artefatos espectrais exportados pelas familias:

- `HELM10`
- `HELMVEC`
- `HELMVEC1`
- `HELMVEC2`
- `HELMVEC3`

O objetivo e deixar audivel, para verificacoes externas, nao apenas a saida
modal resumida, mas tambem:

- as matrizes globais do problema de autovalor;
- os autovalores ordenados;
- os autovetores na mesma ordem.

## Onde esses arquivos aparecem

Em cada caso numerico, os novos artefatos ficam em:

```text
out/<familia>/<caso>/linop/
```

Exemplos:

- `out/helm10/rect/linop/`
- `out/helmvec/circle/linop/`
- `out/helmvec1/coax/linop/`
- `out/helmvec2/rect/linop/`
- `out/helmvec3/rect/linop/`

## 1) Matrizes em CRS

As matrizes sao salvas em CSV no formato CRS, com cabecalho:

```text
section,index,value
```

As linhas aparecem em blocos:

- `nrows`
  - numero de linhas da matriz.
- `ncols`
  - numero de colunas da matriz.
- `nnz`
  - numero de entradas nao nulas exportadas.
- `row_ptr`
  - vetor `row_ptr` do CRS, com tamanho `nrows + 1`.
- `col_idx`
  - vetor de indices de coluna dos nao nulos.
- `values`
  - vetor de valores correspondente a `col_idx`.

Em outras palavras, um arquivo como:

```text
section,index,value
nrows,0,25
ncols,0,25
nnz,0,105
row_ptr,0,0
row_ptr,1,3
...
col_idx,0,0
...
values,0,1.2345
...
```

permite reconstruir a matriz global sem depender do binario do projeto.

Observacao:

- o exportador considera nao nulo todo valor com `|a_ij| > 1e-14`;
- isso evita poluir o arquivo com zeros numericos residuais.

## 2) Autovalores ordenados

Os autovalores sao salvos em:

- `*_eigenvalues.csv`

### Problemas simetricos

Cabecalho:

```text
ordered_rank,solver_index,lambda
```

Interpretacao:

- `ordered_rank`
  - posicao do autopar depois da ordenacao.
- `solver_index`
  - indice original devolvido pelo solver LAPACK.
- `lambda`
  - autovalor.

Ordenacao usada:

- ordem crescente de `lambda`;
- em empate numerico, desempate por `solver_index`.

### Problemas gerais nao simetricos

Cabecalho:

```text
ordered_rank,solver_index,lambda_real,lambda_imag,beta_denominator
```

Interpretacao:

- `lambda_real`, `lambda_imag`
  - partes real e imaginaria do autovalor generalizado.
- `beta_denominator`
  - denominador retornado por `dggev`, util para auditoria numerica.

Ordenacao usada:

- crescente em `lambda_real`;
- depois crescente em `lambda_imag`;
- depois `solver_index`.

## 3) Autovetores ordenados

Os autovetores sao salvos em:

- `*_eigenvectors.csv`

Cabecalho comum:

```text
ordered_rank,solver_index,dof_index,component_real,component_imag
```

Interpretacao:

- `ordered_rank`
  - aponta para a mesma linha de `*_eigenvalues.csv`.
- `solver_index`
  - indice original do solver.
- `dof_index`
  - grau de liberdade da componente.
- `component_real`
  - parte real do autovetor naquele DOF.
- `component_imag`
  - parte imaginaria do autovetor naquele DOF.

### Problemas simetricos

Neles, `component_imag = 0` para todas as linhas.

### Problemas gerais nao simetricos

Para os casos resolvidos com `dggev`, o projeto reconstrui a forma complexa dos
autovetores a partir do armazenamento real do LAPACK:

- autovalor real
  - `component_imag = 0`;
- par complexo conjugado
  - o CSV ja sai com as partes real e imaginaria separadas por DOF.

Isso evita que a verificacao externa precise reinterpretar manualmente a
convencao de empacotamento do `LAPACK`.

## 4) Convencoes por familia

### `HELM10`

Cada execucao exporta um pacote por ramo:

- `..._te_S_crs.csv`, `..._te_T_crs.csv`
- `..._te_eigenvalues.csv`, `..._te_eigenvectors.csv`
- `..._tm_S_crs.csv`, `..._tm_T_crs.csv`
- `..._tm_eigenvalues.csv`, `..._tm_eigenvectors.csv`

### `HELMVEC`

Segue a mesma ideia do `HELM10`, tambem separando:

- ramo `TE`
- ramo `TM`

### `HELMVEC1`

Exporta:

- o problema global `S/T` da formulacao `E`;
- o problema global `S/T` da formulacao `H`;
- os blocos nomeados da Eq. `(92)`:
  - `St`
  - `Tt`
  - `Sz`
  - `Tz`

### `HELMVEC2`

Exporta:

- o problema global `A/B` da Eq. `(119)`;
- os blocos nomeados:
  - `A_tt`, `A_tz`, `A_zt`, `A_zz`
  - `B_tt`, `B_zz`

### `HELMVEC3`

Como uma execucao do `HELMVEC3` resolve varios problemas ao longo das curvas de
dispersao, os artefatos sao exportados por ponto amostrado.

Exemplos de prefixos:

- `helmvec3_rect_figure12_br0_4_*`
- `helmvec3_rect_preview_da0_2_br0_5_*`
- `helmvec3_rect_table10_da0_375_br0_8_*`

Em cada prefixo, o par global sai como:

- `..._P_crs.csv`
- `..._Q_crs.csv`

e os blocos nomeados da Eq. `(136)` saem como:

- `P_tt`, `P_zz`
- `Q_tt`, `Q_tz`, `Q_zt`, `Q_zz`

## 5) Relacao com os outros CSVs do projeto

Esses artefatos nao substituem:

- `modes.csv`
- `fields_<modo>.csv`
- `candidates.csv`
- `run_timing.csv`

Eles complementam essas saidas. A ideia didatica fica:

- `modes.csv`
  - responde "quais modos foram identificados?".
- `fields_<modo>.csv`
  - responde "como o campo foi reconstruido?".
- `linop/*.csv`
  - responde "qual problema matricial foi efetivamente resolvido?".

## Referencias cruzadas

- [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
- [Rastreabilidade_Equacoes_Artigo_Codigo.md](Rastreabilidade_Equacoes_Artigo_Codigo.md)
- [src/helmvec3/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec3/README.md)
