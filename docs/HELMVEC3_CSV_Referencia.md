# HELMVEC3: Referencia dos CSVs de Saida

Este arquivo documenta os CSVs gerados hoje pelos executaveis:

- `helmvec3_fig12_rect`
- `helmvec3_fig13_rect`

Eles devem ser lidos junto com:

- a traducao da Secao 2.2.4 em
  [docs/traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
- a documentacao do modulo em
  [src/helmvec3/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec3/README.md)

## Arquivos gerados

Ao rodar:

```bash
./build/helmvec3_fig12_rect --nx 10 --ny 5 --backend closed-form
./build/helmvec3_fig13_rect --d-over-a-preview 0.20 --nx 10 --ny 5 --backend closed-form
```

o projeto grava em:

- `out/helmvec3/fig12_rect/run.log`
- `out/helmvec3/fig12_rect/run_timing.csv`
- `out/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_table9.csv`
- `out/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_figure12_br*_Et_fields.csv`
- `out/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_figure12_br*_Ez_fields.csv`
- `out/helmvec3/fig13_rect/run.log`
- `out/helmvec3/fig13_rect/run_timing.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_preview.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_preview_da*_br*_Et_fields.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_preview_da*_br*_Ez_fields.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10_da*_br*_Et_fields.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10_da*_br*_Ez_fields.csv`
- `out/helmvec3/fig12_rect/vtk/*.vtk`
- `out/helmvec3/fig13_rect/vtk/*.vtk`

## 1) `helmvec3_fig12_rect_table9.csv`

Este e o CSV principal da Figura 12 / Tabela 9.

Cada linha representa um ponto da curva de dispersao do caso com interface
horizontal, ja depois do matching entre:

- a razao FEM `beta/k0`;
- a referencia analitica;
- a referencia publicada pelo `HELMVEC3`.

Cabecalho atual:

```text
br_over_lambda0,beta_over_k0_fem,beta_over_k0_analytic,beta_over_k0_helmvec3,selected_candidate_rank,selected_eig_index,ez_ratio,error_percent_analytic,error_percent_helmvec3,match_status,field_status,et_fields_csv_file,ez_fields_csv_file,et_vtk_file,ez_vtk_file
```

### Significado das colunas

- `br_over_lambda0`
  - ponto amostrado do eixo horizontal da Figura 12.

- `beta_over_k0_fem`
  - valor FEM casado em forma adimensional.

- `beta_over_k0_analytic`
  - referencia analitica publicada no artigo.

- `beta_over_k0_helmvec3`
  - referencia publicada na coluna `HELMVEC3`.

- `selected_candidate_rank`
  - posicao do candidato escolhido dentro da lista fisica deduplicada daquele
    ponto amostrado.

- `selected_eig_index`
  - indice do autovalor/autovetor no solve generalizado `P x = beta^2 Q x`
    que representou aquele ponto exportado.

- `ez_ratio`
  - fracao de energia Euclidiana do autovetor no bloco longitudinal `Ez`:

```text
ez_ratio = ||Ez||^2 / (||Et||^2 + ||Ez||^2)
```

- `error_percent_analytic`
  - erro percentual relativo absoluto entre `beta_over_k0_fem` e
    `beta_over_k0_analytic`.

- `error_percent_helmvec3`
  - erro percentual relativo absoluto entre `beta_over_k0_fem` e
    `beta_over_k0_helmvec3`.

- `match_status`
  - status do matching.
  - no estado atual, o valor esperado e `matched`.

- `field_status`
  - descricao curta do tipo de normalizacao usada nos artefatos espaciais.

- `et_fields_csv_file`, `ez_fields_csv_file`
  - nomes dos CSVs espaciais por ponto para `Et` e `Ez`.

- `et_vtk_file`, `ez_vtk_file`
  - nomes dos VTKs espaciais por ponto para `Et` e `Ez`.

## 2) `helmvec3_fig13_rect_preview.csv`

Este CSV registra o preview de ramo da Figura 13 para um unico `d/a` escolhido
na linha de comando.

Cabecalho atual:

```text
d_over_a_preview,br_over_lambda0,beta_over_k0_fem_branch,selected_candidate_rank,selected_eig_index,ez_ratio,branch_status,field_status,et_fields_csv_file,ez_fields_csv_file,et_vtk_file,ez_vtk_file
```

### Significado das colunas

- `d_over_a_preview`
  - valor de `d/a` usado na pre-visualizacao.

- `br_over_lambda0`
  - ponto amostrado da curva.

- `beta_over_k0_fem_branch`
  - valor FEM do ramo rastreado por continuidade modal.

- `selected_candidate_rank`
  - posicao do candidato escolhido dentro da lista fisica deduplicada daquele
    ponto do preview.

- `selected_eig_index`
  - indice do autovetor usado para seguir o ramo naquele ponto.

- `ez_ratio`
  - fracao de energia do bloco `Ez` no autovetor acompanhado.

- `branch_status`
  - status do rastreamento.
  - no estado atual, o valor esperado e `tracked_branch`.

- `field_status`, `et_fields_csv_file`, `ez_fields_csv_file`,
  `et_vtk_file`, `ez_vtk_file`
  - metadados e artefatos espaciais do mesmo ponto do ramo.

## 3) `helmvec3_fig13_rect_table10.csv`

Este CSV guarda a validacao completa da Figura 13 / Tabela 10.

Cada linha representa um ponto de uma curva correspondente a um bloco `d/a`.

Cabecalho atual:

```text
d_over_a,br_over_lambda0,beta_over_k0_fem_matched,beta_over_k0_analytic,beta_over_k0_helmvec3,selected_candidate_rank,selected_eig_index,ez_ratio,error_percent_analytic,error_percent_helmvec3,match_status,field_status,et_fields_csv_file,ez_fields_csv_file,et_vtk_file,ez_vtk_file
```

### Significado das colunas

- `d_over_a`
  - parametro geometrico do bloco da Tabela 10.

- `br_over_lambda0`
  - ponto amostrado da curva naquele bloco.

- `beta_over_k0_fem_matched`
  - valor FEM casado com a referencia analitica correspondente.

- `beta_over_k0_analytic`
  - referencia analitica publicada na tabela.

- `beta_over_k0_helmvec3`
  - referencia publicada na coluna `HELMVEC3`.

- `selected_candidate_rank`
  - posicao do candidato escolhido dentro da lista fisica deduplicada daquele
    ponto da Tabela 10.

- `selected_eig_index`
  - indice do autovetor usado para representar aquele ponto.

- `ez_ratio`
  - fracao de energia do bloco `Ez` no autovetor escolhido.

- `error_percent_analytic`
  - erro percentual relativo absoluto entre `beta_over_k0_fem_matched` e
    `beta_over_k0_analytic`.

- `error_percent_helmvec3`
  - erro percentual relativo absoluto entre `beta_over_k0_fem_matched` e
    `beta_over_k0_helmvec3`.

- `match_status`
  - status do matching.
  - no estado atual, o valor esperado e `matched`.

- `field_status`, `et_fields_csv_file`, `ez_fields_csv_file`,
  `et_vtk_file`, `ez_vtk_file`
  - metadados e artefatos espaciais do mesmo ponto exportado.

## 4) CSVs espaciais `Et` e `Ez`

Cada ponto exportado do `HELMVEC3` agora gera dois CSVs complementares, porque
a Eq. `(136)` usa `x = [Et ; Ez]`.

### 4.1) `*_Et_fields.csv`

Cabecalho atual:

```text
cell_id,xc_m,yc_m,Ex,Ey,Emag
```

- `cell_id`
  - indice da celula triangular.
- `xc_m`, `yc_m`
  - coordenadas do centroide da celula.
- `Ex`, `Ey`
  - componentes cartesianas do campo transversal `Et`.
- `Emag`
  - magnitude local `sqrt(Ex^2 + Ey^2)`.

### 4.2) `*_Ez_fields.csv`

Cabecalho atual:

```text
node_id,x_m,y_m,Ez
```

- `node_id`
  - indice do no.
- `x_m`, `y_m`
  - coordenadas do no.
- `Ez`
  - componente longitudinal escalar no no.

## 5) `run_timing.csv`

Esse arquivo registra:

- configuracao da rodada;
- tamanho da malha;
- quantidade de pontos resolvidos nas Tabelas 9 e 10;
- tempos agregados de montagem, solve e pos-processamento.

Cabecalho atual:

```text
case_label,geometry,backend,debug_local_blocks,debug_candidates,a,b,d12,d13_preview_over_a,eps_fill,nx,ny,mesh_nodes,mesh_tris,table9_sample_count,preview_sample_count,table10_block_count,table10_row_count,assembly_ms,solve_ms,post_ms,total_ms
```

## Como os CSVs se complementam

- `helmvec3_fig12_rect_table9.csv`
  - explica a comparacao estatica da Figura 12 / Tabela 9.

- `helmvec3_fig13_rect_preview.csv`
  - explica o preview de um ramo continuo da Figura 13 para um `d/a`
    escolhido.

- `helmvec3_fig13_rect_table10.csv`
  - explica a validacao completa por blocos da Tabela 10.

- `*_Et_fields.csv` e `*_Ez_fields.csv`
  - mostram a decomposicao espacial do mesmo ponto exportado em `Et` e `Ez`.

Em outras palavras:

- `table9.csv` mostra o caso horizontal;
- `preview.csv` mostra um ramo de inspecao;
- `table10.csv` mostra a familia completa de comparacoes da Figura 13.
- os CSVs de `Et`/`Ez` mostram como cada ponto publicado se materializa no
  dominio espacial.

## Observacao importante

Nesta etapa do repositorio, o `HELMVEC3` ja gera um script proprio de imagens
a partir desses CSVs:

- [HELMVEC3_Imagens_Referencia.md](HELMVEC3_Imagens_Referencia.md)

A saida principal da familia passa a ser:

- textual no terminal e em `run.log`
- tabular em `table9.csv`, `preview.csv`, `table10.csv` e nos CSVs espaciais
  `Et`/`Ez`
- temporal em `run_timing.csv`
- grafica em `img/`, por meio de `python3 scripts/helmvec3.py`
- vetorial/escalares em `vtk/` para inspecao externa

## Referencias cruzadas

- teoria traduzida:
  [2.2.4 Caracteristicas de Dispersao](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
- implementacao:
  [src/helmvec3/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec3/README.md)
- imagens:
  [HELMVEC3_Imagens_Referencia.md](HELMVEC3_Imagens_Referencia.md)
- visao operacional:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
