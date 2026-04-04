# HELMVEC2: Referencia dos CSVs de Saida

Este arquivo documenta os CSVs gerados hoje pelo executavel:

- `helmvec2_rect`

Eles devem ser lidos junto com:

- a traducao da Secao 2.2.3 em
  [docs/traducao/2.2.3_Determinação do número de onda.md](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md)
- a documentacao do modulo em
  [src/helmvec2/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec2/README.md)

## Arquivos gerados

Ao rodar:

```bash
./build/helmvec2_rect 10 6 6 --backend closed-form
```

o projeto grava em:

- `out/helmvec2/rect/run.log`
- `out/helmvec2/rect/run_timing.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_modes.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_candidates.csv`

## 1) `helmvec2_rect_modes.csv`

Este e o CSV principal da Tabela 8.

Cada linha representa um modo publicado na validacao da Figura 11 / Tabela 8,
ja depois do matching entre:

- a raiz FEM calculada;
- a referencia `HELMVEC2(ref)`;
- a referencia `Hayata(ref)`.

Cabecalho atual:

```text
mode,matched_candidate_rank,k0_fem_matched,k0l_fem_matched,ref_helmvec2,ref_hayata,error_percent_helmvec2,error_percent_hayata,match_status
```

### Significado das colunas

- `mode`
  - indice didatico da linha da Tabela 8.

- `matched_candidate_rank`
  - posicao, dentro da lista de candidatos fisicos filtrados, da raiz que foi
    escolhida para casar com aquela linha da tabela.

- `k0_fem_matched`
  - valor FEM casado para `k0`.
  - no caso padrao atual, como `L=1`, ele coincide numericamente com `k0L`.

- `k0l_fem_matched`
  - valor FEM casado em forma adimensional `k0 L`.

- `ref_helmvec2`
  - valor de referencia publicado na coluna `HELMVEC2` do artigo.

- `ref_hayata`
  - valor de referencia publicado para `Hayata et al.`.

- `error_percent_helmvec2`
  - erro percentual relativo entre `k0l_fem_matched` e `ref_helmvec2`.

- `error_percent_hayata`
  - erro percentual relativo entre `k0l_fem_matched` e `ref_hayata`.

- `match_status`
  - status do matching.
  - no estado atual, o valor esperado e `matched`.

## 2) `helmvec2_rect_candidates.csv`

Este CSV guarda a lista de candidatos espectrais do solve, antes da etapa
final de matching com a Tabela 8.

Ele e importante porque mostra:

- quais autovalores positivos reais apareceram;
- qual foi a ordem deles no espectro;
- quais passaram no filtro fisico de propagacao;
- quanto cada candidato concentra energia no bloco longitudinal `Ez`.

Cabecalho atual:

```text
candidate_rank,eig_index,k0,k0l,ez_ratio,passes_physical_filter
```

### Significado das colunas

- `candidate_rank`
  - ordem do candidato na lista espectral ordenada por `k0`.

- `eig_index`
  - indice original do autovalor/autovetor retornado pelo solver generalizado.

- `k0`
  - raiz positiva real correspondente ao autovalor filtrado.

- `k0l`
  - forma adimensional `k0 L`.

- `ez_ratio`
  - fracao de energia Euclidiana do autovetor no bloco `Ez`:

```text
ez_ratio = ||Ez||^2 / (||Et||^2 + ||Ez||^2)
```

- `passes_physical_filter`
  - indica se o candidato passou no filtro fisico:

```text
k0 > beta / sqrt(eps_max)
```

## 3) `run_timing.csv`

Esse arquivo registra:

- configuracao da rodada;
- dados geometricos e materiais do caso;
- quantidade de candidatos crus e fisicos;
- tempos de montagem, solve e pos-processamento.

Cabecalho atual:

```text
case_label,geometry,backend,debug_local_blocks,debug_candidates,beta,L,betaL,nx,ny,mesh_nodes,mesh_tris,eps_top,eps_bottom,mu_r,k0_min_phys,candidate_count_raw,candidate_count_phys,matched_mode_count,assembly_ms,solve_ms,post_ms,total_ms
```

## Como os dois CSVs se complementam

- `helmvec2_rect_candidates.csv`
  - explica de onde vieram os candidatos do problema espectral da Eq. `(119)`.

- `helmvec2_rect_modes.csv`
  - explica como esses candidatos foram usados para reproduzir a Tabela 8 do
    artigo.

Em outras palavras:

- `candidates.csv` mostra o espectro fisico disponivel;
- `modes.csv` mostra a tabela final de validacao.

## Observacao importante

Nesta etapa do repositorio, o `HELMVEC2` ainda nao gera script proprio de
imagens a partir desses CSVs. A saida principal da familia e:

- textual no terminal e em `run.log`
- tabular em `modes.csv` e `candidates.csv`
- temporal em `run_timing.csv`

## Referencias cruzadas

- teoria traduzida:
  [2.2.3 Determinação do Número de Onda](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md)
- implementacao:
  [src/helmvec2/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec2/README.md)
- visao operacional:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
