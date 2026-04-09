# `helmvec3` - Calculo de `beta` para `k0` dado (Sec. 2.2.4)

Este modulo resolve o problema complementar ao `helmvec2`:
- entrada: `k0` (frequencia),
- saida: `beta` (constante de propagacao).

## 1) Objetivo do modulo

Montar e resolver o EVP acoplado em `beta^2` para guia retangular
parcialmente preenchido, reproduzindo os estudos da Sec. 2.2.4
(Figuras 12 e 13, Tabelas 9 e 10).

Entradas publicas:
- `main_helmvec3_fig12_rect.cpp`
- `main_helmvec3_fig13_rect.cpp`

Implementacao compartilhada:
- `main_helmvec3_rect.cpp`

Montagem compartilhada:
- `src/helmvec2/helmvec2_coupled_system.cpp`
- utilitarios `src/helmvec2/helmvec23_shared.hpp`

## 2) Sistema matematico

Com `x = [Et; Ez]`:

`P x = beta^2 Q x`

Blocos implementados (arranjo de 2.2.4):

- `P_tt = St - k0^2 Tt`
- `P_zz = k0^2 Tz`

- `Q_tt = -Mt_(1/mu)`
- `Q_tz = +C`
- `Q_zt = +C_orient^T`
- `Q_zz = Sz`

Elementos:
- `St`, `Tt`: blocos edge,
- `Sz`, `Tz`: blocos escalares,
- `Mt_(1/mu)`: massa vetorial com peso `1/mu_r`,
- `C(m,j) = int (1/mu_r) W_m . grad(N_j) dA`.

Na notacao local do artigo:

```text
Sel(tt)  -> Eq. (137)
Tel(tz)  -> Eq. (138)
Tel(zt)  -> Eq. (139)
Sel(zz)  -> Eq. (140)
Tel(tt)  -> Eq. (141)
Tel(zz)  -> Eq. (142)
```

No arranjo numerico implementado:

```text
P_tt = St - k0^2 Tt
P_zz = k0^2 Tz

Q_tt = -(1/mu_r) Mt
Q_tz = (1/mu_r) C
Q_zt = (1/mu_r) C^T
Q_zz = (1/mu_r) Gz
```

onde:

```text
Gz(i,j) = int_T grad(N_i) . grad(N_j) dA
Mt(m,n) = int_T W_m . W_n dA
Tz(i,j) = int_T N_i N_j dA
```

## 2.1) Forma fraca por bloco

Partindo das Eq. `(131)` a `(135)` do artigo, os blocos locais podem ser lidos
como:

```text
Sel(tt)(m,n) =
    (1/mu_r) int_T curl(W_m) curl(W_n) dA
  - k0^2 eps_r int_T W_m . W_n dA

Tel(tz)(m,j) =
    (1/mu_r) int_T W_m . grad(N_j) dA

Tel(zt)(i,n) =
    (1/mu_r) int_T grad(N_i) . W_n dA

Sel(zz)(i,j) =
    (1/mu_r) int_T grad(N_i) . grad(N_j) dA

Tel(tt)(m,n) =
    eps_r int_T W_m . W_n dA

Tel(zz)(i,j) =
    (1/mu_r) int_T grad(N_i) . grad(N_j) dA
  + k0^2 eps_r int_T N_i N_j dA
```

As Eq. `(137)` a `(142)` sao justamente as versoes `closed-form` desses blocos.

## 2.2) Forma fechada implementada

No caminho `--backend closed-form`, o modulo nao depende apenas de reuso
implicito dos montadores base. A formulacao local de `2.2.4` foi deixada
explicitamente nomeada em:

- `src/explicit/tri2d_coupled_explicit.hpp`

Funcoes principais:

- `tri2d_beta_closed_form_eq_137_142(...)`
  - monta diretamente os seis blocos locais do artigo:
  - `Sel(tt)`, `Tel(tz)`, `Tel(zt)`, `Sel(zz)`, `Tel(tt)`, `Tel(zz)`.

- `tri2d_beta_rearranged_closed_form_eq_136(...)`
  - rearranja esses blocos para a forma usada no codigo:
  - `P x = beta^2 Q x`.

Essa segunda funcao e a mais importante para a montagem global, porque o
repositorio resolve a Secao `2.2.4` ja na forma rearranjada da Eq. `(136)`.
Assim, o codigo fica ao mesmo tempo:

- fiel ao artigo no nivel dos blocos locais `(137)` a `(142)`;
- fiel ao fluxo numerico validado do repositorio no nivel da montagem global.

Observacao importante:
- a forma impressa do artigo nem sempre e a melhor representacao para a
  implementacao global;
- por isso, o codigo privilegia a equivalencia algebraica validada e deixa os
  comentarios apontando onde cada equacao entra.
- nesta etapa, essa rastreabilidade tambem passou a existir diretamente na
  estrutura do sistema global, antes do fechamento final de `P` e `Q`.

## 2.3) Expressoes locais explicitas

No backend `closed-form`, os blocos locais ficam:

```text
Sel(tt) = (1/mu_r) * curlcurl_t - k0^2 * eps_r * mass_t
Tel(tz) = (1/mu_r) * C
Tel(zt) = (1/mu_r) * C^T
Sel(zz) = (1/mu_r) * gradgrad_z
Tel(tt) = eps_r * mass_t
Tel(zz) = (1/mu_r) * gradgrad_z + k0^2 * eps_r * mass_z
```

Depois, esses blocos sao rearranjados para a Eq. `(136)`.

## 2.4) Onde a Eq. (136) e montada no codigo

A montagem global da Eq. `(136)` ocorre em:

- `src/helmvec2/helmvec2_coupled_system.cpp`
- funcao `build_coupled_beta_system_E(...)`

No backend `closed-form`, a montagem elemento a elemento passa por:

- `assemble_beta_system_closed_form(...)`

e cada triangulo usa:

- `tri2d_beta_rearranged_closed_form_eq_136(...)`

Ou seja, a cadeia didatica fica:

1. `Eq. (137)` a `Eq. (142)` -> helper local explicito
2. rearranjo para `Eq. (136)` -> helper local rearranjado
3. blocos globais nomeados -> `P_tt`, `P_zz`, `Q_tt`, `Q_tz`, `Q_zt`, `Q_zz`
4. fechamento global explicito -> `assemble_eq136_global_from_named_blocks(...)`
5. entrada publica didatica -> `tp3485::build_eq136_helmvec3_system_E(...)`
6. assembleia global subjacente -> `build_coupled_beta_system_E(...)`

Os blocos nomeados ficam preservados em:

- `CoupledBetaSystem::{P_tt, P_zz, Q_tt, Q_tz, Q_zt, Q_zz}`
- arquivo `src/helmvec2/helmvec2_coupled_system.hpp`

O fechamento global explicito passa por:

- `initialize_named_eq136_global_blocks(...)`
- `assemble_beta_system_closed_form(...)`
- `assemble_coupling_block_C_named(...)`
- `assemble_eq136_global_from_named_blocks(...)`

## 3) Geometrias/casos de validacao

Parametros base:
- `a = 1.0`, `b/a = 0.45`,
- `eps_fill = 2.45`.

### Figura 12 / Tabela 9

Preenchimento horizontal inferior:
- interface em `y = d = 0.5*b`.

Compara `beta/k0` para:
- `b/lambda0 = {0.2, 0.3, 0.4, 0.5, 0.6}`.

### Figura 13 / Tabela 10

Preenchimento vertical a esquerda:
- interface em `x = d`, com `d/a` variavel.

Compara blocos em funcao de `d/a` e `b/lambda0`.

## 4) Solver e filtro fisico

Solver:
- `generalized_eigs_real_vec` (`dggev`).

Filtro no pos-processamento:
- parte imaginaria pequena,
- `beta^2 > 0`,
- faixa fisica de `beta/k0` limitada por `sqrt(eps_max)` com margem numerica.

## 5) Matching modal

Estrategias implementadas:
- matching por proximidade ao valor de referencia analitica,
- opcao de rastreamento de ramo continuo (`trace_ratio_branch`) para inspecao.

Motivacao:
- reduzir erro de associacao quando ha troca de ordem entre ramos numericos.

## 6) Uso

```bash
./build/helmvec3_fig12_rect 10 5
./build/helmvec3_fig12_rect 10 5 1
./build/helmvec3_fig12_rect 10 5 --debug-local-blocks --backend closed-form
./build/helmvec3_fig12_rect 10 5 --debug-candidates
# args: nx ny [debug]

./build/helmvec3_fig13_rect 0.20 10 5
./build/helmvec3_fig13_rect 0.20 10 5 1
./build/helmvec3_fig13_rect 0.20 10 5 --debug-local-blocks --backend closed-form
./build/helmvec3_fig13_rect 0.20 10 5 --debug-candidates
# args: d_over_a_preview nx ny [debug]
```

Saida textual:

- `helmvec3_fig12_rect`
  - Tabela 9 / Figura 12.
- `helmvec3_fig13_rect`
  - preview de ramo da Figura 13 para um `d/a` escolhido;
  - validacao completa da Tabela 10.

Saida estruturada em disco:

```text
out/helmvec3/fig12_rect/
  run.log
  run_timing.csv
  csv/
    helmvec3_fig12_rect_table9.csv
    helmvec3_fig12_rect_figure12_br*_Et_fields.csv
    helmvec3_fig12_rect_figure12_br*_Ez_fields.csv
  vtk/
    helmvec3_fig12_rect_figure12_br*_Et.vtk
    helmvec3_fig12_rect_figure12_br*_Ez.vtk
  linop/
    helmvec3_fig12_rect_figure12_*.csv

out/helmvec3/fig13_rect/
  run.log
  run_timing.csv
  csv/
    helmvec3_fig13_rect_preview.csv
    helmvec3_fig13_rect_table10.csv
    helmvec3_fig13_rect_preview_da*_br*_Et_fields.csv
    helmvec3_fig13_rect_preview_da*_br*_Ez_fields.csv
    helmvec3_fig13_rect_table10_da*_br*_Et_fields.csv
    helmvec3_fig13_rect_table10_da*_br*_Ez_fields.csv
  vtk/
    helmvec3_fig13_rect_preview_da*_br*_Et.vtk
    helmvec3_fig13_rect_preview_da*_br*_Ez.vtk
    helmvec3_fig13_rect_table10_da*_br*_Et.vtk
    helmvec3_fig13_rect_table10_da*_br*_Ez.vtk
  linop/
    helmvec3_fig13_rect_preview_da*_*.csv
    helmvec3_fig13_rect_table10_da*_*.csv
```

O `run.log` espelha tudo que o executavel imprime em tela. O `run_timing.csv`
registra a configuracao da rodada, o tamanho da malha, a quantidade de pontos
resolvidos nas Tabelas 9 e 10 e os tempos agregados de montagem, solve e
pos-processamento.

Os CSVs principais se dividem assim:

- `helmvec3_fig12_rect_table9.csv`
  - curva casada da Figura 12 / Tabela 9;
  - cada linha agora registra tambem o candidato/autovetor escolhido e a
    ponte para os artefatos espaciais `Et`/`Ez`.
- `helmvec3_fig13_rect_preview.csv`
  - ramo rastreado por continuidade para um `d/a` escolhido;
  - cada linha agora registra o candidato acompanhado no ramo e a ponte para
    os artefatos espaciais `Et`/`Ez`.
- `helmvec3_fig13_rect_table10.csv`
  - validacao completa por blocos da Figura 13 / Tabela 10;
  - cada linha agora registra o candidato/autovetor escolhido e a ponte para
    os artefatos espaciais `Et`/`Ez`.

Como a Eq. `(136)` usa `x = [Et ; Ez]`, cada ponto exportado gera:

- um CSV e um VTK de `Et` por celula;
- um CSV e um VTK de `Ez` por no.

Os artefatos em `linop/` seguem outra granularidade:

- um pacote espectral por ponto resolvido da curva;
- cada pacote guarda `P`, `Q`, os blocos nomeados da Eq. `(136)`,
  `eigenvalues.csv` e `eigenvectors.csv`;
- isso permite auditar posteriormente cada problema generalizado resolvido
  ao longo da Figura 12, do preview e da Tabela 10.

Para gerar as imagens a partir desses CSVs:

```bash
python3 scripts/helmvec3.py
python3 scripts/helmvec3.py --case fig12_rect --dpi 120
python3 scripts/helmvec3.py --case fig13_rect --dpi 120 --show-mesh
```

Observacao:

- o script ja acompanha as novas entradas publicas
  `helmvec3_fig12_rect` e `helmvec3_fig13_rect`;
- o nucleo numerico, os CSVs/VTKs e a camada de imagens ficaram separados.

O script gera:

- os graficos-resumo da Figura 12 / Tabela 9, do preview e da Tabela 10;
- mapas de magnitude de `Et` em `img/magnitude/`;
- diagramas `quiver` de `Et` em `img/quiver/`;
- mapas escalares de `Ez` em `img/scalar/`.

Nas tres arvores espaciais, os arquivos ficam organizados em subpastas por
caso geometrico `d/a`, por exemplo:

- `img/magnitude/figure12_da_0_225/`
- `img/quiver/preview_da_0_2/`
- `img/scalar/table10_da_0_5/`

Quando `debug=1`, o driver tambem imprime os blocos locais do primeiro
triangulo para inspecao didatica:

- blocos do artigo `Eq. (137)` a `Eq. (142)`;
- blocos rearranjados usados no sistema global `Eq. (136)`.

Flags de depuracao:
- `--debug-local-blocks`: imprime apenas os blocos locais e o rearranjo da
  Eq. `(136)`
- `--debug-candidates`: imprime os candidatos modais usados no matching
- `debug=1` legado: ativa os dois comportamentos

Para acompanhar a saida bruta via script:

```bash
python3 scripts/validate_2d_22.py --backend closed-form --show-output --debug-candidates
```

## 7) Integracao com scripts

`scripts/validate_2d_22.py` ja foi adaptado para os binarios
`helmvec3_fig12_rect` e `helmvec3_fig13_rect`.
Ele grava em:
- `out/validation/validation_2d_22.csv`

As figuras consolidadas desta secao sao geradas depois por:
- `python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22`

Campos relacionados:
- `section=2.2.4`
- `case=Figure12_Table9` e `case=Figure13_Table10`
- erro relativo contra referencia analitica e HELMVEC3.

## 8) Relacao com a sequencia 2D

- `helmvec1`: `kc` no cutoff
- `helmvec2`: `k0` com `beta` dado
- `helmvec3`: `beta` com `k0` dado

A base de montagem e blocos FEM permanece a mesma; muda apenas o rearranjo
algebrico entre operador da esquerda/direita no EVP generalizado.

## 9) Comparacao entre `gauss` e `closed-form`

O executavel aceita:

- `--backend closed-form`
- `--backend gauss`

Sem `--backend`, o padrao publico do executavel e `closed-form`.

Isso permite comparar:

- autovalores;
- ramos `beta/k0`;
- tempo de montagem e tempo total;
- consistencia entre a integracao numerica e a forma fechada.

Na pratica:

- `closed-form` e util para rastreabilidade matematica e comparacao direta com
  as equacoes do artigo;
- `gauss` e util para manter continuidade com a montagem numerica tradicional
  como referencia auxiliar.

## 10) Referencias didaticas desta familia

- [docs/HELMVEC3_CSV_Referencia.md](/home/sperotto/tp3485-fem-eigen-em/docs/HELMVEC3_CSV_Referencia.md)
- [docs/HELMVEC3_Imagens_Referencia.md](/home/sperotto/tp3485-fem-eigen-em/docs/HELMVEC3_Imagens_Referencia.md)
- [docs/Artefatos_Espectrais_CSV_Referencia.md](/home/sperotto/tp3485-fem-eigen-em/docs/Artefatos_Espectrais_CSV_Referencia.md)
- [docs/traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](/home/sperotto/tp3485-fem-eigen-em/docs/traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
