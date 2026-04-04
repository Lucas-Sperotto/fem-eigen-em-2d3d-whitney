# HELMVEC1: Referencia dos Campos do CSV de Modos

Este arquivo documenta o significado didatico das colunas do CSV de modos
gerado pelos executaveis da familia `HELMVEC1`:

- `mixed_rect`
- `mixed_circle`
- `mixed_coax`

Ele deve ser lido junto com:

- a traducao da Secao 2.2.2 do artigo em
  [docs/traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
- a documentacao de implementacao em
  [src/helmvec1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/README.md)

## Cabecalho atual

O CSV de modos usa hoje o cabecalho:

```text
formulation,dominant_block,component_label,family,mode_label,positive_rank,eig_index,m,n,p,ar_m,b_m,r_m,r1_m,r2_m,kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,kc_r_fem,kc_r_ana,kc_r1_fem,kc_r1_ana,error_percent,rho_abs,edge_energy,scalar_energy,dominant_energy_ratio,match_space,match_method,mode_status
```

As colunas geometricas e analiticas que nao se aplicam ao caso ficam vazias.

Exemplos:

- no caso retangular, aparecem `m`, `n`, `ar_m`, `b_m`, `kc_ana`,
  `kc_ar_fem`, `kc_ar_ana` e `error_percent`
- no caso circular, aparecem `m`, `p`, `r_m`, `kc_r_fem` e `kc_r_ana`
- no caso coaxial, aparecem `m`, `p`, `r1_m`, `r2_m`, `kc_r1_fem` e
  `kc_r1_ana`

## Significado de cada coluna

- `formulation`
  - indica se a linha veio da formulacao `E` ou da formulacao `H`.

- `dominant_block`
  - bloco que dominou energeticamente o autovetor.
  - valores atuais:
    - `edge`
    - `scalar`

- `component_label`
  - componente efetivamente associada ao bloco dominante.
  - no estado atual:
    - formulacao `E` + bloco `edge` -> `Et`
    - formulacao `E` + bloco `scalar` -> `Ez`
    - formulacao `H` + bloco `edge` -> `Ht`
    - formulacao `H` + bloco `scalar` -> `Hz`

- `family`
  - familia modal interpretada a partir da formulacao e do bloco dominante.
  - no estado atual:
    - `E + edge -> TE`
    - `E + scalar -> TM`
    - `H + edge -> TM`
    - `H + scalar -> TE`

- `mode_label`
  - rotulo textual do modo identificado por matching.
  - exemplos:
    - `TE_m1_n0`
    - `TM_m1_p1`

- `positive_rank`
  - ordem didatica do modo dentro do subconjunto exportado daquela combinacao
    de formulacao e bloco dominante.

- `eig_index`
  - posicao do autovalor/autovetor no espectro retornado pelo solver
    generalizado.

- `m`
  - primeiro indice modal analitico do modo identificado.

- `n`
  - segundo indice modal analitico do caso retangular.

- `p`
  - segundo indice modal analitico dos casos circular e coaxial.

- `ar_m`
  - dimensao `a` do guia retangular, em metros.

- `b_m`
  - dimensao `b` do guia retangular, em metros.

- `r_m`
  - raio do guia circular, em metros.

- `r1_m`
  - raio interno do guia coaxial, em metros.

- `r2_m`
  - raio externo do guia coaxial, em metros.

- `kc_fem`
  - numero de onda de corte calculado numericamente pelo FEM.

- `kc_ana`
  - numero de onda de corte analitico associado ao modo identificado pelo
    matching por correlacao de massa.

- `kc_ar_fem`
  - forma normalizada do cutoff numerico no caso retangular:

```text
kc_ar_fem = kc_fem * a
```

- `kc_ar_ana`
  - forma normalizada do cutoff analitico no caso retangular.

- `kc_r_fem`
  - forma normalizada do cutoff numerico no caso circular:

```text
kc_r_fem = kc_fem * r
```

- `kc_r_ana`
  - forma normalizada do cutoff analitico no caso circular.

- `kc_r1_fem`
  - forma normalizada do cutoff numerico no caso coaxial:

```text
kc_r1_fem = kc_fem * r1
```

- `kc_r1_ana`
  - forma normalizada do cutoff analitico no caso coaxial.

- `error_percent`
  - erro percentual entre `kc_fem` e `kc_ana`:

```text
error_percent = 100 * (kc_fem - kc_ana) / kc_ana
```

- `rho_abs`
  - correlacao modal em modulo entre o modo numerico dominante e a referencia
    analitica escolhida.
  - varia tipicamente entre `0` e `1`.
  - quanto mais perto de `1`, melhor o casamento modal.

- `edge_energy`
  - energia numerica do bloco de aresta do autovetor:

```text
edge_energy = ||x_edge||^2
```

- `scalar_energy`
  - energia numerica do bloco escalar do autovetor:

```text
scalar_energy = ||x_scalar||^2
```

- `dominant_energy_ratio`
  - fracao da energia total pertencente ao bloco dominante:

```text
dominant_energy_ratio =
    max(edge_energy, scalar_energy) / (edge_energy + scalar_energy)
```

- `match_space`
  - espaco em que o matching por correlacao de massa foi feito.
  - valores atuais:
    - `edge_Tt`
    - `scalar_Tz`

- `match_method`
  - metodo usado para identificar o rotulo analitico.
  - valor atual:
    - `mass_correlation`

- `mode_status`
  - status textual resumido do modo.
  - valores tipicos:
    - `edge_dominant`
    - `scalar_dominant`

## Como a classificacao por energia funciona

No `HELMVEC1`, o solve e feito sobre o sistema misto da Eq. `(92)`:

```text
[ St   0 ] [et] = kc^2 [ Tt   0 ] [et]
[  0  Sz ] [ez]        [  0  Tz ] [ez]
```

Depois do solve, cada autovetor e lido como:

```text
x = [x_edge ; x_scalar]
```

O projeto calcula:

```text
edge_energy   = ||x_edge||^2
scalar_energy = ||x_scalar||^2
```

e classifica o modo pelo bloco dominante.

No repositorio, isso aparece em:

- [mixed_mode_utils.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/mixed_mode_utils.hpp)

### Por que isso funciona

- No cutoff da Eq. `(92)`, o sistema global e bloco-diagonal.
- Portanto, os autovetores fisicos tendem naturalmente a se concentrar em um
  dos dois blocos.
- A comparacao entre `edge_energy` e `scalar_energy` permite separar
  numericamente os modos dominados pelo bloco transversal daqueles dominados
  pelo bloco longitudinal.
- Isso fornece uma classificacao didatica robusta sem precisar forcar uma
  separacao artificial durante o solve.

## Como `rho_abs` e calculado e por que funciona

No `HELMVEC1`, o matching nao e feito no autovetor misto completo. Ele e feito
no bloco dominante do modo:

- se `dominant_block = edge`, o codigo extrai `x_edge` e compara no espaco
  vetorial com a matriz de massa `Tt`
- se `dominant_block = scalar`, o codigo extrai `x_scalar` e compara no espaco
  escalar com a matriz de massa `Tz`

Em notacao compacta:

```text
rho_abs =
    | x_dom^T M_dom u_ref |
    ------------------------------------------
    sqrt(x_dom^T M_dom x_dom) sqrt(u_ref^T M_dom u_ref)
```

onde:

- `x_dom` e o bloco dominante do autovetor numerico
- `u_ref` e o modo analitico candidato
- `M_dom` e `Tt` ou `Tz`, dependendo de `match_space`

### Como isso aparece no codigo

- extracao do bloco dominante:
  [mixed_mode_match.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/mixed_mode_match.hpp)
- matching retangular:
  [main_mixed_rect.cpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/main_mixed_rect.cpp)
- matching circular:
  [main_mixed_circle.cpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/main_mixed_circle.cpp)
- matching coaxial:
  [main_mixed_coax.cpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/main_mixed_coax.cpp)

### Por que essa correlacao funciona

- ela respeita o produto interno natural do bloco em que o modo vive
- ela e independente do sinal global do autovetor
- ela e independente da escala arbitraria do autovetor
- o modo correto tende a maximizar essa correlacao entre os candidatos

Em outras palavras: em vez de perguntar apenas "qual `kc` parece proximo?",
o repositorio pergunta "qual referencia analitica tem a forma modal mais
parecida com este bloco numerico?".

### Observacao pratica importante

Em malhas muito grosseiras, especialmente em blocos escalares com modos
proximos ou degenerados, diferentes candidatos podem apresentar `rho_abs`
alto. Nesses casos:

- o matching continua sendo util como criterio principal
- mas a leitura correta precisa considerar tambem `kc_fem`, `kc_ana` e o
  refinamento da malha

Por isso, o `HELMVEC1` salva no mesmo CSV:

- a classificacao por energia
- o `rho_abs`
- a referencia analitica escolhida
- o erro percentual em `kc`

## O que este CSV responde

O CSV de modos do `HELMVEC1` responde:

- de qual formulacao (`E` ou `H`) veio o modo
- qual bloco foi dominante
- qual familia modal foi atribuida (`TE` ou `TM`)
- qual componente do problema misto esta sendo destacada (`Et`, `Ez`, `Ht`,
  `Hz`)
- qual foi o `kc` numerico obtido
- qual referencia analitica foi associada ao modo
- qual foi o `rho_abs` desse matching
- em qual espaco o matching foi feito (`edge_Tt` ou `scalar_Tz`)
- qual foi a energia relativa dos blocos

## Observacao importante

Ao contrario do `HELM10` e do `HELMVEC`, o `HELMVEC1` ainda nao exporta, por
padrao, um campo espacial por modo. Nesta etapa do repositorio, a saida
principal da familia e:

- tabular no terminal
- tabular no `modes.csv`
- temporal no `run_timing.csv`
- textual completa no `run.log`
- grafica em `img/`, por meio do script `scripts/helmvec1.py`, com resumos do
  cutoff normalizado, do `rho_abs`, da energia dominante, das energias de
  bloco e do erro quando houver referencia analitica

Para manter a saida didatica e o custo de pos-processamento sob controle, o
export atual do `modes.csv` fica limitado aos primeiros `20` modos de cada
subgrupo (`E/TE`, `E/TM`, `H/TE`, `H/TM`).

## Referencias cruzadas

- teoria traduzida:
  [2.2.2 Três Componentes](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
- implementacao:
  [src/helmvec1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/README.md)
- imagens de resumo:
  [HELMVEC1_Imagens_Referencia.md](HELMVEC1_Imagens_Referencia.md)
- visao operacional dos executaveis:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
