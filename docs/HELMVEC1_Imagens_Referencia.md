# HELMVEC1: Referencia das Imagens

Este arquivo documenta o script:

- [scripts/helmvec1.py](/home/sperotto/tp3485-fem-eigen-em/scripts/helmvec1.py)

Ele le os CSVs modais produzidos pelos executaveis da familia `HELMVEC1` e
gera imagens do sistema misto da Eq. `(92)`.

## Ideia principal

O `HELMVEC1` agora exporta dois tipos de imagem espacial, de acordo com o
bloco dominante do modo:

- `edge` dominante -> campo transversal vetorial
- `scalar` dominante -> componente longitudinal escalar

Isso preserva a leitura fisica correta da Eq. `(92)`.

## Como executar

Todos os casos:

```bash
python3 scripts/helmvec1.py
```

Apenas um caso:

```bash
python3 scripts/helmvec1.py --case rect
python3 scripts/helmvec1.py --case circle
python3 scripts/helmvec1.py --case coax
```

Com outra resolucao:

```bash
python3 scripts/helmvec1.py --dpi 120
```

## Imagens geradas

Para cada caso, o script gera:

- `helmvec1_<caso>_kc_by_mode.png`
  - resumo do cutoff normalizado por modo.
  - no retangular, usa `kc_ar_fem = kc_fem * a`.
  - no circular, usa `kc_r_fem = kc_fem * r`.
  - no coaxial, usa `kc_r1_fem = kc_fem * r1`.

- `helmvec1_<caso>_rho_by_mode.png`
  - resumo da correlacao modal em modulo `|rho|`.
  - usa a coluna `rho_abs` do `modes.csv`.
  - mostra o quao parecido o bloco dominante do autovetor numerico ficou da
    referencia analitica escolhida no matching.

- `helmvec1_<caso>_dominant_ratio_by_mode.png`
  - resumo da razao de energia dominante:

```text
dominant_energy_ratio =
    max(edge_energy, scalar_energy) / (edge_energy + scalar_energy)
```

- `helmvec1_<caso>_block_energy_by_mode.png`
  - comparacao direta entre `edge_energy` e `scalar_energy` por modo.

- `helmvec1_<caso>_error_by_mode.png`
  - resumo do erro percentual por modo.
  - como o matching analitico agora existe em `rect`, `circle` e `coax`,
    esse grafico fica disponivel nos tres casos.

- `img/magnitude/*.png`
  - para modos `edge`, mapa de magnitude do campo transversal por celula.

- `img/quiver/*.png`
  - para modos `edge`, diagrama de setas do campo transversal por centroide.

- `img/scalar/*.png`
  - para modos `scalar`, mapa escalar nodal com isolinhas da componente
    longitudinal (`Ez` ou `Hz`).

## Organizacao dos subgraficos

Cada figura e separada em subgraficos pelas combinacoes didaticas:

- `E / TE`
- `E / TM`
- `H / TE`
- `H / TM`

Dentro de cada subgrafico, os modos sao ordenados por `kc_fem` e, em empate,
por `positive_rank`.

No eixo `x`, o script usa um rotulo compacto do modo:

- `TE_m1_n0` vira `m1,n0`
- `TM_rank03` vira `rank03`

O titulo de cada subgrafico explicita:

- o caso (`RECT`, `CIRCLE`, `COAX`)
- a formulacao (`E` ou `H`)
- a familia modal (`TE` ou `TM`)
- a componente destacada (`Et`, `Ez`, `Ht`, `Hz`)
- o bloco dominante (`edge` ou `scalar`)

## Como ler essas figuras

### 1) Cutoff normalizado por modo

Este grafico responde:

- quais modos aparecem primeiro em cada formulacao;
- como o espectro se organiza dentro de cada familia;
- como comparar os valores numericos com as tabelas do artigo sem misturar
  unidades.

### 2) Razao de energia dominante

Este grafico responde:

- quao "puro" esta cada modo em relacao ao bloco dominante;
- se o autovetor esta quase todo no bloco de aresta ou no bloco escalar.

Nos casos homogeneos atuais do repositorio, a Eq. `(92)` aparece bloco-diagonal
e a razao dominante tende a `1.0`. Isso e um resultado esperado e didatico:
mostra que os modos ficam praticamente separados por bloco.

Quando isso acontece, o grafico registra explicitamente que todos os modos do
subgrupo ficaram com a mesma razao dominante, em vez de repetir dezenas de
anotacoes iguais sobre os pontos.

### 2.1) Correlacao modal `|rho|`

Este grafico responde:

- quao bem o modo numerico foi identificado pelo matching analitico;
- em quais subgrupos a correlacao ficou proxima de `1`;
- onde pode haver ambiguidade numerica em malhas mais grosseiras.

No `HELMVEC1`, esse `rho` nao vem do vetor misto completo. Ele vem do bloco
dominante:

- `edge` -> correlacao com `Tt`
- `scalar` -> correlacao com `Tz`

Isso faz o grafico dialogar diretamente com:

- `dominant_block`
- `match_space`
- `match_method`

### 3) Energias de bloco

Este grafico responde:

- quanto de energia caiu no bloco `edge`;
- quanto de energia caiu no bloco `scalar`;
- por que o modo foi classificado como `edge_dominant` ou
  `scalar_dominant`.

### 4) Erro por modo

Este grafico responde:

- quao perto o `kc_fem` ficou da referencia analitica;
- como o erro evolui ao longo da ordenacao modal em cada geometria.

### 5) Magnitude e quiver dos modos `edge`

Essas figuras aparecem quando `field_data_kind = edge_vector_cell`.

Elas mostram:

- a distribuicao espacial do campo transversal dominante;
- a direcao do campo reconstruido por elemento;
- se o modo identificado como `Et` ou `Ht` faz sentido fisico.

### 6) Mapa escalar dos modos `scalar`

Essas figuras aparecem quando `field_data_kind = scalar_nodal`.

Elas mostram:

- a componente longitudinal nodal (`Ez` ou `Hz`);
- suas isolinhas sobre a secao transversal;
- a forma espacial do bloco escalar que dominou o autovetor misto.

## Relacao com a documentacao teorica

Este material deve ser lido junto com:

- [HELMVEC1_CSV_Modos_Referencia.md](HELMVEC1_CSV_Modos_Referencia.md)
- [HELMVEC1_CSV_Campos_Referencia.md](HELMVEC1_CSV_Campos_Referencia.md)
- [src/helmvec1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/README.md)
- [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
