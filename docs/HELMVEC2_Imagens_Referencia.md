# HELMVEC2: Referencia das Imagens

Este arquivo documenta o script:

- [scripts/helmvec2.py](/home/sperotto/tp3485-fem-eigen-em/scripts/helmvec2.py)

Ele le os CSVs produzidos hoje pelo executavel:

- `helmvec2_rect`

e gera imagens da Secao 2.2.3 a partir da Eq. `(119)`.

## Ideia principal

No `HELMVEC2`, o foco didatico continua sendo:

- o espectro de candidatos `k0`
- o filtro fisico de propagacao
- o matching final com a Figura 11 / Tabela 8

Mas, como a Eq. `(119)` usa o vetor de estado `x = [Et ; Ez]`, os modos
casados agora tambem exportam artefatos espaciais por modo. Por isso, o
script trabalha sobre:

- `out/helmvec2/rect/csv/helmvec2_rect_modes.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_candidates.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_modeXX_candYY_Et_fields.csv`
- `out/helmvec2/rect/csv/helmvec2_rect_modeXX_candYY_Ez_fields.csv`
- `out/helmvec2/rect/vtk/helmvec2_rect_modeXX_candYY_Et.vtk`
- `out/helmvec2/rect/vtk/helmvec2_rect_modeXX_candYY_Ez.vtk`
- `out/helmvec2/rect/run_timing.csv`

e salva as imagens em:

- `out/helmvec2/rect/img/`

## Como executar

```bash
python3 scripts/helmvec2.py
python3 scripts/helmvec2.py --case rect
python3 scripts/helmvec2.py --dpi 120
```

## Imagens geradas

- `helmvec2_rect_table8_k0l_by_mode.png`
  - compara `k0L(FEM matched)` com as referencias `HELMVEC2(ref)` e
    `Hayata(ref)`.
  - e a leitura grafica direta da Figura 11 / Tabela 8.

- `helmvec2_rect_error_by_mode.png`
  - mostra o erro relativo absoluto por modo contra as duas referencias.

- `helmvec2_rect_candidates_k0l_by_rank.png`
  - mostra o espectro de candidatos positivos reais ordenado por `k0L`.
  - destaca em vermelho quais candidatos foram usados no matching da Tabela 8.

- `img/magnitude/helmvec2_rect_modeXX_candYY_Et_magnitude.png`
  - mapa espacial da magnitude do campo transversal `Et` do modo casado.
  - o contorno do dominio aparece desenhado por padrao.
  - a interface horizontal entre materiais tambem aparece desenhada.

- `img/quiver/helmvec2_rect_modeXX_candYY_Et_quiver.png`
  - diagrama de setas do campo transversal `Et`.
  - o contorno do dominio aparece desenhado por padrao.
  - a interface horizontal entre materiais tambem aparece desenhada.

- `img/scalar/helmvec2_rect_modeXX_candYY_Ez_scalar.png`
  - mapa escalar da componente longitudinal `Ez`.
  - o contorno do dominio aparece desenhado por padrao.
  - a interface horizontal entre materiais tambem aparece desenhada.

## Como ler essas figuras

### 1) Figura 11 / Tabela 8

Este grafico responde:

- quao perto o valor FEM ficou das referencias publicadas;
- se a sequencia modal casada pelo codigo acompanha a ordenacao do artigo.

### 2) Erro por modo

Este grafico responde:

- quais linhas da Tabela 8 estao mais proximas da referencia `HELMVEC2`;
- onde a diferenca para `Hayata` fica maior ou menor.

### 3) Espectro de candidatos

Este grafico responde:

- quantos candidatos fisicos apareceram antes do matching final;
- em que posicoes do espectro estavam as raizes escolhidas para a Tabela 8.

Em termos didaticos, ele mostra que a tabela final nao nasce diretamente do
solver. Primeiro o codigo produz um espectro de candidatos, depois filtra e
finalmente casa cada referencia com um candidato disponivel.

### 4) Magnitude e quiver de `Et`

Essas figuras respondem:

- como a parte transversal `Et` se distribui espacialmente;
- em que regiao o modo concentra mais magnitude;
- qual a orientacao local dominante do campo transversal.

### 5) Mapa escalar de `Ez`

Essa figura responde:

- como a parte longitudinal `Ez` se distribui no dominio;
- em que regioes a componente longitudinal muda de sinal;
- como o mesmo modo casado combina `Et` e `Ez` dentro do vetor de estado.

## Relacao com a documentacao teorica

Este material deve ser lido junto com:

- [HELMVEC2_CSV_Referencia.md](HELMVEC2_CSV_Referencia.md)
- [src/helmvec2/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec2/README.md)
- [2.2.3_Determinação do número de onda.md](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md)

## Observacao importante

O `HELMVEC2` agora exporta, para cada modo casado da Tabela 8:

- um CSV e um VTK de `Et`;
- um CSV e um VTK de `Ez`;
- um mapa de magnitude e um `quiver` para `Et`;
- um mapa escalar para `Ez`.

Mesmo assim, a leitura principal da familia continua sendo a da Secao 2.2.3:

- primeiro o espectro de candidatos;
- depois o filtro fisico;
- por fim o matching final com a Tabela 8.

O `ez_ratio` continua salvo no `candidates.csv` e no `modes.csv`, mas segue
sem figura propria porque nao estava ajudando a leitura principal desta
familia.

A malha triangular interna continua opcional e pode ser sobreposta com
`--show-mesh`; o contorno do dominio, por outro lado, aparece sempre nas
figuras espaciais. No caso retangular atual do `HELMVEC2`, a interface entre
os dois materiais (`eps_bottom=1.5` e `eps_top=1.0`) tambem aparece sempre
como uma linha horizontal tracejada em `y = L/2`.
