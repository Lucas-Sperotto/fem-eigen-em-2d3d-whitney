# HELMVEC: Referencia dos Campos do CSV de Modos

Este arquivo documenta o significado didatico das colunas do CSV de modos
gerado pelos executaveis da familia `HELMVEC`:

- `edge_rect`
- `edge_circle`
- `edge_coax`

Ele deve ser lido junto com:

- a traducao da Secao 2.2.1 do artigo em
  [docs/traducao/2.2.1_Guia de Onda Campos Vetoriais.md](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md)
- a documentacao de implementacao em
  [src/helmvec/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec/README.md)

## Cabecalho-base do caso retangular

O CSV de modos do caso retangular usa hoje o cabecalho:

```text
family,transverse_label,mode_label,positive_rank,eig_index,m,n,ar_m,b_m,kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,error_percent,rho_abs,field_status,fields_csv_file,vtk_file
```

Nos casos circular e coaxial, a ideia fisica e a mesma, mas os campos
geometricos mudam:

- circular: aparece `r_m` e a forma normalizada `kc_r_fem`, `kc_r_ana`
- coaxial: aparecem `r1_m,r2_m` e a forma normalizada `kc_r1_fem`, `kc_r1_ana`
- circular e coaxial usam `p` no lugar de `n`

## Significado de cada coluna

- `family`
  - familia modal identificada no problema vetorial.
  - no uso atual do projeto, os valores tipicos sao `TE` e `TM`.

- `transverse_label`
  - campo transversal associado ao ramo calculado.
  - `TE -> Et`
  - `TM -> Ht`

- `mode_label`
  - rotulo textual do modo identificado.
  - exemplos: `TE_m1_n0`, `TM_m1_n1`, `TE_m2_p1`.

- `positive_rank`
  - ordem didatica do modo dentro do conjunto realmente exportado da familia.

- `eig_index`
  - posicao do autovalor/autovetor no espectro retornado pelo solver.

- `m`
  - primeiro indice modal associado ao modo identificado por matching com a
    referencia analitica.

- `n`
  - segundo indice modal do caso retangular.
  - nos casos circular e coaxial, esse papel e exercido por `p`.

- `ar_m`
  - dimensao `a` do guia retangular, em metros.

- `b_m`
  - dimensao `b` do guia retangular, em metros.

- `r_m`
  - raio do guia circular, em metros.

- `r1_m`
  - raio interno da linha coaxial, em metros.

- `r2_m`
  - raio externo da linha coaxial, em metros.

- `kc_fem`
  - numero de onda de corte calculado numericamente pelo FEM.
  - como o autovalor do problema de aresta e `lambda = kc^2`, este valor vem
    de `sqrt(lambda)`.

- `kc_ana`
  - numero de onda de corte analitico associado ao modo identificado.
  - serve como referencia para comparacao e erro.

- `kc_ar_fem`
  - forma adimensional do cutoff numerico no caso retangular.
  - no projeto atual, `kc_ar_fem = kc_fem * a`.

- `kc_ar_ana`
  - forma adimensional do cutoff analitico no caso retangular.
  - no projeto atual, `kc_ar_ana = kc_ana * a`.

- `kc_r_fem`
  - forma adimensional do cutoff numerico no caso circular.
  - no projeto atual, `kc_r_fem = kc_fem * r`.

- `kc_r_ana`
  - forma adimensional do cutoff analitico no caso circular.

- `kc_r1_fem`
  - forma adimensional do cutoff numerico no caso coaxial usando `r1`.
  - no projeto atual, `kc_r1_fem = kc_fem * r1`.

- `kc_r1_ana`
  - forma adimensional do cutoff analitico no caso coaxial usando `r1`.

- `error_percent`
  - erro percentual entre `kc_fem` e `kc_ana`.
  - usado como indicador numerico de proximidade com a referencia analitica.

- `rho_abs`
  - valor absoluto da correlacao usada no casamento modal.
  - mede o quao bem o autovetor numerico de aresta combina com o modo
    analitico candidato.
  - valores mais proximos de `1` indicam casamento mais limpo.

- `field_status`
  - estado textual da reconstrucao exportada.
  - no estado atual do projeto, indica que o campo foi reconstruido nos
    centroides e normalizado por pico unitario para visualizacao.

- `fields_csv_file`
  - nome do arquivo CSV de campos por celula associado a esse modo.

- `vtk_file`
  - nome do arquivo VTK do mesmo modo.
  - o VTK preserva a conectividade triangular e o campo vetorial em
    `CELL_DATA`.

## Como `rho_abs` e calculado e por que funciona

No `HELMVEC`, o casamento modal tambem nao compara apenas autovalores. Ele
compara a forma espacial do autovetor FEM de aresta com a forma espacial de um
candidato analitico, ambos escritos sobre os mesmos graus de liberdade de
aresta.

No codigo, a correlacao usa a matriz global de massa `T` da formulacao
vetorial transversal:

- Eq. `(67)`: massa vetorial local
- Eq. `(65)`: problema global `S e = kc^2 T e`

Essa mesma matriz `T` define o produto interno discreto usado no matching:

```text
<x, y>_T = x^T T y
||x||_T = sqrt(x^T T x)
rho     = |<x, y>_T| / (||x||_T ||y||_T)
```

No repositório, isso aparece diretamente em:

- [mode_match_rect_edge.hpp](/home/sperotto/tp3485-fem-eigen-em/src/edge/mode_match_rect_edge.hpp)
- [mode_match_circle_edge.hpp](/home/sperotto/tp3485-fem-eigen-em/src/edge/mode_match_circle_edge.hpp)
- [mode_match_coax_edge.hpp](/home/sperotto/tp3485-fem-eigen-em/src/edge/mode_match_coax_edge.hpp)

### Leitura matematica

- `x` e o autovetor numerico produzido pelo solve FEM.
- `y` e o modo analitico candidato, projetado nos mesmos DOFs de aresta.
- `x^T T y` aproxima a sobreposicao energetica dos dois campos no espaco
  vetorial discreto.
- `||x||_T` e `||y||_T` removem a dependencia da escala dos vetores.
- o valor absoluto remove a ambiguidade de sinal global do autovetor.

### Por que isso funciona

- A matriz de massa vetorial `T` representa, no espaco discreto de arestas, o
  produto interno natural do problema transversal.
- Portanto, `rho` funciona como um "cosseno" entre dois modos no espaco
  induzido por `T`.
- Se o autovetor FEM e o candidato analitico tiverem a mesma distribuicao
  espacial tangencial, `rho` fica proximo de `1`.
- Se as formas espaciais forem diferentes, cancelamentos no produto interno
  fazem `rho` diminuir.

Nos casos circular e coaxial, familias degeneradas podem produzir mais de um
candidato com o mesmo cutoff analitico. Nesses casos, o matching testa as
variantes compativeis e guarda a maior correlacao.

## Como ler junto com o CSV de campos

O CSV de modos responde:

- qual modo foi identificado
- qual e o seu cutoff
- qual foi o erro contra a referencia
- quao limpo foi o matching modal
- qual arquivo contem o campo por celula daquele modo

O CSV de campos responde:

- qual e o valor reconstruido do campo transversal no centro de cada triangulo
- quais sao as componentes `Ex,Ey` no ramo `TE`
- quais sao as componentes `Hx,Hy` no ramo `TM`
- qual e a magnitude correspondente (`Emag` ou `Hmag`)

Para o detalhe coluna por coluna do arquivo por celula, ver:

- [HELMVEC_CSV_Campos_Referencia.md](HELMVEC_CSV_Campos_Referencia.md)

## Referencias cruzadas

- teoria traduzida:
  [2.2.1 Guia de Onda Campos Vetoriais](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md)
- implementacao:
  [src/helmvec/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec/README.md)
- visao operacional dos executaveis:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
