# HELMVEC: Referencia dos Campos do CSV por Celula

Este arquivo documenta o significado didatico das colunas do CSV de campos por
celula gerado pelos executaveis da familia `HELMVEC`.

Ele deve ser lido junto com:

- [HELMVEC_CSV_Modos_Referencia.md](HELMVEC_CSV_Modos_Referencia.md)
- [src/helmvec/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec/README.md)

## Cabecalho atual

O CSV por modo usa hoje um cabecalho dependente do ramo reconstruido.

Para modos `TE`:

```text
cell_id,xc_m,yc_m,Ex,Ey,Emag
```

Para modos `TM`:

```text
cell_id,xc_m,yc_m,Hx,Hy,Hmag
```

O proprio nome do arquivo liga esse campo a um modo especifico, por exemplo:

- `edge_rect_fields_TE_m1_n0_rank01.csv`
- `edge_circle_fields_TM_m1_p1_rank02.csv`
- `edge_coax_fields_TE_m2_p1_rank03.csv`

## Significado de cada coluna

- `cell_id`
  - identificador do triangulo na numeracao interna da malha.
  - a ordem coincide com a ordem dos triangulos exportados no VTK do mesmo
    modo.

- `xc_m`
  - coordenada `x` do centroide do triangulo, em metros.

- `yc_m`
  - coordenada `y` do centroide do triangulo, em metros.

- `Ex`
  - componente `x` do campo eletrico transversal reconstruido no centroide.
  - aparece quando o ramo exportado e `TE`, isto e, quando o problema vetorial
    foi resolvido para `Et`.

- `Ey`
  - componente `y` do campo eletrico transversal reconstruido no centroide.
  - aparece junto com `Ex` no ramo `TE`.

- `Emag`
  - magnitude do vetor eletrico transversal reconstruido:

```text
Emag = sqrt(Ex^2 + Ey^2)
```

- `Hx`
  - componente `x` do campo magnetico transversal reconstruido no centroide.
  - aparece quando o ramo exportado e `TM`, isto e, quando o problema dual foi
    resolvido para `Ht`.

- `Hy`
  - componente `y` do campo magnetico transversal reconstruido no centroide.
  - aparece junto com `Hx` no ramo `TM`.

- `Hmag`
  - magnitude do vetor magnetico transversal reconstruido:

```text
Hmag = sqrt(Hx^2 + Hy^2)
```

## Leitura fisica

No estado atual do projeto:

- para o ramo `TE`, o campo reconstruido representa `Et`
- para o ramo `TM`, o campo reconstruido representa `Ht`

Por isso, os nomes das colunas agora seguem explicitamente o campo
reconstruido:

- `TE` usa `Ex`, `Ey`, `Emag`
- `TM` usa `Hx`, `Hy`, `Hmag`

Essa escolha evita ambiguidade quando o leitor abre um `fields_<modo>.csv`
diretamente, sem passar antes pelo `modes.csv`.

## Como o campo e obtido

O campo e reconstruido no centro de cada triangulo a partir dos coeficientes
de aresta do autovetor FEM:

```text
F_t(x_c) = sum_{m=1..3} e_m W_m(x_c)
```

onde:

- `e_m` sao os coeficientes do modo no espaco de arestas
- `W_m` sao as bases locais de Whitney do triangulo
- `x_c` e o centroide do elemento

No repositorio, isso aparece em:

- [edge_mode_post.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec/edge_mode_post.hpp)

## O que o quiver mostra

As imagens de `quiver` do `HELMVEC` usam exatamente essas colunas do CSV por
celula:

- no ramo `TE`, as setas sao desenhadas com `Ex` e `Ey`
- no ramo `TM`, as setas sao desenhadas com `Hx` e `Hy`

Isto significa que o `quiver` nao mostra gradientes de potencial escalar, como
acontece no `HELM10`. Aqui ele mostra o proprio campo vetorial transversal
reconstruido no centroide de cada triangulo:

```text
F_t(x_c) = sum_{m=1..3} e_m W_m(x_c)
```

onde `F_t` deve ser lido como:

- `E_t` no ramo `TE`
- `H_t` no ramo `TM`

## Observacao importante sobre normalizacao

No estado atual do `HELMVEC`, o campo salvo neste CSV ja esta normalizado por
amplitude maxima unitária dentro do modo exportado.

Isso significa:

- a forma espacial do campo e preservada
- a comparacao visual entre triangulos do mesmo modo fica direta
- a escala absoluta do autovetor nao e usada, o que e coerente com o fato de
  autovetores de problemas de autovalor terem escala arbitraria

Essa politica aparece no `modes.csv` pela coluna `field_status`.

## Leitura conjunta com o VTK

O `CSV` por celula e a versao tabular minima do campo reconstruido.

O `VTK` do mesmo modo preserva:

- os nos da malha
- a conectividade triangular
- o vetor por celula em `CELL_DATA`

Em outras palavras:

- `fields_<modo>.csv` e o artefato mais didatico para inspeção tabular
- `vtk_file` e o artefato mais natural para malha + visualizacao geometrica

## Referencias cruzadas

- glossario do modo:
  [HELMVEC_CSV_Modos_Referencia.md](HELMVEC_CSV_Modos_Referencia.md)
- guia do script de imagens:
  [HELMVEC_Imagens_Referencia.md](HELMVEC_Imagens_Referencia.md)
