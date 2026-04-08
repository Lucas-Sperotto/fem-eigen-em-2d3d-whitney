# HELMVEC1: Referencia dos CSVs de Campos

Este arquivo documenta os CSVs espaciais produzidos pela familia `HELMVEC1`.

Ele deve ser lido junto com:

- [HELMVEC1_CSV_Modos_Referencia.md](HELMVEC1_CSV_Modos_Referencia.md)
- [HELMVEC1_Imagens_Referencia.md](HELMVEC1_Imagens_Referencia.md)
- [src/helmvec1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/README.md)

## Ideia principal

No `HELMVEC1`, o modo exportado nao e sempre do mesmo tipo fisico.

Pela Eq. `(92)`, o autovetor misto e lido como:

```text
x = [x_edge ; x_scalar]
```

Depois da classificacao por energia dominante:

- modos com `dominant_block = edge` sao exportados como campo vetorial
  transversal
- modos com `dominant_block = scalar` sao exportados como componente
  longitudinal escalar

Por isso, o `fields_<modo>.csv` do `HELMVEC1` tem dois formatos validos.

## Como localizar o arquivo certo

Cada linha do `modes.csv` aponta para seu arquivo espacial correspondente por:

- `field_data_kind`
- `fields_csv_file`
- `vtk_file`

Valores atuais de `field_data_kind`:

- `edge_vector_cell`
- `scalar_nodal`

## 1) CSV vetorial por celula

Quando `field_data_kind = edge_vector_cell`, o modo foi dominado pelo bloco de
aresta e o CSV representa o campo transversal reconstruido por triangulo.

### Cabecalho

No ramo `Et`:

```text
cell_id,xc_m,yc_m,Ex,Ey,Emag
```

No ramo `Ht`:

```text
cell_id,xc_m,yc_m,Hx,Hy,Hmag
```

### Significado

- `cell_id`
  - indice do triangulo na malha.

- `xc_m`, `yc_m`
  - coordenadas do centroide do triangulo, em metros.

- `Ex`, `Ey`
  - componentes cartesianas do campo transversal eletrico reconstruido.
  - aparecem quando `component_label = Et`.

- `Hx`, `Hy`
  - componentes cartesianas do campo transversal magnetico reconstruido.
  - aparecem quando `component_label = Ht`.

- `Emag` ou `Hmag`
  - magnitude euclidiana do vetor transversal por celula.

### Como esse campo e reconstruido

No bloco de aresta, o repositorio reaproveita a mesma leitura didatica do
`HELMVEC`:

```text
F_t(x_c) = sum_m e_m W_m(x_c)
```

onde:

- `e_m` sao os coeficientes do autovetor no bloco `edge`
- `W_m` sao as funcoes de Whitney locais
- `x_c` e o centroide do triangulo

No codigo, isso passa por:

- [mixed_field_output.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/mixed_field_output.hpp)
- [edge_mode_post.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec/edge_mode_post.hpp)

## 2) CSV escalar nodal

Quando `field_data_kind = scalar_nodal`, o modo foi dominado pelo bloco escalar
e o CSV representa a componente longitudinal nodal.

### Cabecalho

No ramo `Ez`:

```text
node_id,x_m,y_m,Ez
```

No ramo `Hz`:

```text
node_id,x_m,y_m,Hz
```

### Significado

- `node_id`
  - indice do no na malha.

- `x_m`, `y_m`
  - coordenadas do no, em metros.

- `Ez` ou `Hz`
  - valor nodal da componente longitudinal dominante do modo.

### Como esse campo e reconstruido

No bloco escalar, o repositorio extrai apenas a parte escalar do autovetor
misto e reinsere os graus de liberdade eliminados pela condicao de contorno:

```text
u_full(node) =
    autovetor(dof(node)), se o no tiver DOF ativo
    0,                   se o no tiver sido eliminado
```

O vetor nodal e normalizado por pico unitario antes da exportacao espacial.

No codigo, isso passa por:

- [mixed_field_output.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helmvec1/mixed_field_output.hpp)
- [scalar_mode_post.hpp](/home/sperotto/tp3485-fem-eigen-em/src/helm10/scalar_mode_post.hpp)

## 3) VTK associado

Cada `fields_<modo>.csv` tambem tem um `vtk_file` correspondente:

- modos `edge` -> VTK vetorial por celula
- modos `scalar` -> VTK escalar nodal

O script `scripts/helmvec1.py` usa esse `vtk_file` para recuperar a
conectividade da malha e gerar as imagens espaciais.
