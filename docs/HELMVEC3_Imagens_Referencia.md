# HELMVEC3: Referencia das Imagens

Este arquivo documenta o script:

- [scripts/helmvec3.py](/home/sperotto/tp3485-fem-eigen-em/scripts/helmvec3.py)

Ele le os CSVs produzidos hoje pelos executaveis:

- `helmvec3_fig12_rect`
- `helmvec3_fig13_rect`

e gera imagens da Secao 2.2.4 a partir da Eq. `(136)`.

## Ideia principal

No `HELMVEC3`, o foco didatico nao e um mapa espacial por modo. O que importa
primeiro e:

- a curva da Figura 12 / Tabela 9
- o preview de continuidade modal da Figura 13
- a familia completa de curvas da Tabela 10

Por isso, o script trabalha sobre:

- `out/helmvec3/fig12_rect/csv/helmvec3_fig12_rect_table9.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_preview.csv`
- `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv`
- `out/helmvec3/fig12_rect/csv/*_Et_fields.csv`
- `out/helmvec3/fig12_rect/csv/*_Ez_fields.csv`
- `out/helmvec3/fig13_rect/csv/*_Et_fields.csv`
- `out/helmvec3/fig13_rect/csv/*_Ez_fields.csv`
- `out/helmvec3/fig12_rect/vtk/*_Et.vtk`
- `out/helmvec3/fig12_rect/vtk/*_Ez.vtk`
- `out/helmvec3/fig13_rect/vtk/*_Et.vtk`
- `out/helmvec3/fig13_rect/vtk/*_Ez.vtk`
- `out/helmvec3/fig12_rect/run_timing.csv`
- `out/helmvec3/fig13_rect/run_timing.csv`

e salva as imagens em:

- `out/helmvec3/fig12_rect/img/`
- `out/helmvec3/fig13_rect/img/`

## Como executar

```bash
python3 scripts/helmvec3.py
python3 scripts/helmvec3.py --case fig12_rect
python3 scripts/helmvec3.py --case fig13_rect
python3 scripts/helmvec3.py --dpi 120
python3 scripts/helmvec3.py --case fig13_rect --dpi 120 --show-mesh
```

## Imagens geradas

- `helmvec3_fig12_rect_table9_beta_over_k0.png`
  - compara `beta/k0(FEM)` com as referencias analitica e `HELMVEC3` da
    Figura 12 / Tabela 9.

- `helmvec3_fig12_rect_table9_error_by_point.png`
  - mostra o erro relativo absoluto por ponto contra as duas referencias da
    Tabela 9.

- `helmvec3_fig13_rect_preview_branch.png`
  - mostra o ramo rastreado por continuidade para o `d/a` escolhido na linha
    de comando.

- `helmvec3_fig13_rect_table10_fem_branches.png`
  - mostra, no mesmo plano, as curvas FEM da Tabela 10 para cada bloco `d/a`.

- `helmvec3_fig13_rect_table10_error_by_branch.png`
  - mostra o erro relativo absoluto contra a referencia analitica para cada
    bloco `d/a`.

- `img/magnitude/*_Et_magnitude.png`
  - mostra a magnitude do campo transversal `Et` para cada ponto exportado.
  - os arquivos ficam organizados em subpastas por `d/a`.

- `img/quiver/*_Et_quiver.png`
  - mostra o diagrama de setas de `Et` para cada ponto exportado.
  - os arquivos ficam organizados em subpastas por `d/a`.

- `img/scalar/*_Ez_scalar.png`
  - mostra o mapa escalar da componente longitudinal `Ez`.
  - os arquivos ficam organizados em subpastas por `d/a`.

## Como ler essas figuras

### 1) Figura 12 / Tabela 9

Este grafico responde:

- se a curva FEM acompanha a referencia analitica;
- quao perto a implementacao ficou do valor publicado pelo `HELMVEC3`.

### 2) Erro da Tabela 9

Este grafico responde:

- em quais pontos da Figura 12 a discrepancia ficou maior;
- como o erro varia ao longo da curva.

### 3) Preview de ramo

Este grafico responde:

- como um ramo continuo de `beta/k0` se comporta para um `d/a` escolhido;
- se a continuidade modal parece plausivel antes da comparacao completa da
  Tabela 10.

### 4) Curvas da Tabela 10

Este grafico responde:

- como os ramos FEM mudam com `d/a`;
- como comparar rapidamente os blocos de dispersao da Figura 13.

### 5) Erro da Tabela 10

Este grafico responde:

- quais blocos `d/a` ficaram mais proximos da referencia analitica;
- como o erro varia ao longo de cada ramo.

### 6) Magnitude e quiver de `Et`

Essas figuras respondem:

- como a parte transversal `Et` se distribui no dominio;
- em que regiao o ponto exportado concentra maior magnitude;
- qual a orientacao local dominante do campo transversal.

### 7) Mapa escalar de `Ez`

Essas figuras respondem:

- como a parte longitudinal `Ez` se distribui no dominio;
- onde a componente longitudinal muda de sinal;
- como o mesmo ponto exportado combina `Et` e `Ez`.

## Relacao com a documentacao teorica

Este material deve ser lido junto com:

- [HELMVEC3_CSV_Referencia.md](HELMVEC3_CSV_Referencia.md)
- [src/helmvec3/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helmvec3/README.md)
- [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)

## Observacao importante

A separacao do `HELMVEC3` em dois executaveis publicos mudou a origem
dos CSVs e VTKs:

- `helmvec3_fig12_rect` produz a parte da Figura 12 / Tabela 9;
- `helmvec3_fig13_rect` produz o preview e a Tabela 10.

O `scripts/helmvec3.py` ja foi adaptado para essa nova casca publica. A leitura
conceitual das imagens, porem, continua a mesma.

O `HELMVEC3` agora exporta, para cada ponto publicado da Figura 12, do preview
e da Tabela 10:

- um CSV e um VTK de `Et`;
- um CSV e um VTK de `Ez`;
- um mapa de magnitude e um `quiver` para `Et`;
- um mapa escalar para `Ez`.

Mesmo assim, a leitura principal da familia continua sendo a da Secao 2.2.4:

- a curva da Figura 12 / Tabela 9;
- o preview de continuidade modal;
- a validacao completa da Figura 13 / Tabela 10.

Nas figuras espaciais, o contorno do dominio aparece sempre. A interface de
material tambem aparece sempre:

- horizontal em `y = d12` no caso da Figura 12;
- vertical em `x = d` no preview e na Tabela 10.

Para facilitar a navegacao, os diretorios espaciais ficam organizados em
subpastas por `d/a`, por exemplo:

- `img/magnitude/figure12_da_0_225/`
- `img/quiver/preview_da_0_2/`
- `img/scalar/table10_da_0_5/`
