# `HELMVEC`: geração de imagens a partir dos CSVs

Este arquivo documenta o script:

- `python3 scripts/helmvec.py`

Ele foi criado para transformar as saídas tabulares do `HELMVEC` em figuras
didáticas, preservando a ligação:

`modes.csv -> fields_<modo>.csv -> imagem do modo`

## 1) Objetivo

Gerar, para cada modo exportado pelos casos:

- `rect`
- `circle`
- `coax`

as seguintes imagens:

- mapa de magnitude do campo transversal por célula
- diagrama de setas do campo transversal reconstruído
- gráfico-resumo de erro por modo
- gráfico-resumo de `rho` por modo

## 2) Fontes de entrada

O script lê, por caso:

- `out/helmvec/<caso>/csv/edge_<caso>_modes.csv`
- `out/helmvec/<caso>/csv/edge_<caso>_fields_<...>.csv`
- `out/helmvec/<caso>/vtk/*.vtk`

Observação importante:

- os CSVs são a fonte principal dos valores físicos e dos metadados;
- o VTK é usado para reaproveitar a conectividade triangular real da malha;
- o campo vetorial tabular é lido do CSV por célula, enquanto a geometria da
  malha vem do VTK do mesmo modo.

## 3) Estrutura de saída

Para cada caso, o script escreve em:

- `out/helmvec/<caso>/img/magnitude/`
- `out/helmvec/<caso>/img/quiver/`
- `out/helmvec/<caso>/img/`

Nomes típicos:

- `out/helmvec/rect/img/magnitude/helmvec_rect_TE_m1_n0_rank01_magnitude.png`
- `out/helmvec/rect/img/quiver/helmvec_rect_TE_m1_n0_rank01_quiver.png`
- `out/helmvec/rect/img/helmvec_rect_error_by_mode.png`
- `out/helmvec/rect/img/helmvec_rect_rho_by_mode.png`

## 4) Uso básico

Gerar todas as imagens dos três casos:

```bash
python3 scripts/helmvec.py
```

Gerar apenas o caso retangular:

```bash
python3 scripts/helmvec.py --case rect
```

Gerar circular e coaxial com a malha sobreposta:

```bash
python3 scripts/helmvec.py --case circle --case coax --show-mesh
```

Ajustar resolução e densidade máxima de setas:

```bash
python3 scripts/helmvec.py --dpi 240 --max-arrows 400
```

## 5) Parâmetros do script

- `--root`
  - diretório raiz das saídas do `HELMVEC`
  - padrão: `out/helmvec`
- `--case rect|circle|coax|all`
  - seleciona quais casos serão processados
  - pode ser repetido
  - padrão: `all`
- `--dpi`
  - resolução das imagens salvas
- `--show-mesh`
  - sobrepõe a malha triangular nas figuras
- `--max-arrows`
  - limita aproximadamente o número de setas do quiver
  - o script faz subamostragem dos centróides quando necessário

## 6) Relação com os CSVs

O vínculo entre a linha modal e o arquivo por célula é:

- coluna `fields_csv_file` em `edge_<caso>_modes.csv`

Além disso, a própria linha modal também aponta para:

- coluna `vtk_file`, usada para recuperar a conectividade triangular da malha

Assim:

- `modes.csv` concentra os metadados do modo
- `fields_<modo>.csv` concentra o campo reconstruído por célula
- a imagem usa ambos

## 7) O que cada figura representa

Mapa de magnitude:

- usa `Emag` no ramo `TE` e `Hmag` no ramo `TM`
- colore cada triângulo pela magnitude do campo transversal reconstruído
- preserva a conectividade triangular do VTK do mesmo modo

Diagrama de setas:

- usa `Ex,Ey` no ramo `TE`
- usa `Hx,Hy` no ramo `TM`
- desenha o vetor transversal no centroide de cada triângulo selecionado
- a orientação das setas vem diretamente do campo reconstruído no pós-processamento

Gráfico de erro por modo:

- usa `error_percent` do `modes.csv`
- na imagem-resumo, o valor é mostrado em módulo:
  `|error_percent|`
- dentro de cada família modal, os pontos são ordenados por `kc_fem`
- em casos degenerados, o desempate usa `positive_rank`

Gráfico de `rho` por modo:

- usa `rho_abs` do `modes.csv`
- como o CSV já salva o valor em módulo, a figura mostra diretamente `|rho|`
- dentro de cada família modal, os pontos são ordenados por `kc_fem`
- em casos degenerados, o desempate usa `positive_rank`
- o intervalo vertical é fixado em `[0, 1.05]` para facilitar comparação visual

## 8) Observação didática

O script foi pensado como ferramenta de leitura do problema em C++ e não apenas
como pós-processamento gráfico. Por isso, a navegação recomendada é:

1. abrir `edge_<caso>_modes.csv`
2. escolher uma linha modal
3. abrir o `fields_<modo>.csv` indicado em `fields_csv_file`
4. conferir a imagem correspondente em `img/magnitude/` e `img/quiver/`

Esse fluxo ajuda a rastrear com clareza:

- o modo identificado
- seus parâmetros numéricos
- os valores por célula exportados
- a figura final usada para inspeção visual
