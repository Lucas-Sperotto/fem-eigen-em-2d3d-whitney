# `HELM10`: geração de imagens a partir dos CSVs

Este arquivo documenta o script:

- `python3 scripts/helm10.py`

Ele foi criado para transformar as saídas tabulares do `HELM10` em figuras
didáticas, preservando a ligação:

`modes.csv -> fields_<modo>.csv -> imagem do modo`

## 1) Objetivo

Gerar, para cada modo exportado pelos casos:

- `rect`
- `circle`
- `coax`

as seguintes imagens:

- diagrama de isolinhas do potencial escalar `psi`
- diagrama de setas do campo transversal `(Ex, Ey)`
- gráfico-resumo de erro por modo
- gráfico-resumo de `rho` por modo

## 2) Fontes de entrada

O script lê, por caso:

- `out/helm10/<caso>/csv/helm10_<caso>_modes.csv`
- `out/helm10/<caso>/csv/helm10_<caso>_fields_<...>.csv`

e, quando disponível, reutiliza também:

- `out/helm10/<caso>/vtk/*.vtk`

Observação importante:

- os CSVs são a fonte principal dos valores físicos e dos metadados;
- o VTK é usado apenas para reaproveitar a conectividade triangular real da
  malha, especialmente útil no caso coaxial.

## 3) Estrutura de saída

Para cada caso, o script escreve em:

- `out/helm10/<caso>/img/isopotential/`
- `out/helm10/<caso>/img/quiver/`
- `out/helm10/<caso>/img/`

Nomes típicos:

- `out/helm10/rect/img/isopotential/helm10_rect_TE_m1_n0_rank01_psi.png`
- `out/helm10/rect/img/quiver/helm10_rect_TE_m1_n0_rank01_field.png`
- `out/helm10/rect/img/quiver/helm10_rect_TE_m1_n0_rank01_gradient.png`
- `out/helm10/rect/img/helm10_rect_error_by_mode.png`
- `out/helm10/rect/img/helm10_rect_rho_by_mode.png`

## 4) Uso básico

Gerar todas as imagens dos três casos:

```bash
python3 scripts/helm10.py
```

Nesse modo padrão, o script salva os dois quivers de cada modo:

- `..._field.png`
- `..._gradient.png`

Gerar apenas o caso retangular:

```bash
python3 scripts/helm10.py --case rect
```

Gerar circular e coaxial com malha sobreposta:

```bash
python3 scripts/helm10.py --case circle --case coax --show-mesh
```

Ajustar resolução e densidade das setas:

```bash
python3 scripts/helm10.py --dpi 240 --quiver-nx 20 --quiver-ny 10
```

Gerar o `quiver` com as derivadas do potencial, em vez de `Ex,Ey`:

```bash
python3 scripts/helm10.py --quiver-source gradient
```

## 5) Parâmetros do script

- `--root`
  - diretório raiz das saídas do `HELM10`
  - padrão: `out/helm10`
- `--case rect|circle|coax|all`
  - seleciona quais casos serão processados
  - pode ser repetido
  - padrão: `all`
- `--dpi`
  - resolução das imagens salvas
- `--levels`
  - número de níveis de isolinhas para `psi`
- `--quiver-nx`
  - número de pontos na direção `x` da grade regular usada no quiver
- `--quiver-ny`
  - número de pontos na direção `y` da grade regular usada no quiver
- `--quiver-source`
  - `all`: salva `field` e `gradient` na mesma execução
  - `field`: desenha setas com `Ex` e `Ey`
  - `gradient`: desenha setas com `dpsi_dx` e `dpsi_dy`
- `--show-mesh`
  - sobrepõe a malha triangular na figura
- `--no-vtk`
  - força a triangulação apenas a partir dos nós dos CSVs

Observação geométrica atual do quiver:

- `rect` usa grade regular no plano
- `circle` usa uma amostragem polar didática com `13` círculos concêntricos e
  `51` pontos em cada círculo
- `coax` usa uma amostragem polar didática com `13` círculos concêntricos e
  `51` pontos em cada círculo

## 6) Relação com os CSVs

O vínculo entre a linha modal e o arquivo nodal é:

- coluna `fields_csv_file` em `helm10_<caso>_modes.csv`

Assim:

- `modes.csv` concentra os metadados do modo
- `fields_<modo>.csv` concentra os valores nodais
- a imagem usa ambos

## 7) O que cada figura representa

Diagrama de isolinhas:

- usa `psi` do `fields_<modo>.csv`
- mostra o potencial longitudinal do modo
- em TE, esse potencial está associado a `Hz`
- em TM, esse potencial está associado a `Ez`

Diagrama de setas:

- pode usar `Ex` e `Ey` do `fields_<modo>.csv`
- ou `dpsi_dx` e `dpsi_dy`, dependendo da fonte escolhida
- no projeto, `Ex` e `Ey` são sempre exportados sem multiplicação por `Ztm`
- isso vale tanto para TE quanto para TM
- para visualização, os vetores do quiver são normalizados antes do desenho

Gráfico de erro por modo:

- usa `error_percent` do `modes.csv`
- na imagem-resumo, o valor é mostrado em módulo:
  `|error_percent|`
- dentro de cada familia modal, os pontos sao ordenados por `kc_fem`
- em casos degenerados, o desempate usa `positive_rank`

Gráfico de `rho` por modo:

- usa `rho_abs` do `modes.csv`
- como o CSV ja salva o valor em modulo, a figura mostra diretamente `|rho|`
- dentro de cada familia modal, os pontos sao ordenados por `kc_fem`
- em casos degenerados, o desempate usa `positive_rank`
- o intervalo vertical e fixado em `[0, 1.05]` para facilitar comparacao visual

## 8) Observação didática

O script foi pensado como ferramenta de leitura do problema em C++ e não apenas
como pós-processamento gráfico. Por isso, a navegação recomendada é:

1. abrir `helm10_<caso>_modes.csv`
2. escolher uma linha modal
3. abrir o `fields_<modo>.csv` indicado em `fields_csv_file`
4. conferir a imagem correspondente em `img/isopotential/` e `img/quiver/`

Esse fluxo ajuda a rastrear com clareza:

- o modo identificado
- seus parâmetros numéricos
- os valores nodais exportados
- a figura final usada para inspeção visual
