# Índice Geral de Documentação — NASA TP-3485

Mapa completo de todos os arquivos `.md` do repositório, organizados por tema.
Para a entrada mais rápida, vá ao [README raiz](../README.md).

---

## Resultados Numéricos

Páginas individuais por caso, com problema, condições de cálculo, tabelas de resultados e figuras de campos:

| # | Arquivo | Conteúdo |
| --- | --- | --- |
| 01 | [caso_01_tab1_helm10_rect.md](results/caso_01_tab1_helm10_rect.md) | Guia retangular escalar — Tabela 1 |
| 02 | [caso_02_tab2_helm10_circle.md](results/caso_02_tab2_helm10_circle.md) | Guia circular escalar — Tabela 2 |
| 03 | [caso_03_tab3_helm10_coax.md](results/caso_03_tab3_helm10_coax.md) | Linha coaxial escalar — Tabela 3 |
| 04 | [caso_04_tab4_helmvec_rect.md](results/caso_04_tab4_helmvec_rect.md) | Guia retangular edge 2D — Tabela 4 |
| 05 | [caso_05_tab5_helmvec_circle.md](results/caso_05_tab5_helmvec_circle.md) | Guia circular edge 2D — Tabela 5 |
| 06 | [caso_06_tab6_helmvec1_rect.md](results/caso_06_tab6_helmvec1_rect.md) | Retangular misto 3 comp. — Tabela 6 |
| 07 | [caso_07_tab7_helmvec1_circle.md](results/caso_07_tab7_helmvec1_circle.md) | Circular misto 3 comp. — Tabela 7 |
| 08 | [caso_08_fig11_tab8_helmvec2.md](results/caso_08_fig11_tab8_helmvec2.md) | Guia parcial `k0(β)` — Fig. 11 / Tabela 8 |
| 09 | [caso_09_fig12_tab9_helmvec3.md](results/caso_09_fig12_tab9_helmvec3.md) | Guia parcial `β(k0)`, ex. 1 — Fig. 12 / Tabela 9 |
| 10 | [caso_10_fig13_tab10_helmvec3.md](results/caso_10_fig13_tab10_helmvec3.md) | Guia parcial `β(k0)`, ex. 2 — Fig. 13 / Tabela 10 |
| 11 | [caso_11_tab12_fem3d_air.md](results/caso_11_tab12_fem3d_air.md) | Cavidade 3D ar — Tabela 12 |
| 12 | [caso_12_tab13_fem3d_half.md](results/caso_12_tab13_fem3d_half.md) | Cavidade 3D semi-preenchida — Tabela 13 |
| 13 | [caso_13_tab14_fem3d_cyl.md](results/caso_13_tab14_fem3d_cyl.md) | Cavidade cilíndrica 3D — Tabela 14 |
| 14 | [caso_14_tab15_fem3d_sphere.md](results/caso_14_tab15_fem3d_sphere.md) | Cavidade esférica 3D — Tabela 15 |
| — | [fem_vs_efgmi.md](results/fem_vs_efgmi.md) | Comparação consolidada FEM × EFGMI |
| — | [results/README.md](results/README.md) | Índice de resultados por seção |

---

## Teoria — Tradução Comentada do Artigo

| Arquivo | Conteúdo |
| --- | --- |
| [traducao/0.0_Capa.md](traducao/0.0_Capa.md) | Capa — metadados do NASA TP-3485 |
| [traducao/0.1_Sumario.md](traducao/0.1_Sumario.md) | Sumário |
| [traducao/0.2_Tabelas.md](traducao/0.2_Tabelas.md) | Lista de tabelas |
| [traducao/0.3_Figuras.md](traducao/0.3_Figuras.md) | Lista de figuras |
| [traducao/0.4_Simbolos.md](traducao/0.4_Simbolos.md) | Notação e símbolos (140+ entradas) |
| [traducao/0.5_Abstract.md](traducao/0.5_Abstract.md) | Resumo |
| [traducao/1_Introducao.md](traducao/1_Introducao.md) | Introdução — modos espúrios, motivação FEM |
| [traducao/2_Problemas_Bidimensionais.md](traducao/2_Problemas_Bidimensionais.md) | Visão geral dos problemas 2D |
| [traducao/2.1_Guias de Onda Homogêneos.md](<traducao/2.1_Guias de Onda Homogêneos.md>) | Sec. 2.1 — Formulação escalar, Eq. (43) |
| [traducao/2.2_Guia_Nao_Homogeneo.md](traducao/2.2_Guia_Nao_Homogeneo.md) | Sec. 2.2 — Guias não homogêneos (intro) |
| [traducao/2.2.1_Guia de Onda Campos Vetoriais.md](<traducao/2.2.1_Guia de Onda Campos Vetoriais.md>) | Sec. 2.2.1 — Edge 2D transversal, Eq. (65) |
| [traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md) | Sec. 2.2.2 — Sistema misto, Eq. (92) |
| [traducao/2.2.3_Determinação do número de onda.md](<traducao/2.2.3_Determinação do número de onda.md>) | Sec. 2.2.3 — k0 dado β, Eq. (119) |
| [traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md) | Sec. 2.2.4 — β dado k0, Eq. (136) |
| [traducao/3_Problemas_Tridimensionais.md](traducao/3_Problemas_Tridimensionais.md) | Sec. 3.1 — Cavidades 3D, Eq. (178) |
| [traducao/4_Consideracoes_Finais_Apendice_Referencias.md](traducao/4_Consideracoes_Finais_Apendice_Referencias.md) | Conclusões, apêndice FORTRAN, referências |
| [traducao/5_Report.md](traducao/5_Report.md) | Página de relatório técnico |

---

## Referência de Código

| Arquivo | Conteúdo |
| --- | --- |
| [Rastreabilidade_Equacoes_Artigo_Codigo.md](Rastreabilidade_Equacoes_Artigo_Codigo.md) | Mapeamento equação → arquivo C++ → função |
| [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md) | Todos os 21 executáveis: entradas, saídas, CLIs |
| [Casos_de_Teste_do_Artigo.md](Casos_de_Teste_do_Artigo.md) | 14 casos publicados vs casos extras do repositório |
| [Matriz_Casos_Executaveis_Arvore_de_Chamada.md](Matriz_Casos_Executaveis_Arvore_de_Chamada.md) | Árvore de chamada por executável |
| [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md) | Formato dos CSVs em `linop/` (matrizes CRS, autovalores) |

---

## Referência de Artefatos por Família

| Família | Modos CSV | Campos CSV | Imagens |
| --- | --- | --- | --- |
| HELM10 | [HELM10_CSV_Modos_Referencia.md](HELM10_CSV_Modos_Referencia.md) | [HELM10_CSV_Campos_Referencia.md](HELM10_CSV_Campos_Referencia.md) | [HELM10_Imagens_Referencia.md](HELM10_Imagens_Referencia.md) |
| HELMVEC | [HELMVEC_CSV_Modos_Referencia.md](HELMVEC_CSV_Modos_Referencia.md) | [HELMVEC_CSV_Campos_Referencia.md](HELMVEC_CSV_Campos_Referencia.md) | [HELMVEC_Imagens_Referencia.md](HELMVEC_Imagens_Referencia.md) |
| HELMVEC1 | [HELMVEC1_CSV_Modos_Referencia.md](HELMVEC1_CSV_Modos_Referencia.md) | [HELMVEC1_CSV_Campos_Referencia.md](HELMVEC1_CSV_Campos_Referencia.md) | [HELMVEC1_Imagens_Referencia.md](HELMVEC1_Imagens_Referencia.md) |
| HELMVEC2 | [HELMVEC2_CSV_Referencia.md](HELMVEC2_CSV_Referencia.md) | — | [HELMVEC2_Imagens_Referencia.md](HELMVEC2_Imagens_Referencia.md) |
| HELMVEC3 | [HELMVEC3_CSV_Referencia.md](HELMVEC3_CSV_Referencia.md) | — | [HELMVEC3_Imagens_Referencia.md](HELMVEC3_Imagens_Referencia.md) |
| FEM3D | [FEM3D_CSV_Referencia.md](FEM3D_CSV_Referencia.md) | — | [FEM3D_Imagens_Referencia.md](FEM3D_Imagens_Referencia.md) |

---

## Diagramas de Execução

| Arquivo | Conteúdo |
| --- | --- |
| [diagramas_execucao/README.md](diagramas_execucao/README.md) | Índice dos diagramas |
| [diagramas_execucao/Diagrama_Mestre_R1_a_R7.md](diagramas_execucao/Diagrama_Mestre_R1_a_R7.md) | Diagrama mestre de todos os solvers |
| [diagramas_execucao/Diagrama_R1_HELM10_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R1_HELM10_Sequencia_e_Arvore.md) | HELM10 |
| [diagramas_execucao/Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md) | HELMVEC |
| [diagramas_execucao/Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md) | HELMVEC1 |
| [diagramas_execucao/Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md) | HELMVEC2 |
| [diagramas_execucao/Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md) | HELMVEC3 |
| [diagramas_execucao/Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md) | FEM3D0 (denso) |
| [diagramas_execucao/Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md) | FEM3D1 (esparso) |

---

## READMEs por módulo de código (`src/`)

| Módulo | README | Formulação |
| --- | --- | --- |
| `src/helm10` | [src/helm10/README.md](../src/helm10/README.md) | Escalar 2D, Eq. (43) |
| `src/helmvec` | [src/helmvec/README.md](../src/helmvec/README.md) | Edge 2D transversal, Eq. (65) |
| `src/helmvec1` | [src/helmvec1/README.md](../src/helmvec1/README.md) | Misto 3 comp., Eq. (92) |
| `src/helmvec2` | [src/helmvec2/README.md](../src/helmvec2/README.md) | k0 dado β, Eq. (119) |
| `src/helmvec3` | [src/helmvec3/README.md](../src/helmvec3/README.md) | β dado k0, Eq. (136) |
| `src/fem3d0` | [src/fem3d0/README.md](../src/fem3d0/README.md) | 3D denso, Eq. (178) |
| `src/fem3d1` | [src/fem3d1/README.md](../src/fem3d1/README.md) | 3D esparso, Eq. (178) |
| `src/explicit` | [src/explicit/README.md](../src/explicit/README.md) | Formas fechadas locais |
| `src/article` | [src/article/README.md](../src/article/README.md) | API pública por equação global |

---

→ [README do repositório](../README.md)
