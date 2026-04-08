# Documentação Teórica

Esta pasta reúne a trilha principal de leitura do artigo NASA TP-3485 dentro do projeto. Os arquivos numerados em `docs/traducao/` preservam a tradução-base revisada manualmente e receberam notas de leitura, passagens intermediárias entre equações e observações editoriais pontuais.

## Fontes usadas na revisão

- PDF local: [19950011772.pdf](refs/19950011772.pdf)
- Cópia oficial no NTRS da NASA: `https://ntrs.nasa.gov/citations/19950011772`
- Documentação técnica dos módulos em `src/*/README.md`, usada apenas para checagem de coerência matemática com a implementação

## Critérios desta revisão

- manter a numeração original das equações do artigo
- separar a tradução por bloco lógico, e não por fragmentos misturados
- apontar quando a documentação antiga misturava seções, pulava numeração ou deformava fórmulas
- deixar claro o que é teoria do artigo e o que já é interpretação voltada ao código

## Trilha principal

1. [0.0 Capa](traducao/0.0_Capa.md)
2. [0.1 Sumário](traducao/0.1_Sumario.md)
3. [0.2 Tabelas](traducao/0.2_Tabelas.md)
4. [0.3 Figuras](traducao/0.3_Figuras.md)
5. [0.4 Lista de Símbolos](traducao/0.4_Simbolos.md)
6. [0.5 Abstract e Escopo do Artigo](traducao/0.5_Abstract.md)
7. [1. Introdução](traducao/1_Introducao.md)
8. [2. Problemas Bidimensionais](traducao/2_Problemas_Bidimensionais.md)
9. [2.1 Guias de Onda Homogêneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md)
10. [2.2 Guias de Onda Não Homogêneos](traducao/2.2_Guia_Nao_Homogeneo.md)
11. [2.2.1 Campos Vetoriais Transversais](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md)
12. [2.2.2 Três Componentes](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
13. [2.2.3 Determinação do Número de Onda](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md)
14. [2.2.4 Características de Dispersão](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
15. [3. Problemas Tridimensionais](traducao/3_Problemas_Tridimensionais.md)
16. [4. Considerações Finais, Apêndice e Referências](traducao/4_Consideracoes_Finais_Apendice_Referencias.md)
17. [5. Página de Documentação do Relatório](traducao/5_Report.md)

Todos os arquivos da trilha principal agora terminam com links de `Anterior` e `Próximo`, para permitir leitura linear sem voltar ao índice a cada etapa.

## Mapa rápido

| Parte | Arquivo | Função no conjunto |
| --- | --- | --- |
| Material preliminar | [0.0 Capa](traducao/0.0_Capa.md) a [0.5 Abstract](traducao/0.5_Abstract.md) | Identificação do relatório, sumário, figuras, tabelas, símbolos e escopo do artigo |
| Abertura conceitual | [1. Introdução](traducao/1_Introducao.md) | Contexto histórico, motivação e problema dos modos espúrios |
| Visão 2D | [2. Problemas Bidimensionais](traducao/2_Problemas_Bidimensionais.md) | Ponte para as formulações 2D escalar e vetorial |
| Formulação escalar | [2.1 Guias de Onda Homogêneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md) | Forma fraca, discretização nodal e reconstrução de campos |
| Visão vetorial | [2.2 Guias de Onda Não Homogêneos](traducao/2.2_Guia_Nao_Homogeneo.md) | Introdução aos elementos de aresta e ao problema vetorial |
| Formulação vetorial transversal | [2.2.1 Campos Vetoriais Transversais](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md) | Base algébrica com elementos de aresta em 2D |
| Formulação mista no cutoff | [2.2.2 Três Componentes](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md) | Acoplamento entre arestas e nós sem termos fora da diagonal |
| `k_0` com `\beta` fixado | [2.2.3 Determinação do Número de Onda](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md) | Problema generalizado para obter `k_0` |
| `\beta` com `k_0` fixado | [2.2.4 Características de Dispersão](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md) | Problema de dispersão e curvas de propagação |
| Formulação 3D | [3. Problemas Tridimensionais](traducao/3_Problemas_Tridimensionais.md) | Cavidades, tetraedros e bases de aresta tridimensionais |
| Fechamento | [4. Considerações Finais, Apêndice e Referências](traducao/4_Consideracoes_Finais_Apendice_Referencias.md) | Síntese, programas do apêndice e bibliografia |
| Ficha do relatório | [5. Página de Documentação do Relatório](traducao/5_Report.md) | Metadados institucionais do NASA TP-3485 |

## Mapa compacto do artigo

| Parte | Tema | Equação global principal | Programa do artigo |
| --- | --- | --- | --- |
| 2.1 | Guia homogêneo, formulação escalar | Eq. (34) / Eq. (43) no repositório | `HELM10` |
| 2.2.1 | Guia homogêneo, formulação vetorial transversal | Eq. (65) | `HELMVEC` |
| 2.2.2 | Guia inhomogêneo, três componentes no cutoff | Eq. (92) | `HELMVEC1` |
| 2.2.3 | `k_0` para `\beta` dado | Eq. (119) | `HELMVEC2` |
| 2.2.4 | `\beta` para `k_0` dado | Eq. (136) | `HELMVEC3` |
| 3.1 | Cavidades 3D com elementos de aresta | Eq. (178) | `FEM3D0` / `FEM3D1` |

## Materiais de apoio

- [19950011772.pdf](refs/19950011772.pdf): PDF original do artigo usado como referência de confronto.
- [Rastreabilidade_Equacoes_Artigo_Codigo.md](Rastreabilidade_Equacoes_Artigo_Codigo.md): mapa central `equação -> função -> arquivo` para auditoria didática da implementação em C++.
- [Casos_de_Teste_do_Artigo.md](Casos_de_Teste_do_Artigo.md): tabela consolidada com os casos numéricos do artigo.
- [Matriz_Casos_Executaveis_Arvore_de_Chamada.md](Matriz_Casos_Executaveis_Arvore_de_Chamada.md): ponte entre casos do artigo, executáveis atuais e raízes da árvore de chamada.
- [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md): tabela direta `executável -> entrada -> saída`, útil para operação e estudo.
- [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md): referencia central dos arquivos `linop/*.csv`, cobrindo matrizes em CRS, autovalores ordenados e autovetores ordenados das familias `HELM10` a `FEM3D1`.
- [FEM3D_CSV_Referencia.md](FEM3D_CSV_Referencia.md): referencia didatica das saidas 3D de `FEM3D0` e `FEM3D1`, cobrindo `run.log`, `run_timing.csv`, `modes.csv`, `fields.csv`, `vtk/` e `linop/` por caso.
- [FEM3D_Imagens_Referencia.md](FEM3D_Imagens_Referencia.md): guia do script `scripts/fem3d.py`, que gera resumos modais e projecoes ortogonais `XY/XZ/YZ` dos campos 3D exportados por modo casado.
- [HELM10_CSV_Modos_Referencia.md](HELM10_CSV_Modos_Referencia.md): glossario didatico das colunas do CSV de modos do `HELM10`, ligado à Secao 2.1, incluindo a explicacao de `rho_abs` e do matching por correlacao de massa.
- [HELM10_CSV_Campos_Referencia.md](HELM10_CSV_Campos_Referencia.md): glossario didatico das colunas do CSV de campos nodais do `HELM10`, ligado à Secao 2.1.
- [HELM10_Imagens_Referencia.md](HELM10_Imagens_Referencia.md): guia do script `scripts/helm10.py`, que lê os CSVs do `HELM10` e gera isolinhas, diagramas de setas e resumo de erro por modo.
- [HELMVEC_CSV_Modos_Referencia.md](HELMVEC_CSV_Modos_Referencia.md): glossario didatico das colunas do CSV de modos do `HELMVEC`, ligado à Secao 2.2.1, incluindo a explicacao de `rho_abs` e do matching por correlacao de massa vetorial.
- [HELMVEC_CSV_Campos_Referencia.md](HELMVEC_CSV_Campos_Referencia.md): glossario didatico das colunas do CSV de campos por celula do `HELMVEC`, ligado à Secao 2.2.1.
- [HELMVEC_Imagens_Referencia.md](HELMVEC_Imagens_Referencia.md): guia do script `scripts/helmvec.py`, que lê os CSVs do `HELMVEC` e gera mapas de magnitude, diagramas de setas e resumos de erro e `rho`.
- [HELMVEC1_CSV_Modos_Referencia.md](HELMVEC1_CSV_Modos_Referencia.md): glossario didatico das colunas do CSV de modos do `HELMVEC1`, ligado à Secao 2.2.2, incluindo a explicacao da classificacao por energia de bloco e do `rho_abs` por correlacao de massa no bloco dominante.
- [HELMVEC1_CSV_Campos_Referencia.md](HELMVEC1_CSV_Campos_Referencia.md): glossario didatico dos CSVs espaciais do `HELMVEC1`, separando os casos `edge_vector_cell` e `scalar_nodal`.
- [HELMVEC1_Imagens_Referencia.md](HELMVEC1_Imagens_Referencia.md): guia do script `scripts/helmvec1.py`, que lê o `modes.csv` do `HELMVEC1`, os CSVs espaciais e os VTKs por modo para gerar imagens-resumo e imagens espaciais coerentes com o bloco dominante.
- [HELMVEC2_CSV_Referencia.md](HELMVEC2_CSV_Referencia.md): referencia didatica dos CSVs do `HELMVEC2`, ligando o espectro de candidatos da Eq. `(119)` ao matching final da Tabela 8 e aos artefatos espaciais `Et`/`Ez` de cada modo casado.
- [HELMVEC2_Imagens_Referencia.md](HELMVEC2_Imagens_Referencia.md): guia do script `scripts/helmvec2.py`, que lê os CSVs do `HELMVEC2` e gera tanto o resumo grafico da Tabela 8 quanto as imagens espaciais `Et`/`Ez` por modo casado.
- [HELMVEC3_CSV_Referencia.md](HELMVEC3_CSV_Referencia.md): referencia didatica dos CSVs do `HELMVEC3`, ligando a Eq. `(136)` aos dados estruturados da Figura 12 / Tabela 9, da Figura 13 / Tabela 10 e aos artefatos espaciais `Et`/`Ez` por ponto exportado.
- [HELMVEC3_Imagens_Referencia.md](HELMVEC3_Imagens_Referencia.md): guia do script `scripts/helmvec3.py`, que lê os CSVs do `HELMVEC3` e gera tanto os resumos graficos quanto as imagens espaciais `Et`/`Ez` por ponto exportado.
- [diagramas_execucao/README.md](diagramas_execucao/README.md): caderno de navegação da coleção de diagramas de execução.
- [diagramas_execucao/Diagrama_Mestre_R1_a_R7.md](diagramas_execucao/Diagrama_Mestre_R1_a_R7.md): visão global resumida da coleção `R1 -> R7`.
- [Diagrama_R1_HELM10_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R1_HELM10_Sequencia_e_Arvore.md) a [Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md): diagramas-base por família de executáveis.

## Relação com o código

- [2.1 Guias de Onda Homogêneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md) -> `src/helm10`
  - apoio operacional: [HELM10_CSV_Modos_Referencia.md](HELM10_CSV_Modos_Referencia.md)
  - apoio operacional nodal: [HELM10_CSV_Campos_Referencia.md](HELM10_CSV_Campos_Referencia.md)
  - apoio grafico: [HELM10_Imagens_Referencia.md](HELM10_Imagens_Referencia.md)
- [2.2.1 Campos Vetoriais Transversais](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md) -> `src/helmvec`
  - apoio operacional: [HELMVEC_CSV_Modos_Referencia.md](HELMVEC_CSV_Modos_Referencia.md)
  - apoio operacional por celula: [HELMVEC_CSV_Campos_Referencia.md](HELMVEC_CSV_Campos_Referencia.md)
  - apoio grafico: [HELMVEC_Imagens_Referencia.md](HELMVEC_Imagens_Referencia.md)
- [2.2.2 Três Componentes](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md) -> `src/helmvec1`
  - apoio operacional: [HELMVEC1_CSV_Modos_Referencia.md](HELMVEC1_CSV_Modos_Referencia.md)
  - apoio operacional espacial: [HELMVEC1_CSV_Campos_Referencia.md](HELMVEC1_CSV_Campos_Referencia.md)
  - apoio grafico: [HELMVEC1_Imagens_Referencia.md](HELMVEC1_Imagens_Referencia.md)
  - apoio operacional: [Rastreabilidade_Equacoes_Artigo_Codigo.md](Rastreabilidade_Equacoes_Artigo_Codigo.md)
  - apoio de navegacao: [Matriz_Casos_Executaveis_Arvore_de_Chamada.md](Matriz_Casos_Executaveis_Arvore_de_Chamada.md)
  - apoio de uso: [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
- [2.2.3 Determinação do Número de Onda](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md) -> `src/helmvec2`
  - apoio operacional: [HELMVEC2_CSV_Referencia.md](HELMVEC2_CSV_Referencia.md)
  - apoio grafico: [HELMVEC2_Imagens_Referencia.md](HELMVEC2_Imagens_Referencia.md)
  - apoio espacial: campos `Et`/`Ez` exportados por modo casado
  - apoio de uso: [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
- [2.2.4 Características de Dispersão](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md) -> `src/helmvec3`
  - apoio operacional: [HELMVEC3_CSV_Referencia.md](HELMVEC3_CSV_Referencia.md)
  - apoio grafico: [HELMVEC3_Imagens_Referencia.md](HELMVEC3_Imagens_Referencia.md)
  - apoio espacial: campos `Et`/`Ez` exportados por ponto da Figura 12, do preview e da Tabela 10
  - apoio espectral: [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md)
  - apoio de uso: [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
- [3. Problemas Tridimensionais](traducao/3_Problemas_Tridimensionais.md) -> `src/fem3d0` e `src/fem3d1`
  - apoio operacional: [FEM3D_CSV_Referencia.md](FEM3D_CSV_Referencia.md)
  - apoio grafico: [FEM3D_Imagens_Referencia.md](FEM3D_Imagens_Referencia.md)
  - apoio espectral: [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md)

Se a ideia for estudar a teoria antes de mergulhar na implementação, esta pasta deve vir antes dos READMEs em `src/`.
