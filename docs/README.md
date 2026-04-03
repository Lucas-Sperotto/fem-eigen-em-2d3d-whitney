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
- [Casos_de_Teste_do_Artigo.md](Casos_de_Teste_do_Artigo.md): tabela consolidada com os casos numéricos do artigo.
- [Matriz_Casos_Executaveis_Arvore_de_Chamada.md](Matriz_Casos_Executaveis_Arvore_de_Chamada.md): ponte entre casos do artigo, executáveis atuais e raízes da árvore de chamada.
- [diagramas_execucao/README.md](diagramas_execucao/README.md): caderno de navegação da coleção de diagramas de execução.
- [diagramas_execucao/Diagrama_Mestre_R1_a_R7.md](diagramas_execucao/Diagrama_Mestre_R1_a_R7.md): visão global resumida da coleção `R1 -> R7`.
- [Diagrama_R1_HELM10_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R1_HELM10_Sequencia_e_Arvore.md) a [Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md](diagramas_execucao/Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md): diagramas-base por família de executáveis.

## Relação com o código

- [2.1 Guias de Onda Homogêneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md) -> `src/helm10`
- [2.2.1 Campos Vetoriais Transversais](traducao/2.2.1_Guia%20de%20Onda%20Campos%20Vetoriais.md) -> `src/helmvec`
- [2.2.2 Três Componentes](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md) -> `src/helmvec1`
- [2.2.3 Determinação do Número de Onda](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md) -> `src/helmvec2`
- [2.2.4 Características de Dispersão](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md) -> `src/helmvec3`
- [3. Problemas Tridimensionais](traducao/3_Problemas_Tridimensionais.md) -> `src/fem3d0` e `src/fem3d1`

Se a ideia for estudar a teoria antes de mergulhar na implementação, esta pasta deve vir antes dos READMEs em `src/`.
