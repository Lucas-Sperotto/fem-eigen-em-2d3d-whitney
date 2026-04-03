# Colecao de Diagramas de Execucao

Esta pasta funciona como caderno de navegacao da colecao de diagramas-base do projeto. Ela organiza a leitura deles como uma trilha unica, do problema mais simples ao mais complexo.

## 1) Como ler esta colecao

Os diagramas foram construidos com uma logica progressiva:

1. primeiro, entender a espinha minima de montagem, solve e matching em 2D;
2. depois, observar como o espaco discreto muda de nodos para arestas;
3. em seguida, acompanhar o aparecimento dos sistemas mistos e acoplados;
4. por fim, transportar essa mesma filosofia para o dominio 3D, primeiro em denso e depois em esparso.

Se a ideia for estudar o repositorio como um curso de implementacao FEM, esta e a ordem de leitura recomendada.

## 2) Ordem recomendada dos diagramas

| Raiz | Arquivo principal | Casos cobertos | O que aprender aqui |
| --- | --- | --- | --- |
| `R1` | [Diagrama_R1_HELM10_Sequencia_e_Arvore.md](Diagrama_R1_HELM10_Sequencia_e_Arvore.md) | Casos 1 a 3 | forma escalar 2D, BCs TE/TM, matching analitico e reconstrução de campo |
| `R2` | [Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md](Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md) | Casos 4 e 5 | bases de Whitney 2D, orientacao de arestas e campo por celula |
| `R3` | [Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md](Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md) | Casos 6 e 7 | sistema misto no cutoff e separacao modal por energia de bloco |
| `R4` | [Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md](Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md) | Caso 8 | problema acoplado para `k0` com `beta` dado, filtro fisico e matching da Tabela 8 |
| `R5` | [Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md](Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md) | Casos 9 e 10 | dispersao, solve repetido por ponto amostrado e rastreamento de ramos |
| `R6` | [Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md](Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md) | Casos 11 a 14 | formulacao vetorial 3D densa, tetraedros e matching com degenerescencia |
| `R7` | [Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md](Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md) | Casos 11 a 14 | mesma fisica 3D com armazenamento esparso simetrico |

## 3) Mapa resumido da progressao

```text
R1: escalar 2D nodal
 -> R2: vetorial 2D com arestas
 -> R3: misto em blocos no cutoff
 -> R4: acoplado para k0 com beta dado
 -> R5: acoplado para beta com k0 dado
 -> R6: vetorial 3D denso
 -> R7: vetorial 3D esparso
```

## 4) Guia rapido por tipo de pergunta

Se a pergunta for "onde nasce o FEM mais basico do projeto?":

- comece por [Diagrama_R1_HELM10_Sequencia_e_Arvore.md](Diagrama_R1_HELM10_Sequencia_e_Arvore.md)

Se a pergunta for "onde entram os elementos de aresta?":

- leia [Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md](Diagrama_R2_HELMVEC_Sequencia_e_Arvore.md)

Se a pergunta for "onde o codigo passa a acoplar `Et` e `Ez`?":

- leia [Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md](Diagrama_R3_HELMVEC1_Sequencia_e_Arvore.md), [Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md](Diagrama_R4_HELMVEC2_Sequencia_e_Arvore.md) e [Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md](Diagrama_R5_HELMVEC3_Sequencia_e_Arvore.md)

Se a pergunta for "onde o projeto realmente entra em 3D?":

- comece por [Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md](Diagrama_R6_FEM3D0_Sequencia_e_Arvore.md)

Se a pergunta for "onde aparece a estrategia esparsa?":

- leia [Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md](Diagrama_R7_FEM3D1_Sequencia_e_Arvore.md)

## 5) Arquivos de apoio desta colecao

- [Diagrama_Mestre_R1_a_R7.md](Diagrama_Mestre_R1_a_R7.md): visao global resumida da colecao inteira.
- [Matriz_Casos_Executaveis_Arvore_de_Chamada.md](../Matriz_Casos_Executaveis_Arvore_de_Chamada.md): ponte entre artigo, executaveis e raizes reais de execucao.
- [Casos_de_Teste_do_Artigo.md](../Casos_de_Teste_do_Artigo.md): lista consolidada dos casos numericos do autor.

## 6) Papel desta pasta no projeto

Este caderno existe para reduzir a friccao de leitura. Esta pasta concentra a navegacao orientada para implementacao e abriga a colecao completa de diagramas da trilha `R1 -> R7`.
