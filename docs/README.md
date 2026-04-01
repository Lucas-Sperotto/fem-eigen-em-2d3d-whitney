# Documentacao Teorica

Esta pasta concentra a trilha de leitura do artigo NASA TP-3485 dentro do projeto. Os arquivos `.md` na raiz de `docs/` preservam a traducao-base revisada manualmente e agora incluem comentarios complementares, passagens intermediarias entre equacoes, observacoes de consistencia e notas de leitura para estudo.

O principio adotado nesta revisao foi simples: preservar o texto traduzido e apenas acrescentar explicacoes, interpretacoes e alertas onde o desenvolvimento pedia mais clareza.

## Ordem recomendada de leitura

1. [0.0_Lista_de_Simbolos.md](0.0_Lista_de_Simbolos.md)
2. [0_Abstract.md](0_Abstract.md)
3. [1_Introducao.md](1_Introducao.md)
4. [2_Problemas_Bidimensionais.md](2_Problemas_Bidimensionais.md)
5. [2.2.md](2.2.md)
6. [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
7. [2.2.3_Determinação do número de onda.md](2.2.3_Determinação%20do%20número%20de%20onda.md)
8. [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
9. [3_Problemas_Tridimensionais.md](3_Problemas_Tridimensionais.md)
10. [4_Consideracoes_Finais_Apendice_Referencias.md](4_Consideracoes_Finais_Apendice_Referencias.md)

## Mapa rapido do artigo

| Parte do estudo | Arquivo | Papel no conjunto |
| --- | --- | --- |
| Lista de simbolos | [0.0_Lista_de_Simbolos.md](0.0_Lista_de_Simbolos.md) | Convenioes, notacao e blocos matriciais usados ao longo da teoria |
| Abstract | [0_Abstract.md](0_Abstract.md) | Escopo do artigo e promessa metodologica |
| Introducao | [1_Introducao.md](1_Introducao.md) | Contexto historico, modos espurios e motivacao para elementos de aresta |
| Secao 2.1 | [2_Problemas_Bidimensionais.md](2_Problemas_Bidimensionais.md) | Formulacao escalar 2D, forma fraca, discretizacao nodal e campos TE/TM |
| Secao 2.2.1 | [2.2.md](2.2.md) | Formulacao vetorial transversal 2D com elementos de aresta |
| Secao 2.2.2 | [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md) | Sistema misto vetorial + escalar no cutoff |
| Secao 2.2.3 | [2.2.3_Determinação do número de onda.md](2.2.3_Determinação%20do%20número%20de%20onda.md) | Problema acoplado para determinar `k0` com `beta` especificado |
| Secao 2.2.4 | [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md) | Problema de dispersao para determinar `beta` com `k0` especificado |
| Secao 3.1 | [3_Problemas_Tridimensionais.md](3_Problemas_Tridimensionais.md) | Cavidades 3D, tetraedros e bases de aresta tridimensionais |
| Encerramento | [4_Consideracoes_Finais_Apendice_Referencias.md](4_Consideracoes_Finais_Apendice_Referencias.md) | Consideracoes finais, apendice dos programas e referencias |

## Materiais de apoio

- [19950011772.pdf](19950011772.pdf): PDF original do artigo usado como referencia de confronto.
- [Casos_de_Teste_do_Artigo.md](Casos_de_Teste_do_Artigo.md): tabela consolidada com todos os casos numéricos apresentados pelo autor.
- [results/README.md](results/README.md): resultados curados, validacoes e figuras do projeto.
- Arquivos `figura*.png` nesta pasta: imagens extraidas/associadas a partes da documentacao original.

## Estado atual da revisao

- A traducao-base foi preservada nos arquivos principais.
- Foram adicionadas explicacoes de leitura, comentarios didaticos e passagens intermediarias entre equacoes.
- Irregularidades provaveis de OCR, notacao e placeholders de figuras foram apontadas dentro dos arquivos, sem apagar o texto original.
- Os nomes de alguns arquivos foram ajustados para refletir melhor seu papel no conjunto.

## Relacao com o codigo

Esta trilha teorica conversa diretamente com os modulos do projeto:

- `2_Problemas_Bidimensionais.md` -> `src/helm10`
- `2.2.md` -> `src/helmvec`
- `2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md` -> `src/helmvec1`
- `2.2.3_Determinação do número de onda.md` -> `src/helmvec2`
- `2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md` -> `src/helmvec3`
- `3_Problemas_Tridimensionais.md` -> `src/fem3d0` e `src/fem3d1`

Se a ideia for estudar a teoria antes de mergulhar na implementacao, esta pasta deve ser lida primeiro e o `README.md` da raiz deve entrar em seguida como ponte para os executaveis, scripts e validacoes.
