# Casos de Teste do Artigo

Este arquivo resume apenas os casos de teste explicitamente apresentados pelo autor nos exemplos numéricos do artigo NASA TP-3485. A ideia aqui é separar com clareza:

- os casos efetivamente publicados no artigo;
- o programa histórico associado a cada caso;
- a tabela ou figura onde o caso aparece;
- o executável atual do repositório que reproduz ou aproxima esse caso.

Não entram aqui casos extras criados pelo repositório para validação ampliada, comparação entre backends, varreduras de malha ou figuras auxiliares que não fazem parte do conjunto original de testes do autor.

## Resumo rápido

- Total de casos numéricos apresentados pelo autor: `14`
- Casos bidimensionais: `10`
- Casos tridimensionais: `4`

## Tabela geral

| ID | Seção do artigo | Programa do artigo | Caso de teste | Geometria / configuração | Saída comparada pelo autor | Figura / tabela | Malha citada no texto | Executável atual do repositório |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 2.1 | `HELM10` | Guia de onda retangular homogêneo | Retangular com `a_r/b_r = 2` | Números de onda de corte `k_c` e campos modais | Figura 4, Tabela 1, Figura 5 | `400` triângulos | `./build/helm10_rect 14 14 8` |
| 2 | 2.1 | `HELM10` | Guia de onda circular homogêneo | Circular com raio unitário | Números de onda de corte `k_c` e campos modais | Figura 6, Tabela 2, Figura 7 | `200` triângulos | `./build/helm10_circle 10 48 8` |
| 3 | 2.1 | `HELM10` | Linha coaxial homogênea | Coaxial com `r_2/r_1 = 4` | Números de onda de corte `k_c` e campos TE/TM superiores | Figura 8, Tabela 3, Figura 9 | `340` triângulos | `./build/helm10_coax 10 48 8` |
| 4 | 2.2.1 | `HELMVEC` | Guia de onda retangular vetorial | Retangular com ar e `a_r/b_r = 2` | Números de onda de corte `k_c` para modos TE/TM | Figura 4, Tabela 4 | `400` triângulos | `./build/edge_rect 14 14 8` |
| 5 | 2.2.1 | `HELMVEC` | Guia de onda circular vetorial | Circular com ar e raio unitário | Números de onda de corte `k_c` para modos TE/TM | Tabela 5 | `200` triângulos | `./build/edge_circle 10 48 8` |
| 6 | 2.2.2 | `HELMVEC1` | Guia de onda retangular misto | Retangular com ar e `a_r/b_r = 2` | Números de onda de corte `k_c` nas formulações `E` e `H` | Figura 4, Tabela 6 | Não explicitada nesse trecho | `./build/mixed_rect 12 6` |
| 7 | 2.2.2 | `HELMVEC1` | Guia de onda circular misto | Circular com ar e raio unitário | Números de onda de corte `k_c` nas formulações `E` e `H` | Tabela 7 | Não explicitada nesse trecho | `./build/mixed_circle 10 48` |
| 8 | 2.2.3 | `HELMVEC2` | Guia de onda quadrado parcialmente preenchido | Guia quadrado parcialmente preenchido com `beta = 10` | Primeiros `10` números de onda | Figura 11, Tabela 8 | Malha atual documentada no repositório: `6x6` (`72` triângulos) | `./build/helmvec2_rect 10 6 6` |
| 9 | 2.2.4 | `HELMVEC3` | Dispersão em guia retangular parcialmente preenchido, exemplo 1 | `b_r/a_r = 0.45`, `d/a_r = 0.5`, `\varepsilon_r = 2.45` | Razão `beta/k_0` | Figura 12, Tabela 9 | Não explicitada no artigo traduzido | `./build/helmvec3_rect 0.20 10 5` |
| 10 | 2.2.4 | `HELMVEC3` | Dispersão em guia retangular parcialmente preenchido, exemplo 2 | `b_r/a_r = 0.45`, `\varepsilon_r = 2.45`, vários valores de `d/a_r` | Ramos de dispersão `beta/k_0` | Figura 13, Tabela 10 | `100` triângulos | `./build/helmvec3_rect 0.20 10 5` |
| 11 | 3.1 | `FEM3D1` | Cavidade retangular preenchida com ar | Cavidade retangular homogênea | Autovalores de ressonância `k_0` | Figura 15, Tabela 12 | `343` tetraedros | `./build/fem3d1_rect --air` ou `./build/fem3d0_rect --air` |
| 12 | 3.1 | `FEM3D1` | Cavidade retangular semi-preenchida | Preenchimento dielétrico em `z = 0.5` a `1` cm, `\varepsilon_r = 2.0` | Autovalores de ressonância `k_0` | Figura 16, Tabela 13 | `615` tetraedros | `./build/fem3d1_rect --half` ou `./build/fem3d0_rect --half` |
| 13 | 3.1 | `FEM3D1` | Cavidade cilíndrica circular com ar | Cavidade cilíndrica homogênea | Autovalores de ressonância `k_0` | Figura 17, Tabela 14 | `633` tetraedros | `./build/fem3d1_rect --cyl` ou `./build/fem3d0_rect --cyl` |
| 14 | 3.1 | `FEM3D1` | Cavidade esférica | Esfera de raio `1 cm` | Autovalores de ressonância `k_0` | Tabela 15 | `473` tetraedros | `./build/fem3d1_rect --sphere` ou `./build/fem3d0_rect --sphere` |

## Observações importantes

- O artigo reaproveita algumas mesmas geometrias em formulações diferentes. Por isso, um mesmo tipo de guia aparece mais de uma vez na tabela, mas com problema espectral distinto.
- A Seção `2.2.1` reutiliza a geometria retangular e circular da parte escalar para validar a formulação com elementos de aresta.
- A Seção `2.2.2` continua usando os mesmos tipos básicos de guia, agora para validar a formulação mista com três componentes.
- A Seção `2.2.3` introduz o primeiro caso realmente acoplado com material não homogêneo e `beta` especificado.
- A Seção `2.2.4` usa o mesmo tipo geral de estrutura parcialmente preenchida, mas passa a buscar `beta` para frequência dada.
- Na Seção `3.1`, o artigo apresenta quatro casos de cavidades tridimensionais, mas o texto preservado no repositório registra os resultados numéricos pelo `FEM3D1`; o repositório atual também oferece o `FEM3D0` como baseline denso.

## Casos do repositório que não pertencem ao conjunto original do autor

Os itens abaixo existem no repositório por conveniência de validação ou ampliação da cobertura, mas não entram na lista original de casos publicados:

- `edge_coax`
- `mixed_coax`
- snapshots auxiliares de validação 2D
- campanhas de benchmark entre `gauss` e `closed-form`
- varreduras de malha e campanhas completas em `docs/results/`

## Fonte principal usada nesta consolidação

- [2_Problemas_Bidimensionais.md](traducao/2_Problemas_Bidimensionais.md)
- [2.2_Guia_Nao_Homogeneo.md](traducao/2.2_Guia_Nao_Homogeneo.md)
- [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
- [2.2.3_Determinação do número de onda.md](traducao/2.2.3_Determinação%20do%20número%20de%20onda.md)
- [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
- [README.md](/home/sperotto/tp3485-fem-eigen-em/README.md)
