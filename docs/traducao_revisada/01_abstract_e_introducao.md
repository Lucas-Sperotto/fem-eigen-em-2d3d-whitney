# Abstract e Introducao

## Abstract revisado

O artigo apresenta uma trilha completa de formulacoes FEM para problemas de autovalor em eletromagnetismo, com foco em dois grupos:

- problemas 2D em guias de onda
- problemas 3D em cavidades

O objetivo declarado dos autores e construir ferramentas FEM reutilizaveis para eletromagnetismo, com um primeiro passo em problemas de contorno fechado antes de atacar problemas 3D de espalhamento/radiacao.

## Introducao revisada

O argumento central da introducao e:

1. FEM ja era maduro em outras areas de engenharia.
2. Em eletromagnetismo, o uso de formulacoes vetoriais nodais sofria com modos espurios.
3. O uso de elementos de aresta reabre o caminho para formulacoes vetoriais robustas.
4. O artigo organiza uma progressao didatica:
   - primeiro o caso escalar 2D
   - depois os casos vetoriais 2D, do mais simples ao acoplado
   - por fim a extensao 3D em tetraedros

## Leitura guiada da estrutura do artigo

- Secao 2.1:
  - resolve guias homogeneos 2D por potencial escalar
  - separa TE e TM por condicoes de contorno
- Secao 2.2:
  - abandona a formulacao puramente escalar
  - introduz elementos de aresta 2D
  - constroi uma sequencia incremental de quatro problemas
- Secao 3.1:
  - estende a ideia de elementos de aresta para tetraedros 3D
  - formula o problema de autovalor para cavidades

## O que esta secoes motivam no codigo

- `src/helm10`: implementa a trilha escalar da Secao 2.1
- `src/helmvec`: cobre a Secao 2.2.1
- `src/helmvec1`: cobre a Secao 2.2.2
- `src/helmvec2`: cobre a Secao 2.2.3
- `src/helmvec3`: cobre a Secao 2.2.4
- `src/fem3d0` e `src/fem3d1`: cobrem a Secao 3.1

## Diagnostico da traducao antiga

Na documentacao legada, a introducao estava relativamente preservada, mas faltavam tres coisas importantes:

- um mapa explicito de como as secoes conversam entre si
- um mapa claro das equacoes globais principais
- uma separacao entre traducao do artigo e interpretacao voltada ao repositorio

Esta trilha revisada corrige essa falta de estrutura logo no inicio para facilitar a leitura teorica antes da leitura do codigo.
