# Conclusoes, Apendice e Referencias

## Secao 4 - Conclusoes

As mensagens finais do artigo sao:

- o FEM foi formulado de modo consistente para varios problemas de autovalor em eletromagnetismo
- elementos de aresta resolvem o problema historico dos modos espurios nas formulacoes vetoriais
- a mesma filosofia vale em 2D e em 3D
- malhas triangulares e tetraedricas permitem tratar geometrias arbitrarias

## Relacao com a arquitetura do repositorio

O fechamento do artigo combina perfeitamente com a estrutura atual:

- `HELM10`: caminho escalar nodal
- `HELMVEC`, `HELMVEC1`, `HELMVEC2`, `HELMVEC3`: sequencia vetorial 2D
- `FEM3D0`, `FEM3D1`: sequencia vetorial 3D

## Apendice - Programas do artigo

- `HELM10`
  - solver 2D escalar para guias homogeneos
- `HELMVEC`
  - solver 2D vetorial transversal
- `HELMVEC1`
  - solver 2D misto para guias inhomogeneos no cutoff
- `HELMVEC2`
  - solver 2D para `k0` com `beta` dado
- `HELMVEC3`
  - solver 2D para `beta` com `k0` dado
- `FEM3D0`
  - solver 3D denso
- `FEM3D1`
  - solver 3D explorando simetria/esparsidade

## Referencias mais importantes para a leitura do projeto

- ref. 6:
  - matriz universal triangular escalar
- ref. 9:
  - valores analiticos usados nas tabelas de guias de onda
- ref. 10:
  - referencia classica para linha coaxial
- ref. 12:
  - base de Whitney / elementos de aresta
- ref. 13:
  - formulacao full-wave para guias dieletricos com elementos tangenciais
- ref. 15:
  - motivacao para problemas vetoriais 2D sem modos espurios
- ref. 16:
  - modelagem 3D de cavidades inhomogeneas com edge elements

## Fechamento teorico para o projeto

Depois desta trilha, a leitura do codigo pode seguir um criterio simples:

1. entender qual problema de autovalor cada modulo resolve
2. localizar a equacao global principal da secao correspondente
3. localizar os blocos locais fechados que alimentam essa equacao
4. confrontar a teoria do artigo com os detalhes de BC, orientacao de aresta e filtragem numerica do repositorio
