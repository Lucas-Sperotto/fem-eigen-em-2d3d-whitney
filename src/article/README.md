# `src/article` - APIs didaticas nomeadas pelas equacoes do TP-3485

Esta pasta concentra a camada mais diretamente ligada ao objetivo de
reconstruir, em C++, os programas do artigo NASA TP-3485.

## Objetivo

Dar ao projeto uma entrada publica em linguagem do proprio paper:

- `Eq. (43)` para `HELM10`
- `Eq. (65)` para `HELMVEC`
- `Eq. (92)` para `HELMVEC1`
- `Eq. (119)` para `HELMVEC2`
- `Eq. (136)` para `HELMVEC3`
- `Eq. (178)` para `FEM3D0/FEM3D1`

## Papel desta camada

- servir como fachada didatica sobre os montadores reais do repositorio;
- deixar explicito, no proprio codigo C++, onde cada sistema global do artigo
  comeca;
- reduzir a dependencia de nomes historicos de modulo quando a leitura correta
  pede o numero da equacao.

## Regra de arquitetura

- formas fechadas (`closed-form`) sao o fluxo principal;
- quadratura de Gauss fica como caminho auxiliar de verificacao;
- qualquer nova formulacao global do artigo deve ganhar uma API nomeada por
  equacao nesta pasta.

## Arquivo atual

- `tp3485_systems.hpp`
  - wrappers publicos para as equacoes globais principais do artigo.
