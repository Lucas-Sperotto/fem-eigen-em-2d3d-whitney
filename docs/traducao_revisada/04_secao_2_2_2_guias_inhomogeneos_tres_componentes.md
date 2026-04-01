# Secao 2.2.2 - Guias Inhomogeneos com Tres Componentes

## Escopo

Esta secao introduz a formulacao mista para guias inhomogeneos no cutoff. Ela combina:

- campo transversal em elementos de aresta
- componente longitudinal em funcoes nodais escalares

No repositorio, ela corresponde ao programa `HELMVEC1`.

## Ideia central

O artigo desacopla o problema no cutoff em dois blocos compatíveis:

- bloco vetorial transversal
- bloco escalar longitudinal

Isso produz um sistema bloco-diagonal que preserva a estrutura fisica do problema e evita forcar uma discretizacao nodal em todo o campo vetorial.

## Revisao das equacoes

- Eq. (78): equacao vetorial de onda usada como ponto de partida.
- Eq. (79): decomposicao `E = E_t + z_hat E_z`.
- Eq. (80): parte transversal do problema.
- Eq. (81): parte longitudinal escalar.
- Eq. (82): forma fraca da parte transversal.
- Eq. (83): forma fraca da parte longitudinal.
- Eq. (84): definicao da funcao teste total `T = T_t + z_hat T_z`.
- Eq. (85): expansao de `E_t` nas bases de aresta.
- Eq. (86): escolha de Galerkin para `T_t`.
- Eq. (87): expansao de `E_z` nas funcoes nodais.
- Eq. (88): escolha de Galerkin para `T_z`.
- Eq. (89): equacao elemental do bloco transversal.
- Eq. (90): equacao elemental do bloco longitudinal.
- Eq. (91): sistema elemental em blocos.
- Eq. (92): sistema global em blocos `diag(St,Sz)` e `diag(Tt,Tz)`.

## Resultado teorico principal

O sistema final da secao e:

```text
[ St   0 ] [e_t] = kc^2 [ Tt   0 ] [e_t]
[  0  Sz ] [e_z]        [  0  Tz ] [e_z]
```

Ele nao introduz acoplamento cruzado entre `e_t` e `e_z` porque o problema esta no cutoff.

## Reuso das secoes anteriores

- Eq. (88) da secao 2.2.2 reusa a forma vetorial de 2.2.1.
- Eq. (89) reusa a forma escalar de 2.1.
- Eq. (92) e, na pratica, uma justaposicao dos dois problemas anteriores.

## Figuras e tabelas desta secao

- Tabela 6: guia retangular
- Tabela 7: guia circular

## Leitura para o codigo

- Esta secao explica por que `src/helmvec1` pode ser montado em termos de dois blocos independentes.
- O bloco `St/Tt` vem da trilha de `HELMVEC`.
- O bloco `Sz/Tz` vem da trilha de `HELM10`.

## Ajustes em relacao a documentacao legada

- O arquivo antigo desta secao estava estruturalmente aceitavel, mas misturava traducao literal com comentarios soltos.
- Nesta revisao, a ligacao Eq. (88)-(92) foi deixada mais clara porque ela e a base conceitual para as secoes 2.2.3 e 2.2.4.
