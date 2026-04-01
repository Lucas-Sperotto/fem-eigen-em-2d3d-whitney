# Secao 2.2.3 - Determinacao de `k0` para `beta` Dado

## Escopo

Esta secao estende a formulacao mista da Secao 2.2.2 para resolver o numero de onda `k0` quando a constante de propagacao `beta` e conhecida. No repositorio, ela corresponde ao programa `HELMVEC2`.

## Ideia central

O artigo reorganiza a equacao vetorial completa para produzir um problema generalizado:

```text
A x = k0^2 B x
```

com estado misto:

```text
x = [E_t ; E_z]
```

## Revisao das equacoes

- Eq. (93): equacao vetorial de onda completa em 3 componentes.
- Eq. (94): dependencia longitudinal `e^{-j beta z}` do campo.
- Eq. (95): parte transversal apos separar as componentes.
- Eq. (96): parte longitudinal apos separar as componentes.
- Eq. (97): definicao de `E_t`.
- Eq. (98): reescalonamento do campo transversal para separar `beta`.
- Eq. (99): reescalonamento da componente longitudinal.
- Eq. (100): forma reescrita da parte transversal.
- Eq. (101): forma reescrita da parte longitudinal.
- Eq. (102): aplicacao de Galerkin na equacao transversal.
- Eq. (103): aplicacao de Galerkin na equacao longitudinal.
- Eq. (104): identidade vetorial para o termo `curl`.
- Eq. (105): identidade de contorno para produto vetorial.
- Eq. (106): identidade `div(f A)`.
- Eq. (107): teorema da divergencia em 2D.
- Eq. (108): forma fraca final do bloco transversal.
- Eq. (109): forma fraca final do bloco longitudinal.
- Eq. (110): equacao elemental do bloco `tt/tz`.
- Eq. (111): equacao elemental do bloco `zt/zz`.
- Eq. (112): forma matricial elemental acoplada.
- Eq. (113): definicao de `Sel(tt)`.
- Eq. (114): definicao de `Sel(tz)`.
- Eq. (115): definicao de `Sel(zt)`.
- Eq. (116): definicao de `Sel(zz)`.
- Eq. (117): definicao de `Tel(tt)`.
- Eq. (118): definicao de `Tel(zz)`.
- Eq. (119): sistema global `A x = k0^2 B x`.
- Eq. (120): forma fechada de `Sel(tt)`.
- Eq. (121): forma fechada de `Sel(tz)`.
- Eq. (122): forma fechada de `Sel(zt)`.
- Eq. (123): forma fechada de `Sel(zz)`.
- Eq. (124): forma fechada de `Tel(tt)`.
- Eq. (125): forma fechada de `Tel(zz)`.

## Resultado teorico principal

O artigo chega ao sistema:

```text
[ St + beta^2 Mt^(1/mu)     beta^2 C   ] [E_t] = k0^2 [ Tt      0 ] [E_t]
[ beta^2 C^T               beta^2 Sz  ] [E_z]        [ 0   beta^2 Tz ] [E_z]
```

Na notacao do repositorio:

```text
A_tt = St + beta^2 Mt^(1/mu)
A_tz = beta^2 C
A_zt = beta^2 C^T
A_zz = beta^2 Sz
B_tt = Tt
B_zz = beta^2 Tz
```

## Ponto de atencao teorico

A documentacao do repositorio e correta ao destacar que a Eq. (120) impressa no artigo nao deve ser copiada cegamente. A forma numericamente coerente reaproveita:

- o termo `curl-curl` da Eq. (66)
- o termo de massa vetorial da Eq. (67)
- o bloco de acoplamento `C`

Isso e importante porque a Secao 2.2.3 e exatamente o ponto onde a rastreabilidade do artigo encontra a engenharia da implementacao.

## Exemplo numerico

- Figura 11: guia quadrado parcialmente preenchido
- Tabela 8: primeiros 10 modos para `beta = 10`

## Leitura para o codigo

- Eq. (119) explica o papel central de `build_coupled_wavenumber_system_E(...)`.
- Eq. (120) a Eq. (125) justificam os helpers locais em `src/explicit/tri2d_coupled_explicit.hpp`.
- O bloco `C` e a peca conceitual que liga elementos de aresta a funcoes nodais.

## Ajustes em relacao a documentacao legada

- O arquivo antigo desta secao comecava com `2.2.2.5`, o que ja indicava mistura de conteudo.
- O mesmo arquivo tambem continha trechos inteiros de `3.1.3` a `3.1.6`, o que o tornava estruturalmente inviavel como traducao confiavel.
- Havia varias formulas quebradas por OCR, como `W*tm`, `T*t` e matrizes com barras invertidas ausentes.
