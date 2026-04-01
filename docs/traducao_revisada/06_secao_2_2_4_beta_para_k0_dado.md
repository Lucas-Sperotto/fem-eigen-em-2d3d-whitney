# Secao 2.2.4 - Caracteristicas de Dispersao, `beta` para `k0` Dado

## Escopo

Esta secao e a extensao complementar da Secao 2.2.3. Em vez de resolver `k0` para `beta` dado, ela rearranja o problema para resolver `beta` quando `k0` e conhecido. No repositorio, corresponde ao programa `HELMVEC3`.

## Ideia central

O artigo nao troca a discretizacao. O que muda e o rearranjo algebrico:

```text
P x = beta^2 Q x
```

com o mesmo estado misto:

```text
x = [E_t ; E_z]
```

## Revisao das equacoes

- Eq. (126): lado esquerdo reorganizado para o bloco transversal.
- Eq. (127): lado esquerdo reorganizado para o bloco longitudinal.
- Eq. (128): equacao elemental do bloco transversal no novo rearranjo.
- Eq. (129): equacao elemental do bloco longitudinal no novo rearranjo.
- Eq. (130): forma matricial elemental da secao 2.2.4.
- Eq. (131): bloco elemental `Sel(tt)`.
- Eq. (132): bloco elemental `Tel(tt)`.
- Eq. (133): bloco elemental `Tel(tz)`.
- Eq. (134): bloco elemental `Tel(zt)`.
- Eq. (135): bloco elemental `Tel(zz)`.
- Eq. (136): sistema global `P x = beta^2 Q x`.
- Eq. (137): forma fechada de `Sel(tt)`.
- Eq. (138): forma fechada de `Tel(tz)`.
- Eq. (139): forma fechada de `Tel(zt)`.
- Eq. (140): forma fechada do termo grad-grad escalar.
- Eq. (141): forma fechada de `Tel(tt)`.
- Eq. (142): forma fechada de `Tel(zz)`.

## Resultado teorico principal

O rearranjo usado pelo repositorio pode ser lido assim:

```text
P_tt = St - k0^2 Tt
P_zz = k0^2 Tz

Q_tt = -Mt^(1/mu)
Q_tz = C
Q_zt = C^T
Q_zz = Sz
```

Isso explica por que a implementacao de `HELMVEC3` reaproveita praticamente toda a infraestrutura de `HELMVEC2`.

## Exemplo numerico

- Figura 12 / Tabela 9:
  - guia retangular parcialmente preenchido
  - caso com `b_r / a_r = 0.45`, `d / b_r = 0.5` e `eps_r = 2.45`
- Figura 13 / Tabela 10:
  - segundo exemplo de guia parcialmente preenchido
  - variacao com `d / a_r`

## O que esta secao ensina

- a discretizacao nao muda entre 2.2.3 e 2.2.4
- o que muda e a distribuicao dos termos entre os operadores da esquerda e da direita
- a forma local fechada continua reaproveitando `St`, `Tt`, `Sz`, `Tz` e `C`

## Leitura para o codigo

- Eq. (136) e a justificativa matematica de `build_coupled_beta_system_E(...)`.
- Eq. (137) a Eq. (142) justificam os helpers locais fechados usados no backend `closed-form`.
- A comparacao com a Secao 2.2.3 e essencial para entender por que `src/helmvec2/helmvec2_coupled_system.cpp` serve aos dois modulos.

## Ajustes em relacao a documentacao legada

- O arquivo antigo desta secao estava menos contaminado do que o de 2.2.3, mas sofria com varias deformacoes de OCR:
  - `E*t`, `W*tm`, `T*{tm}`
  - matrizes sem quebras de linha validas
  - formulas com barras invertidas ausentes
- A Eq. (142) estava escrita de forma pouco legivel; nesta trilha ela foi tratada conceitualmente como a soma do termo grad-grad e do termo `k0^2 eps_r M_z`.
