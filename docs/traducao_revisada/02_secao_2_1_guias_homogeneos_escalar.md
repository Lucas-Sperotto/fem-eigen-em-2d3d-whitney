# Secao 2.1 - Guias Homogeneos, Formulacao Escalar

## Escopo

Esta secao trata do problema 2D homogeneo em guias de onda fechados usando potencial escalar nodal em triangulos lineares. No repositorio, ela corresponde ao programa `HELM10`.

## Ideia central

O artigo comeca pelo caso mais didatico:

- resolver o problema de Helmholtz no corte transversal
- passar da forma forte para a forma fraca
- discretizar com funcoes nodais `P1`
- montar um problema generalizado de autovalor
- recuperar campos transversais a partir do potencial

## Revisao das equacoes

- Eq. (1): equacao de Helmholtz escalar no dominio transversal.
- Eq. (2): forma integral inicial apos multiplicar pela funcao teste escalar.
- Eq. (3): reescrita do laplaciano como divergente do gradiente.
- Eq. (4): identidade vetorial usada para a integracao por partes.
- Eq. (5): teorema da divergencia em 2D no contorno `dGamma`.
- Eq. (6): forma fraca com termo de contorno ainda explicito.
- Eq. (7): forma fraca final quando as BCs TE/TM fazem o termo de contorno desaparecer.
- Eq. (8): aproximacao linear de `psi` dentro do triangulo.
- Eq. (9): valor de `psi` no no 1.
- Eq. (10): valor de `psi` no no 2.
- Eq. (11): valor de `psi` no no 3.
- Eq. (12): sistema para recuperar `a`, `b` e `c` a partir dos valores nodais.
- Eq. (13): reescrita compacta de `psi` via o vetor de coeficientes.
- Eq. (14): expansao por funcoes de forma nodais.
- Eq. (15): definicao de `alpha_i(x,y)`.
- Eq. (16): definicao de `a_i`.
- Eq. (17): definicao de `b_i`.
- Eq. (18): definicao de `c_i`.
- Eq. (19): area orientada do triangulo.
- Eq. (20): escolha de Galerkin para a funcao teste `T_s = alpha_j`.
- Eq. (21): contribuicao do lado esquerdo da forma fraca em um elemento.
- Eq. (22): contribuicao do lado direito em um elemento.
- Eq. (23): equacao elemental escalar antes da notacao matricial.
- Eq. (24): forma matricial elemental.
- Eq. (25): definicao da matriz elemental de rigidez.
- Eq. (26): definicao da matriz elemental de massa.
- Eq. (27): forma vetorial do gradiente de `alpha_i`.
- Eq. (28): gradiente explicito de `alpha_i` em funcao de `b_i` e `c_i`.
- Eq. (29): matriz de produtos internos `grad alpha_i . grad alpha_j`.
- Eq. (30): forma fechada do bloco grad-grad no triangulo.
- Eq. (31): `S_el` obtida integrando a Eq. (30) no elemento.
- Eq. (32): escrita intermediaria da matriz de massa.
- Eq. (33): matriz de massa escalar linear fechada.
- Eq. (34): problema global de autovalor `S psi = kc^2 T psi`.
- Eq. (35): reconstrucao de `psi` dentro do triangulo.
- Eq. (36): derivada de `psi` em `x`.
- Eq. (37): derivada de `psi` em `y`.
- Eq. (38): `E_x` para modo TE.
- Eq. (39): `E_y` para modo TE.
- Eq. (40): comentario operacional sobre BC de Dirichlet para TM.
- Eq. (41): `E_x` para modo TM.
- Eq. (42): `E_y` para modo TM.
- Eq. (43): fecha a trilha escalar e prepara a comparacao com tabelas/modos.

## Resultado teorico principal

O coracao da secao e a passagem:

```text
nabla_t^2 psi + kc^2 psi = 0
          ->
S_el psi_el = kc^2 T_el psi_el
          ->
S psi = kc^2 T psi
```

No repositorio, a forma fechada usada pelo backend `closed-form` e a mesma destacada no artigo:

```text
S_e(i,j) = (1/mu_r) * (b_i b_j + c_i c_j) / (4A)
T_e(i,j) = eps_r * (A/12) * [2 1 1; 1 2 1; 1 1 2]_(i,j)
```

## Figuras e tabelas desta secao

- Figura 1: geometria generica do problema
- Figura 2: triangulo linear
- Figura 3: fluxograma do `HELM10`
- Figura 4: secao retangular
- Figura 5: campos de modos retangulares
- Figura 6: secao circular
- Figura 7: campos de modos circulares
- Figura 8: linha coaxial
- Figura 9: campos de modos da linha coaxial
- Tabelas 1 a 3: validacoes retangular, circular e coaxial

## Leitura para o codigo

- Eq. (30) e Eq. (33) sao as formulas locais mais importantes para o backend fechado.
- Eq. (34) e o alvo da montagem global.
- As Eq. (38) a (42) explicam por que o repositorio consegue exportar campos a partir do potencial escalar.

## Ajustes em relacao a documentacao legada

- As legendas de figuras 3, 4, 6 e 8 foram mantidas como pendencias da traducao antiga; a trilha revisada trata essas figuras como referencias do artigo, nao como documentos completos dentro dos `.md`.
- A cadeia Eq. (30) -> Eq. (33) -> Eq. (34) foi explicitada porque ela e a ponte principal para `src/helm10`.
