# `src/explicit` - Backend de formas explicitas do artigo

Esta pasta foi idealizada para concentrar a implementacao didatica das formas
explicitas (closed-form) das matrizes elementares do artigo NASA TP-3485,
sem misturar de imediato esse trabalho com os montadores atuais baseados em
reuso de blocos e quadratura.

## 1) Objetivo da pasta

Organizar um backend explicito para:

- codificar diretamente formulas fechadas por elemento;
- manter rastreabilidade por numero de equacao do artigo;
- comparar cada forma explicita com a implementacao atual do repositorio;
- permitir validacao incremental antes de trocar qualquer fluxo principal.

Ou seja: primeiro esta pasta entra como referencia e verificacao; so depois,
se os testes fecharem, ela pode substituir partes dos montadores atuais.

## 1.1) Relacao com os programas do apendice em FORTRAN

O objetivo desta organizacao e deixar o repositorio o mais proximo possivel,
em rastreabilidade matematica, da proposta original dos programas do apendice:

- `HELM10`  -> formulacao escalar 2D (Secao 2.1, sistema global Eq. (43))
- `HELMVEC` -> formulacao vetorial transversal 2D (Secao 2.2.1, Eq. (65))
- `HELMVEC1` -> formulacao mista para `kc` (Secao 2.2.2, Eq. (92))
- `HELMVEC2` -> formulacao para `k0` dado `beta` (Secao 2.2.3, Eq. (119))
- `HELMVEC3` -> formulacao para `beta` dado `k0` (Secao 2.2.4, Eq. (136))
- `FEM3D0/FEM3D1` -> formulacao vetorial 3D (Secao 3.1, Eq. (178))

Por isso, sempre que possivel, os comentarios do codigo indicam:

- onde a matriz local closed-form e calculada;
- em qual funcao a assembleia global correspondente acontece;
- qual numero de equacao do artigo aquela etapa representa.

## 1.2) Estado atual da implementacao

Esta pasta ja contem implementacoes efetivamente usadas quando o usuario
seleciona:

- `--backend closed-form`

Arquivos ja implementados:

- `src/explicit/tri2d_scalar_explicit.hpp`
  - Eq. `(30)` e Eq. `(33)`
- `src/explicit/tri2d_edge_explicit.hpp`
  - Eq. `(66)`, Eq. `(67)` e dependencias `(73)` a `(77)`
- `src/explicit/tri2d_coupled_explicit.hpp`
  - Eq. `(120)` a `(125)`
  - Eq. `(137)` a `(142)`
  - rearranjo local coerente com a Eq. `(136)`
- `src/explicit/tet3d_edge_explicit.hpp`
  - Eq. `(181)` e Eq. `(182)`

Assim, o backend `closed-form` ja esta operacional em 2D e 3D, convivendo com
o backend por quadratura para comparacao, validacao e benchmarking.

## 2) Mapeamento das equacoes do artigo

### 2.1) Escalar 2D

- Eq. (30): forma explicita do produto `grad(a_i) . grad(a_j)` no triangulo.
- Eq. (33): matriz de massa escalar linear no triangulo.

Formas implementadas:

```text
S_e(i,j) = (b_i b_j + c_i c_j) / (4A)
T_e(i,j) = (A/12) * [2 1 1; 1 2 1; 1 1 2]_(i,j)
```

Uso esperado:
- comparar com `src/core/fem_scalar.hpp`
- eventualmente fornecer um montador escalar 100% fechado e nomeado por eq.

### 2.2) Edge 2D transversal

- Eq. (66): termo `curl-curl` local das bases de Whitney no triangulo.
- Eq. (67): termo de massa vetorial local, em funcao de `It1..It5`.

Equacoes dependentes:
- Eq. (68) a (77): definicoes de `It1..It5` e formulas de integracao.

Uso esperado:
- comparar com `src/edge/edge_assembly.cpp`
- explicitar os coeficientes `A_m`, `B_m`, `C_m`, `D_m`, `Xtri`, `Ytri`.

Formas implementadas:

```text
S_e(m,n) = (L_m L_n)/(4 A^3) * D_m D_n
T_e(m,n) = (L_m L_n)/(16 A^3) * (It1 + It2 + It3 + It4 + It5)
```

com:

```text
It1 = A_m A_n + C_m C_n
It2 = (C_m D_n + C_n D_m) x_tri
It3 = (A_m B_n + A_n B_m) y_tri
It4 = B_m B_n * (1/A) int_T y^2 dA
It5 = D_m D_n * (1/A) int_T x^2 dA
```

### 2.3) Acoplado 2.2.3 (`k0` dado `beta`)

- Eq. (120): bloco `Sel(tt)`.
- Eq. (121): bloco `Sel(tz)`.
- Eq. (122): bloco `Sel(zt)`.
- Eq. (123): bloco `Sel(zz)`.
- Eq. (124): bloco `Tel(tt)`.
- Eq. (125): bloco `Tel(zz)`.

Uso esperado:
- comparar com `src/helmvec2/helmvec2_coupled_system.cpp`
- deixar o sistema `A x = k0^2 B x` montado por formulas fechadas.

Formas implementadas:

```text
Sel(tt) = (1/mu_r) * curlcurl_t + beta^2 * (1/mu_r) * mass_t
Sel(tz) = beta^2 * (1/mu_r) * C
Sel(zt) = beta^2 * (1/mu_r) * C^T
Sel(zz) = beta^2 * (1/mu_r) * gradgrad_z
Tel(tt) = eps_r * mass_t
Tel(zz) = beta^2 * eps_r * mass_z
```

### 2.4) Acoplado 2.2.4 (`beta` dado `k0`)

- Eq. (137): bloco `Sel(tt)`.
- Eq. (138): bloco `Tel(tz)`.
- Eq. (139): bloco `Tel(zt)`.
- Eq. (140): bloco `Sel(zz)`.
- Eq. (141): bloco `Tel(tt)`.
- Eq. (142): bloco `Tel(zz)`.

Uso esperado:
- comparar com `src/helmvec2/helmvec2_coupled_system.cpp`
- deixar o sistema `P x = beta^2 Q x` montado por formulas fechadas.

Formas implementadas:

```text
Sel(tt) = (1/mu_r) * curlcurl_t - k0^2 * eps_r * mass_t
Tel(tz) = (1/mu_r) * C
Tel(zt) = (1/mu_r) * C^T
Sel(zz) = (1/mu_r) * gradgrad_z
Tel(tt) = (1/mu_r) * mass_t
Tel(zz) = (1/mu_r) * gradgrad_z - k0^2 * eps_r * mass_z
```

Implementacao atual:
- `tri2d_beta_closed_form_eq_137_142(...)`
- `tri2d_beta_rearranged_closed_form_eq_136(...)`

### 2.5) Edge 3D tetraedrico

- Eq. (181): rigidez local `Sel`.
- Eq. (182): massa vetorial local `Tel`.

Uso esperado:
- comparar com `src/edge3d/edge3d_assembly.cpp`
- substituir quadratura da massa por expressao totalmente explicita.

Formas implementadas:

```text
S_e(m,n) = (L_m L_n)/(324 V^3) * K_mn
T_e(m,n) = (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

com:

```text
K_mn = C_zm C_zn + C_xm C_xn + B_ym B_yn

I1  = A_xm A_xn + A_ym A_yn + A_zm A_zn
I2  = (A_ym B_yn + A_yn B_ym + A_zm B_zn + A_zn B_zm) x_tet
I3  = (A_xm B_xn + A_xn B_xm + A_zm C_zn + A_zn C_zm) y_tet
I4  = (A_xm C_xn + A_xn C_xm + A_ym C_yn + A_yn C_ym) z_tet
I5  = (B_zm C_zn + B_zn C_zm) * (1/V) int_T xy dV
I6  = (B_xm C_xn + B_xn C_xm) * (1/V) int_T yz dV
I7  = (B_ym C_yn + B_yn C_ym) * (1/V) int_T xz dV
I8  = (B_ym B_yn + B_zm B_zn) * (1/V) int_T x^2 dV
I9  = (B_xm B_xn + C_zm C_zn) * (1/V) int_T y^2 dV
I10 = (C_xm C_xn + C_ym C_yn) * (1/V) int_T z^2 dV
```

## 3) Estrutura sugerida da pasta

Arquivos propostos para as proximas etapas:

- `src/explicit/tri2d_scalar_explicit.hpp`
- `src/explicit/tri2d_scalar_explicit.cpp`
- `src/explicit/tri2d_edge_explicit.hpp`
- `src/explicit/tri2d_edge_explicit.cpp`
- `src/explicit/tri2d_coupled_explicit.hpp`
- `src/explicit/tri2d_coupled_explicit.cpp`
- `src/explicit/tet3d_edge_explicit.hpp`
- `src/explicit/tet3d_edge_explicit.cpp`

Motivacao da separacao:

- `tri2d_scalar_*`: Eq. (30) e Eq. (33)
- `tri2d_edge_*`: Eq. (66), Eq. (67) e dependencias
- `tri2d_coupled_*`: Eq. (120) a (125), Eq. (137) a (142)
- `tet3d_edge_*`: Eq. (181) e Eq. (182)

## 4) Estrategia de implementacao

### Fase 1 - coeficientes geometricos explicitos

Implementar estruturas auxiliares para triangulo:

- `A`, `Xtri`, `Ytri`
- `b_i`, `c_i`
- coeficientes `A_m`, `B_m`, `C_m`, `D_m`
- `L_tm`

E, para tetraedro:

- `V`
- coeficientes usados em `I1..I10`
- comprimentos `L_m`

Esta fase e a base de tudo. Sem ela, o resto fica duplicado ou opaco.

### Fase 2 - escalar 2D

Implementar primeiro Eq. (30) e Eq. (33), porque:

- sao pequenas;
- ja conhecemos o resultado esperado;
- ajudam a fixar o padrao da pasta;
- permitem validar a convencao geometrica do triangulo.

### Fase 3 - edge 2D transversal

Implementar Eq. (66), Eq. (67) e `It1..It5`.

Neste ponto, o objetivo nao e trocar o montador principal, mas comparar:

- matriz local explicita
- versus matriz local atual por formulacao existente

### Fase 4 - blocos acoplados de 2.2.3

Implementar Eq. (120) a (125) sobre o mesmo triangulo e as mesmas estruturas
de coeficientes da Fase 3.

Saida desejada:

- funcao que devolve todos os sub-blocos locais:
  - `Sel_tt`, `Sel_tz`, `Sel_zt`, `Sel_zz`
  - `Tel_tt`, `Tel_zz`

### Fase 5 - blocos acoplados de 2.2.4

Implementar Eq. (137) a (142), reaproveitando toda a geometria da Fase 3.

Aqui a validacao precisa ser dupla:

- consistencia algebraica com o artigo;
- consistencia numerica com o backend atual.

### Fase 6 - edge 3D

Implementar Eq. (181) e Eq. (182) para tetraedros.

Essa fase e mais delicada porque:

- envolve muitos coeficientes geometricos;
- a massa 3D atual esta por quadratura;
- o ganho principal sera ter formula fechada totalmente rastreavel.

## 5) Estrategia de validacao

A validacao recomendada e sempre local antes da global.

### 5.1) Validacao local

Para um unico triangulo/tetraedro:

- calcular matriz explicita;
- calcular matriz atual;
- comparar norma da diferenca;
- exigir tolerancia numerica pequena.

### 5.2) Validacao global

Depois de validado localmente:

- montar sistema global explicito;
- comparar autovalores com o montador atual;
- comparar modos em casos simples do artigo.

## 6) Regra importante de integracao

Nesta etapa, a pasta `src/explicit` nao deve alterar os executaveis principais.

A ordem correta e:

1. implementar forma explicita;
2. validar contra o backend atual;
3. criar chave de selecao (`explicit` vs `legacy`) se fizer sentido;
4. so entao considerar troca de fluxo.

## 7) Proxima etapa recomendada

A proxima codificacao deve ser:

1. `tri2d_scalar_explicit.*`
2. `tri2d_edge_explicit.*`
3. um teste/comparador local entre explicito e atual

Esse e o caminho mais seguro, mais didatico e com menor risco de regressao.
