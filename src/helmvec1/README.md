# `helmvec1` - Sistema misto vetorial+escalar para `kc` (Sec. 2.2.2)

Este modulo implementa a formulacao de 3 componentes no cutoff,
combinando:
- base de aresta para campo transversal,
- base nodal `P1` para componente longitudinal.

## 1) Objetivo do modulo

Resolver o problema em `kc` com estado misto:

`x = [transversal_edge ; longitudinal_scalar]`

em geometrias 2D padrao (retangular, circular, coaxial),
com comparacoes tabuladas e separacao modal por energia de bloco.

Arquivos principais:
- `helmvec1_mixed_system.{hpp,cpp}`
- `main_mixed_rect.cpp`
- `main_mixed_circle.cpp`
- `main_mixed_coax.cpp`
- `mixed_mode_utils.hpp`
- `mixed_rect_reference.hpp`

## 2) Sistema global (Eq. 92)

No cutoff (`beta = 0`), a versao implementada usa blocos diagonais:

`S = block_diag(St, Sz)`

`T = block_diag(Tt, Tz)`

`S x = kc^2 T x`

onde:
- `St`, `Tt` vem da montagem edge,
- `Sz`, `Tz` vem da montagem escalar.

Mais explicitamente:

```text
St(m,n) = (1/mu_r) int_T curl(W_m) curl(W_n) dA
Tt(m,n) = eps_r int_T W_m . W_n dA

Sz(i,j) = (1/mu_r) int_T grad(N_i) . grad(N_j) dA
Tz(i,j) = eps_r int_T N_i N_j dA
```

No problema global:

```text
[ St   0 ] [et] = kc^2 [ Tt   0 ] [et]
[  0  Sz ] [ez]        [  0  Tz ] [ez]
```

No backend `closed-form`, esses blocos sao construidos pelas expressoes:

```text
St(m,n) = (1/mu_r) * (L_m L_n)/(4 A^3) * D_m D_n
Tt(m,n) = eps_r * (L_m L_n)/(16 A^3) * (It1 + It2 + It3 + It4 + It5)

Sz(i,j) = (1/mu_r) * (b_i b_j + c_i c_j) / (4A)
Tz(i,j) = eps_r * (A/12) * [2 1 1; 1 2 1; 1 1 2]_(i,j)
```

Ou seja, a Eq. `(92)` e montada apenas por justaposicao de dois problemas ja
presentes nas secoes anteriores:

- o problema edge transversal da Eq. `(65)`;
- o problema escalar nodal da Eq. `(43)`.

## 2.1) Trilha de rastreabilidade

Para conectar o artigo ao codigo, a trilha principal deste modulo e:

1. Eq. `(66)` e Eq. `(67)` -> bloco vetorial transversal em
   `src/explicit/tri2d_edge_explicit.hpp`
2. Eq. `(30)` e Eq. `(33)` -> bloco escalar longitudinal em
   `src/explicit/tri2d_scalar_explicit.hpp`
3. Eq. `(92)` -> blocos nomeados `St`, `Tt`, `Sz`, `Tz` preservados em
   `MixedSystem92`, em `src/helmvec1/helmvec1_mixed_system.hpp`
4. entrada publica didatica da Eq. `(92)` ->
   `tp3485::build_eq92_helmvec1_system_E/H(...)` em
   `src/article/tp3485_systems.hpp`
5. Eq. `(92)` -> sistema global misto montado em
   `src/helmvec1/helmvec1_mixed_system.cpp`
6. funcoes de montagem principais subjacentes ->
   `build_system92_E(...)` e `build_system92_H(...)`

## 3) Formulacao E e formulacao H

### 3.1) `build_system92_E`

- bloco edge com `EdgeBC::TE_PEC_TangentialZero`,
- bloco escalar com `ScalarBC::TM_Dirichlet`.

Interpretacao:
- edge tende a familia TE,
- escalar tende a familia TM.

Em termos de operador:

```text
St <- bloco vetorial transversal com PEC tangencial
Sz <- bloco escalar longitudinal com Dirichlet
```

### 3.2) `build_system92_H`

Operador dual por troca constitutiva:
- `eps_proxy <- mu`,
- `mu_proxy <- eps`.

BCs:
- edge: `EdgeBC::TM_PEC_NormalZero`,
- escalar: `ScalarBC::TE_Neumann`.

Interpretacao:
- edge tende a familia TM,
- escalar tende a familia TE.

Na forma dual:

```text
eps_proxy <- mu_r
mu_proxy  <- eps_r
```

e os blocos passam a ser montados com essa troca constitutiva:

```text
St^H(m,n) = (1/eps_r) int_T curl(W_m) curl(W_n) dA
Tt^H(m,n) = mu_r int_T W_m . W_n dA

Sz^H(i,j) = (1/eps_r) int_T grad(N_i) . grad(N_j) dA
Tz^H(i,j) = mu_r int_T N_i N_j dA
```

Assim, o repositorio consegue representar a formulacao `H` sem duplicar toda a
infraestrutura geometrica: muda apenas o peso material e o tipo de BC.

## 4) Classificacao modal por energia de bloco

`mixed_mode_utils.hpp` implementa:
- `split_modes_by_block_energy(...)`

Para cada autovetor:
- calcula energia em bloco 0 e bloco 1,
- classifica pelo bloco dominante,
- salva `k = sqrt(lambda)` no grupo correspondente.

Ou seja, para `x = [x0; x1]`, o criterio numerico e:

```text
E0 = ||x0||^2
E1 = ||x1||^2
```

e o modo vai para o bloco com maior energia.

Esse criterio evita ambiguidade de ordenacao quando o solve global mistura
modos de familias diferentes.

## 5) Referencias analiticas (retangular)

`mixed_rect_reference.hpp` centraliza:
- geracao de familias analiticas TE/TM para retangulo,
- impressao de tabela de comparacao no formato esperado pelos scripts.

Formula de cutoff:

`kc(m,n) = sqrt((m*pi/a)^2 + (n*pi/b)^2)`

regras:
- TE: `(m,n) != (0,0)`
- TM: `m>=1`, `n>=1`

## 6) Solver

Como os sistemas montados sao simetricos:
- `generalized_eigs_sym_vec` (`dsygv`).

## 7) Uso

Retangular:

```bash
./build/mixed_rect 12 6
```

Circular:

```bash
./build/mixed_circle 10 48
```

Coaxial:

```bash
./build/mixed_coax 10 48
```

## 7.1) Backends e depuracao

Flags disponiveis:

- `--backend closed-form`
- `--backend gauss`
- `--debug-local-blocks`
- `--debug-candidates`
- `--debug` ou `--debug-all`

Sem `--backend`, o padrao publico do executavel e `closed-form`.

Exemplos:

```bash
./build/mixed_rect 12 6 --backend closed-form --debug-local-blocks --debug-candidates
./build/mixed_circle 10 48 --debug-candidates
./build/mixed_coax 10 48 --backend closed-form --debug-local-blocks
```

Interpretacao:

- `--debug-local-blocks`: imprime os blocos locais edge e escalar do primeiro
  triangulo que alimentam a Eq. `(92)`;
- `--debug-candidates`: imprime os candidatos separados por bloco dominante
  antes da tabela final.
- os blocos globais nomeados `St`, `Tt`, `Sz` e `Tz` permanecem disponiveis em
  `MixedSystem92` para leitura didatica direta da Eq. `(92)`.

## 8) Saidas tipicas

- retangular: tabelas TE/TM para formulacao E e H
- circular/coaxial: snapshot de primeiros modos por bloco dominante

As strings dos blocos no `main_mixed_rect` sao mantidas estaveis para parse do
script `scripts/validate_2d_22.py`.

## 9) Relacao com o artigo

Corresponde a Sec. 2.2.2:
- formulacao mista no cutoff,
- base para os sistemas acoplados de `helmvec2` e `helmvec3`.

No repositorio, o modulo ja sai preparado para:
- comparacao automatica,
- documentacao de familias modais,
- extensao para inhomogeneidade por triangulo.
