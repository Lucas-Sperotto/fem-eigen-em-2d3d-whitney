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

## 3) Formulacao E e formulacao H

### 3.1) `build_system92_E`

- bloco edge com `EdgeBC::TE_PEC_TangentialZero`,
- bloco escalar com `ScalarBC::TM_Dirichlet`.

Interpretacao:
- edge tende a familia TE,
- escalar tende a familia TM.

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

## 4) Classificacao modal por energia de bloco

`mixed_mode_utils.hpp` implementa:
- `split_modes_by_block_energy(...)`

Para cada autovetor:
- calcula energia em bloco 0 e bloco 1,
- classifica pelo bloco dominante,
- salva `k = sqrt(lambda)` no grupo correspondente.

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
