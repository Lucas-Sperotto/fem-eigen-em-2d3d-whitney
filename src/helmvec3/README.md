# `helmvec3` - Calculo de `beta` para `k0` dado (Sec. 2.2.4)

Este modulo resolve o problema complementar ao `helmvec2`:
- entrada: `k0` (frequencia),
- saida: `beta` (constante de propagacao).

## 1) Objetivo do modulo

Montar e resolver o EVP acoplado em `beta^2` para guia retangular
parcialmente preenchido, reproduzindo os estudos da Sec. 2.2.4
(Figuras 12 e 13, Tabelas 9 e 10).

Arquivo principal:
- `main_helmvec3_rect.cpp`

Montagem compartilhada:
- `src/helmvec2/helmvec2_coupled_system.cpp`
- utilitarios `src/helmvec2/helmvec23_shared.hpp`

## 2) Sistema matematico

Com `x = [Et; Ez]`:

`P x = beta^2 Q x`

Blocos implementados (arranjo de 2.2.4):

- `P_tt = St - k0^2 Tt`
- `P_zz = k0^2 Tz`

- `Q_tt = -Mt_(1/mu)`
- `Q_tz = +C`
- `Q_zt = +C_orient^T`
- `Q_zz = Sz`

Elementos:
- `St`, `Tt`: blocos edge,
- `Sz`, `Tz`: blocos escalares,
- `Mt_(1/mu)`: massa vetorial com peso `1/mu_r`,
- `C(m,j) = int (1/mu_r) W_m . grad(N_j) dA`.

## 3) Geometrias/casos de validacao

Parametros base:
- `a = 1.0`, `b/a = 0.45`,
- `eps_fill = 2.45`.

### Figura 12 / Tabela 9

Preenchimento horizontal inferior:
- interface em `y = d = 0.5*b`.

Compara `beta/k0` para:
- `b/lambda0 = {0.2, 0.3, 0.4, 0.5, 0.6}`.

### Figura 13 / Tabela 10

Preenchimento vertical a esquerda:
- interface em `x = d`, com `d/a` variavel.

Compara blocos em funcao de `d/a` e `b/lambda0`.

## 4) Solver e filtro fisico

Solver:
- `generalized_eigs_real_vec` (`dggev`).

Filtro no pos-processamento:
- parte imaginaria pequena,
- `beta^2 > 0`,
- faixa fisica de `beta/k0` limitada por `sqrt(eps_max)` com margem numerica.

## 5) Matching modal

Estrategias implementadas:
- matching por proximidade ao valor de referencia analitica,
- opcao de rastreamento de ramo continuo (`trace_ratio_branch`) para inspecao.

Motivacao:
- reduzir erro de associacao quando ha troca de ordem entre ramos numericos.

## 6) Uso

```bash
./build/helmvec3_rect 0.20 10 5
./build/helmvec3_rect 0.20 10 5 1
# args: d_over_a nx ny [debug]
```

Saida textual em 3 blocos:
- Tabela 9 (Figura 12),
- preview de ramo (Figura 13, um `d/a`),
- validacao completa da Tabela 10.

## 7) Integracao com scripts

`scripts/validate_2d_22.py` parseia os blocos de `helmvec3_rect` e grava em:
- `build/validation_2d_22.csv`

Campos relacionados:
- `section=2.2.4`
- `case=Figure12_Table9` e `case=Figure13_Table10`
- erro relativo contra referencia analitica e HELMVEC3.

## 8) Relacao com a sequencia 2D

- `helmvec1`: `kc` no cutoff
- `helmvec2`: `k0` com `beta` dado
- `helmvec3`: `beta` com `k0` dado

A base de montagem e blocos FEM permanece a mesma; muda apenas o rearranjo
algebrico entre operador da esquerda/direita no EVP generalizado.
