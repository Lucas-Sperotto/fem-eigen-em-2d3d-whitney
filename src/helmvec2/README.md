# `helmvec2` - Calculo de `k0` para `beta` dado (Sec. 2.2.3)

Este modulo resolve o problema acoplado vetorial+escalar onde a entrada e
`beta` e a saida e `k0`.

## 1) Objetivo do modulo

Para guia 2D inhomogeneo em permissividade:
- montar sistema acoplado `Et/Ez`,
- resolver EVP generalizado real nao simetrico,
- extrair raizes fisicas de `k0` e comparar com Tabela 8 (Figura 11).

Arquivos principais:
- `helmvec2_coupled_system.{hpp,cpp}`
- `main_helmvec2_rect.cpp`
- `helmvec23_shared.hpp` (utilitarios comuns com `helmvec3`)

## 2) Sistema matematico

Com `x = [Et; Ez]`:

`A x = k0^2 B x`

A forma implementada segue a reorganizacao da Sec. 2.2.3, com blocos:

- `A_tt = St + beta^2 Mt_(1/mu)`
- `A_tz = beta^2 C`
- `A_zt = beta^2 C_orient^T`
- `A_zz = beta^2 Sz`

- `B_tt = Tt`
- `B_zz = beta^2 Tz`

onde:
- `St`, `Tt`: blocos edge,
- `Sz`, `Tz`: blocos escalares,
- `Mt_(1/mu)`: massa vetorial com peso `1/mu_r`,
- `C(m,j) = int (1/mu_r) W_m . grad(N_j) dA`.

## 3) Implementacao (`helmvec2_coupled_system.cpp`)

Pontos importantes:
- validacao de materiais por triangulo (`eps_r_tri`, `mu_r_tri`),
- contexto comum `CoupledContextE` para evitar duplicacao,
- montagem de acoplamento com quadratura triangular `P2`,
- correcoes de orientacao local/global de aresta no bloco cruzado.

BCs usados:
- edge: `EdgeBC::TE_PEC_TangentialZero`
- escalar: `ScalarBC::TM_Dirichlet`

## 4) Solver e filtragem

Solver:
- `generalized_eigs_real_vec` (`LAPACKE_dggev`).

Filtragem no `main_helmvec2_rect.cpp`:
- remove autovalores com parte imaginaria relevante,
- aceita `lambda_re > 0`,
- converte `k0 = sqrt(lambda_re)`,
- aplica filtro fisico de propagacao:

`k0 > beta / sqrt(eps_max)`

## 5) Caso de validacao (Figura 11 / Tabela 8)

Configuracao padrao:
- guia quadrado `L=1`,
- `eps_top=1.0`, `eps_bottom=1.5`,
- `beta*L = 10`,
- malha `6x6` (72 triangulos).

Saida principal (formato estavel para parser):

`mode  k0L(FEM matched)  HELMVEC2(ref)  Hayata(ref)`

Matching:
- escolha gulosa por proximidade ao valor HELMVEC2 de referencia,
- sem reuso de raiz ja casada.

## 6) Uso

```bash
./build/helmvec2_rect 10 6 6
./build/helmvec2_rect 10 6 6 1
# args: beta nx ny [debug]
```

## 7) Integracao com scripts

`scripts/validate_2d_22.py` le este executavel e salva em:
- `build/validation_2d_22.csv`

Campos relacionados:
- `section=2.2.3`
- `case=Figure11_Table8`
- erro relativo contra HELMVEC2 e Hayata.

## 8) Relacao com a sequencia 2D

- `helmvec1`: cutoff `kc` (desacoplado em blocos)
- `helmvec2`: `k0` para `beta` dado
- `helmvec3`: `beta` para `k0` dado

`helmvec2` e `helmvec3` compartilham utilitarios de perfil de material,
raizes reais e deduplicacao numerica em `helmvec23_shared.hpp`.
