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

Na forma local do artigo:

```text
Sel(tt)  -> Eq. (120)
Sel(tz)  -> Eq. (121)
Sel(zt)  -> Eq. (122)
Sel(zz)  -> Eq. (123)
Tel(tt)  -> Eq. (124)
Tel(zz)  -> Eq. (125)
```

e, na implementacao:

```text
Sel(tt) = (1/mu_r) * curlcurl_t + beta^2 * (1/mu_r) * mass_t
Sel(tz) = beta^2 * (1/mu_r) * C
Sel(zt) = beta^2 * (1/mu_r) * C^T
Sel(zz) = beta^2 * (1/mu_r) * gradgrad_z

Tel(tt) = eps_r * mass_t
Tel(zz) = beta^2 * eps_r * mass_z
```

Na forma fraca por bloco, isso equivale a:

```text
Sel(tt)(m,n) =
    (1/mu_r) int_T curl(W_m) curl(W_n) dA
  + (beta^2/mu_r) int_T W_m . W_n dA

Sel(tz)(m,j) =
    (beta^2/mu_r) int_T W_m . grad(N_j) dA

Sel(zt)(i,n) =
    (beta^2/mu_r) int_T grad(N_i) . W_n dA

Sel(zz)(i,j) =
    (beta^2/mu_r) int_T grad(N_i) . grad(N_j) dA

Tel(tt)(m,n) =
    eps_r int_T W_m . W_n dA

Tel(zz)(i,j) =
    beta^2 eps_r int_T N_i N_j dA
```

## 2.1) Forma fechada implementada

No caminho `--backend closed-form`, o modulo deixa explicitamente nomeadas:

- `Eq. (120)` a `Eq. (125)` em
  `src/explicit/tri2d_coupled_explicit.hpp`
- o rearranjo local coerente com a `Eq. (119)` em
  `tri2d_wavenumber_rearranged_closed_form_eq_119(...)`

Assim, a trilha didatica fica:

1. `Eq. (120)` a `Eq. (125)` -> helper local explicito
2. rearranjo para a `Eq. (119)` -> helper local rearranjado
3. blocos globais nomeados -> `A_tt`, `A_tz`, `A_zt`, `A_zz`, `B_tt`, `B_zz`
4. fechamento global explicito -> `assemble_eq119_global_from_named_blocks(...)`
5. entrada publica didatica -> `tp3485::build_eq119_helmvec2_system_E(...)`
6. assembleia global subjacente -> `build_coupled_wavenumber_system_E(...)`

O ponto de montagem global fica em:

- `src/helmvec2/helmvec2_coupled_system.cpp`

Os blocos nomeados ficam preservados em:

- `CoupledWaveNumberSystem::{A_tt, A_tz, A_zt, A_zz, B_tt, B_zz}`
- arquivo `src/helmvec2/helmvec2_coupled_system.hpp`

### Observacao importante sobre a Eq. (120) do artigo

A Eq. `(120)` impressa no artigo nao fecha algebraicamente com as Eq. `(66)`,
`(67)` e `(113)` se tomada ao pe da letra. Ha duas incoerencias no texto
impresso:

1. falta o fator `beta^2` multiplicando o termo de massa vetorial
   `sum_{k=1}^5 I_tk`;
2. o termo `curl-curl`, se derivado diretamente da Eq. `(66)`, deveria trazer
   `4 D_m D_n` dentro da forma fatorada por `1/(16 A^3)`, e nao apenas
   `D_m D_n`.

Assim, a forma coerente para o bloco local `Sel(tt)` fica:

```text
Sel(tt) = (1/mu_r) * (Ltm Ltn)/(16 A^3) * (4 Dm Dn + beta^2 * sum_{k=1}^5 I_tk)
```

ou, de forma equivalente, sem colocar tudo sob o mesmo fator:

```text
Sel(tt) = (1/mu_r) * (Ltm Ltn)/(4 A^3) * Dm Dn
        + (beta^2/mu_r) * (Ltm Ltn)/(16 A^3) * sum_{k=1}^5 I_tk
```

O codigo deste repositorio ja esta consistente, porque a montagem nao depende
da Eq. `(120)` impressa isoladamente. Em vez disso, o bloco `A_tt` e montado
reaproveitando blocos locais ja validados:

- o termo `curl-curl` vem da Eq. `(66)`;
- o termo de massa vetorial vem da Eq. `(67)`;
- a combinacao global e feita como `A_tt = St + beta^2 Mt_(1/mu)`.

Em outras palavras: o artigo, nesse ponto, esta inconsistente na impressao; o
codigo esta correto porque reconstrui `Sel(tt)` a partir dos blocos elementares
ja separados e nao copia cegamente a Eq. `(120)` como aparece na pagina.

## 2.2) Bloco de acoplamento escalar-vetorial

O acoplamento local e baseado em:

```text
C(m,j) = int_T W_m . grad(N_j) dA
```

que no repositorio aparece na forma fechada:

```text
C(m,j) = (L_m / (8 A^2)) *
         [ b_j (A_m + B_m y_tri) + c_j (C_m + D_m x_tri) ]
```

Essa mesma expressao alimenta:

- Eq. `(121)` e Eq. `(122)` em `2.2.3`;
- Eq. `(138)` e Eq. `(139)` em `2.2.4`.

Como `grad(N_j)` e constante em cada triangulo linear,

```text
grad(N_j) = [b_j, c_j] / (2A),
```

o termo de acoplamento e especialmente importante para a implementacao
`closed-form`, porque permite escrever os blocos mistos sem depender de
quadratura por ponto. E justamente esse reuso que faz o codigo ficar coerente
mesmo quando a impressao da Eq. `(120)` no artigo esta incompleta.

## 3) Implementacao (`helmvec2_coupled_system.cpp`)

Pontos importantes:
- validacao de materiais por triangulo (`eps_r_tri`, `mu_r_tri`),
- contexto comum `CoupledContextE` para evitar duplicacao,
- montagem de acoplamento com quadratura triangular `P2`,
- correcoes de orientacao local/global de aresta no bloco cruzado.
- preservacao didatica dos sub-blocos globais da `Eq. (119)` antes do
  fechamento final em `A` e `B`.

Trilha de montagem global mais rastreavel:

- `initialize_named_eq119_global_blocks(...)`
- `assemble_wavenumber_system_closed_form(...)`
- `assemble_coupling_block_C_named(...)`
- `assemble_eq119_global_from_named_blocks(...)`

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
./build/helmvec2_rect 10 6 6 --debug-local-blocks --backend closed-form
./build/helmvec2_rect 10 6 6 --debug-candidates
# args: beta nx ny [debug]
```

Flags de depuracao:
- `--debug-local-blocks`: imprime os blocos locais do primeiro triangulo
  nas Eq. `(120)` a `(125)` e no rearranjo da Eq. `(119)`
- `--debug-candidates`: imprime os candidatos modais apos o filtro fisico
- `debug=1` legado: ativa os dois comportamentos ao mesmo tempo

Backends disponiveis:
- `--backend closed-form`
- `--backend gauss`

Sem `--backend`, o padrao publico do executavel e `closed-form`.

## 6.1) Saidas tipicas

O executavel agora grava em:

- `out/helmvec2/rect/run.log`
  - copia completa do terminal.
- `out/helmvec2/rect/run_timing.csv`
  - configuracao da rodada e tempos de montagem, solve e pos-processamento.
- `out/helmvec2/rect/csv/helmvec2_rect_modes.csv`
  - tabela didatica final da Figura 11 / Tabela 8, com `k0L(FEM matched)`,
    referencias, `ez_ratio`, `eig_index` casado e ponte para os artefatos
    espaciais do modo.
- `out/helmvec2/rect/csv/helmvec2_rect_candidates.csv`
  - lista dos candidatos espectrais positivos reais, com `eig_index`,
    `k0L`, `ez_ratio` e indicacao do filtro fisico.
- `out/helmvec2/rect/csv/helmvec2_rect_modeXX_candYY_Et_fields.csv`
  - campo transversal `Et` reconstruido por celula, com
    `cell_id`, `xc_m`, `yc_m`, `Ex`, `Ey`, `Emag`.
- `out/helmvec2/rect/csv/helmvec2_rect_modeXX_candYY_Ez_fields.csv`
  - componente longitudinal `Ez` reconstruida por no, com
    `node_id`, `x_m`, `y_m`, `Ez`.
- `out/helmvec2/rect/vtk/`
  - VTKs por modo casado:
    - `*_Et.vtk` para o campo transversal vetorial;
    - `*_Ez.vtk` para a componente longitudinal escalar.
- `out/helmvec2/rect/img/`
  - imagens geradas por `python3 scripts/helmvec2.py`
  - o conjunto atual inclui:
    - comparacao grafica da Figura 11 / Tabela 8
    - erro por modo
    - espectro de candidatos por rank
    - mapas de magnitude de `Et` por modo casado
    - diagramas `quiver` de `Et` por modo casado
    - mapas escalares de `Ez` por modo casado
- `out/helmvec2/rect/linop/`
  - problema global `A/B` da Eq. `(119)` em CSV no formato CRS.
  - blocos nomeados `A_tt`, `A_tz`, `A_zt`, `A_zz`, `B_tt`, `B_zz`.
  - `helmvec2_rect_eigenvalues.csv` e `helmvec2_rect_eigenvectors.csv`
    ordenados por autovalor.

Observacao:

- a tabela impressa no terminal continua estavel para o parser legado;
- os CSVs novos registram essa mesma trilha de forma estruturada;
- como a Eq. `(119)` usa `x = [Et ; Ez]`, cada modo casado exporta os dois
  componentes espaciais complementares: `Et` e `Ez`.

## 7) Integracao com scripts

`scripts/validate_2d_22.py` le este executavel e salva em:
- `out/validation/validation_2d_22.csv`

As figuras consolidadas desta secao sao geradas depois por:
- `python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22`

Para inspecao detalhada durante a validacao:

```bash
python3 scripts/validate_2d_22.py --backend closed-form --show-output --debug-local-blocks
```

Campos relacionados:
- `section=2.2.3`
- `case=Figure11_Table8`
- erro relativo contra HELMVEC2 e Hayata.

Para gerar as imagens locais do modulo:

```bash
python3 scripts/helmvec2.py
```

## 7.1) Referencias de saida

- `docs/HELMVEC2_CSV_Referencia.md`
- `docs/HELMVEC2_Imagens_Referencia.md`
- `docs/Artefatos_Espectrais_CSV_Referencia.md`
- `docs/Tabela_Executaveis_Entradas_Saidas.md`

## 8) Relacao com a sequencia 2D

- `helmvec1`: cutoff `kc` (desacoplado em blocos)
- `helmvec2`: `k0` para `beta` dado
- `helmvec3`: `beta` para `k0` dado

`helmvec2` e `helmvec3` compartilham utilitarios de perfil de material,
raizes reais e deduplicacao numerica em `helmvec23_shared.hpp`.
