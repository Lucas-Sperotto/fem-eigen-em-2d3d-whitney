# Rastreabilidade das Equações do Artigo no Código C++

Este arquivo passa a ser a referência central de didática e auditoria do
repositório: para cada equação importante do NASA TP-3485, ele registra
exatamente onde a forma local ou o sistema global correspondente está
codificado em C++.

## Convenção de leitura

- Equações locais do elemento devem morar, sempre que possível, em `src/explicit/`.
- Equações globais dos programas do artigo devem aparecer nas rotinas de
  montagem em `src/core/`, `src/edge/`, `src/helmvec1/`, `src/helmvec2/` e
  `src/edge3d/`.
- Fachadas públicas nomeadas pelas equações globais ficam em
  `src/article/tp3485_systems.hpp`.
- O backend `closed-form` é o caminho principal de reconstrução didática.
- O backend `gauss` permanece como caminho auxiliar de verificação numérica.

## 1) Seção 2.1 - `HELM10`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(30)` | base geométrica de `grad N_i` no triângulo linear | `explicit_tri2d::tri2d_scalar_closed_form_eq_30_33(...)` | `src/explicit/tri2d_scalar_explicit.hpp` |
| Eq. `(31)` | bloco local escalar `Sel` | `explicit_tri2d::tri2d_scalar_closed_form_eq_30_33(...)` | `src/explicit/tri2d_scalar_explicit.hpp` |
| Eq. `(33)` | bloco local escalar `Tel` | `explicit_tri2d::tri2d_scalar_closed_form_eq_30_33(...)` | `src/explicit/tri2d_scalar_explicit.hpp` |
| Eq. `(43)` | montagem interna do sistema global escalar `S u = kc^2 T u` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp` |
| Eq. `(43)` | entrada pública didática do sistema global | `tp3485::build_eq43_helm10_system(...)` | `src/article/tp3485_systems.hpp` |

Observação:

- o dispatcher local entre `closed-form` e `gauss` fica em
  `element_mats_scalar(...)`, em `src/core/fem_scalar.hpp`;
- a forma fechada principal continua sendo chamada via
  `element_mats_scalar_closed_form(...)`, também em `src/core/fem_scalar.hpp`.
- a entrada pública didática da Eq. `(43)` agora é
  `tp3485::build_eq43_helm10_system(...)`, em `src/article/tp3485_systems.hpp`.

## 2) Seção 2.2.1 - `HELMVEC`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(57)` a Eq. `(60)` | coeficientes geométricos `A_m`, `B_m`, `C_m`, `D_m` | `explicit_tri2d::tri2d_edge_closed_form_data(...)` | `src/explicit/tri2d_edge_explicit.hpp` |
| Eq. `(66)` | bloco local `curl-curl` de aresta 2D | `explicit_tri2d::tri2d_edge_curlcurl_entry_eq_66(...)` | `src/explicit/tri2d_edge_explicit.hpp` |
| Eq. `(67)` | bloco local de massa de aresta 2D | `explicit_tri2d::tri2d_edge_mass_entry_eq_67(...)` | `src/explicit/tri2d_edge_explicit.hpp` |
| Eq. `(73)` a Eq. `(77)` | termos auxiliares `It1..It5` | `explicit_tri2d::tri2d_edge_it_terms_eq_73_77(...)` | `src/explicit/tri2d_edge_explicit.hpp` |
| Eq. `(66)` e Eq. `(67)` juntas | bloco local completo `Sel/Tel` do triângulo de aresta | `explicit_tri2d::tri2d_edge_closed_form_eq_66_67(...)` | `src/explicit/tri2d_edge_explicit.hpp` |
| Eq. `(65)` | montagem interna do sistema global transversal `S x = kc^2 T x` | `build_helm10_edge_system(...)` e `assemble_edge_system_with_tri_material(...)` | `src/edge/edge_assembly.cpp` |
| Eq. `(65)` | entrada pública didática do sistema global | `tp3485::build_eq65_helmvec_system(...)` | `src/article/tp3485_systems.hpp` |

Observação:

- o dispatcher local entre `closed-form` e `gauss` fica em
  `element_mats_edge(...)`, em `src/edge/edge_assembly.cpp`.
- a entrada pública didática da Eq. `(65)` agora é
  `tp3485::build_eq65_helmvec_system(...)`, em `src/article/tp3485_systems.hpp`.

## 3) Seção 2.2.2 - `HELMVEC1`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(88)` | fonte interna do bloco vetorial transversal `St` | `build_helm10_edge_system(...)` | `src/edge/edge_assembly.cpp` |
| Eq. `(89)` | fonte interna do bloco escalar longitudinal `Sz` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp` |
| Eq. `(90)` | fonte interna do bloco vetorial de massa `Tt` | `build_helm10_edge_system(...)` | `src/edge/edge_assembly.cpp` |
| Eq. `(91)` | fonte interna do bloco escalar de massa `Tz` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp` |
| Eq. `(92)` | blocos nomeados `St`, `Tt`, `Sz`, `Tz` | `MixedSystem92::{St,Tt,Sz,Tz}` | `src/helmvec1/helmvec1_mixed_system.hpp` |
| Eq. `(92)` | montagem explícita a partir dos blocos nomeados | `load_named_eq92_blocks_from_subsystems(...)` e `assemble_eq92_global_from_named_blocks(...)` | `src/helmvec1/helmvec1_mixed_system.cpp` |
| Eq. `(92)` | entrada pública didática do sistema global | `tp3485::build_eq92_helmvec1_system_E/H(...)` | `src/article/tp3485_systems.hpp` |

Observação:

- nesta etapa do repositório, a Eq. `(92)` já está rastreável, mas ainda é
  construída por composição explícita dos blocos escalares e de aresta;
- os blocos nomeados `St`, `Tt`, `Sz` e `Tz` ficam preservados diretamente em
  `MixedSystem92`, em `src/helmvec1/helmvec1_mixed_system.hpp`;
- a entrada pública didática da Eq. `(92)` agora é
  `tp3485::build_eq92_helmvec1_system_E/H(...)`, em
  `src/article/tp3485_systems.hpp`.

## 4) Seção 2.2.3 - `HELMVEC2`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(120)` a Eq. `(125)` | seis blocos locais do problema `A x = k0^2 B x` | `explicit_tri2d::tri2d_wavenumber_closed_form_eq_120_125(...)` | `src/explicit/tri2d_coupled_explicit.hpp` |
| Eq. `(121)` e Eq. `(122)` | acoplamento local aresta-escalar | `explicit_tri2d::tri2d_edge_scalar_coupling_entry_base(...)` | `src/explicit/tri2d_coupled_explicit.hpp` |
| Eq. `(119)` | rearranjo local para o sistema global usado no código | `explicit_tri2d::tri2d_wavenumber_rearranged_closed_form_eq_119(...)` | `src/explicit/tri2d_coupled_explicit.hpp` |
| Eq. `(119)` | blocos globais nomeados `A_tt`, `A_tz`, `A_zt`, `A_zz`, `B_tt`, `B_zz` | `CoupledWaveNumberSystem::{A_tt,A_tz,A_zt,A_zz,B_tt,B_zz}` | `src/helmvec2/helmvec2_coupled_system.hpp` |
| Eq. `(119)` | montagem interna do sistema global `A x = k0^2 B x` por elemento | `assemble_wavenumber_system_closed_form(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(119)` | inicialização explícita dos sub-blocos globais | `initialize_named_eq119_global_blocks(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(119)` | fechamento explícito da matriz global a partir dos blocos nomeados | `assemble_eq119_global_from_named_blocks(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(119)` | entrada pública didática do sistema global | `tp3485::build_eq119_helmvec2_system_E(...)` | `src/article/tp3485_systems.hpp` |
| Eq. `(119)` | API interna que fecha o sistema global `HELMVEC2` | `build_coupled_wavenumber_system_E(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |

Observação importante:

- a Eq. `(120)` impressa no artigo foi validada no repositório pela
  decomposição consistente entre Eq. `(66)`, Eq. `(67)` e os termos com
  `beta^2`; a forma fechada codificada aqui segue a versão corrigida já
  documentada em `src/helmvec2/README.md` e em
  `docs/traducao/2.2.3_Determinação do número de onda.md`.
- nesta etapa, a Eq. `(119)` também fica visível no nível do código-fonte como
  blocos globais nomeados `A_tt`, `A_tz`, `A_zt`, `A_zz`, `B_tt` e `B_zz`,
  preservados em `CoupledWaveNumberSystem` antes do fechamento final de `A` e
  `B`.
- a entrada pública didática da Eq. `(119)` agora é
  `tp3485::build_eq119_helmvec2_system_E(...)`, em
  `src/article/tp3485_systems.hpp`.

## 5) Seção 2.2.4 - `HELMVEC3`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(137)` a Eq. `(142)` | seis blocos locais do problema em `beta^2` | `explicit_tri2d::tri2d_beta_closed_form_eq_137_142(...)` | `src/explicit/tri2d_coupled_explicit.hpp` |
| Eq. `(136)` | rearranjo local para `P x = beta^2 Q x` | `explicit_tri2d::tri2d_beta_rearranged_closed_form_eq_136(...)` | `src/explicit/tri2d_coupled_explicit.hpp` |
| Eq. `(136)` | blocos globais nomeados `P_tt`, `P_zz`, `Q_tt`, `Q_tz`, `Q_zt`, `Q_zz` | `CoupledBetaSystem::{P_tt,P_zz,Q_tt,Q_tz,Q_zt,Q_zz}` | `src/helmvec2/helmvec2_coupled_system.hpp` |
| Eq. `(136)` | montagem interna do sistema global `P x = beta^2 Q x` por elemento | `assemble_beta_system_closed_form(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(136)` | inicialização explícita dos sub-blocos globais | `initialize_named_eq136_global_blocks(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(136)` | fechamento explícito da matriz global a partir dos blocos nomeados | `assemble_eq136_global_from_named_blocks(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |
| Eq. `(136)` | entrada pública didática do sistema global | `tp3485::build_eq136_helmvec3_system_E(...)` | `src/article/tp3485_systems.hpp` |
| Eq. `(136)` | API interna que fecha o sistema global `HELMVEC3` | `build_coupled_beta_system_E(...)` | `src/helmvec2/helmvec2_coupled_system.cpp` |

Observação:

- embora o executável esteja em `src/helmvec3/`, a montagem global da Eq.
  `(136)` vive hoje em `src/helmvec2/helmvec2_coupled_system.cpp`; isso deve
  permanecer explicitamente documentado para evitar perda de rastreabilidade.
- nesta etapa, a Eq. `(136)` também fica visível no nível do código-fonte como
  blocos globais nomeados `P_tt`, `P_zz`, `Q_tt`, `Q_tz`, `Q_zt` e `Q_zz`,
  preservados em `CoupledBetaSystem` antes do fechamento final de `P` e `Q`.
- a entrada pública didática da Eq. `(136)` agora é
  `tp3485::build_eq136_helmvec3_system_E(...)`, em
  `src/article/tp3485_systems.hpp`.

## 6) Seção 3.1 - `FEM3D0` e `FEM3D1`

| Equação | Papel matemático | Função C++ | Arquivo |
| --- | --- | --- | --- |
| Eq. `(162)` a Eq. `(172)` | coeficientes geométricos e momentos tetraédricos | `explicit_tet3d::tet3d_edge_closed_form_data(...)` | `src/explicit/tet3d_edge_explicit.hpp` |
| Eq. `(181)` | bloco local `curl-curl` 3D | `explicit_tet3d::tet3d_edge_curlcurl_entry_eq_181(...)` | `src/explicit/tet3d_edge_explicit.hpp` |
| Eq. `(182)` | bloco local de massa vetorial 3D | `explicit_tet3d::tet3d_edge_mass_entry_eq_182(...)` | `src/explicit/tet3d_edge_explicit.hpp` |
| Eq. `(181)` e Eq. `(182)` juntas | bloco local completo `Sel/Tel` do tetraedro | `explicit_tet3d::tet3d_edge_closed_form_eq_181_182(...)` | `src/explicit/tet3d_edge_explicit.hpp` |
| Eq. `(178)` | blocos locais nomeados `Sel` e `Tel` do tetraedro | `Eq178TetLocalBlocks3D::{Sel,Tel}` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | montagem local do tetraedro a partir das Eq. `(181)` e `(182)` | `build_eq178_local_tet_blocks(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | montagem global generica compartilhada entre denso e esparso | `assemble_eq178_global_generic(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | inicializacao explicita dos operadores globais densos `S` e `T` | `initialize_eq178_dense_global_system(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | fechamento explicito do sistema global denso | `assemble_eq178_global_dense(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | inicializacao explicita dos operadores globais esparsos `S` e `T` | `initialize_eq178_sparse_global_system(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | fechamento explicito do sistema global esparso | `assemble_eq178_global_sparse(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | entrada pública didática do sistema global denso | `tp3485::build_eq178_fem3d_system_dense(...)` | `src/article/tp3485_systems.hpp` |
| Eq. `(178)` | entrada pública didática do sistema global esparso | `tp3485::build_eq178_fem3d_system_sparse(...)` | `src/article/tp3485_systems.hpp` |
| Eq. `(178)` | montagem interna do sistema global vetorial 3D denso | `build_helm3d_edge_system(...)` | `src/edge3d/edge3d_assembly.cpp` |
| Eq. `(178)` | montagem interna do sistema global vetorial 3D esparso | `build_helm3d_edge_system_sparse(...)` | `src/edge3d/edge3d_assembly.cpp` |

Observação:

- a montagem local densa e a esparsa compartilham a mesma rotina-base
  `assemble_eq178_global_generic(...)`, em `src/edge3d/edge3d_assembly.cpp`;
- nesta etapa, a Eq. `(178)` tambem fica rastreavel no proprio codigo por
  tres degraus explicitos: `build_eq178_local_tet_blocks(...)`,
  `assemble_eq178_global_generic(...)` e os fechamentos
  `assemble_eq178_global_dense(...)` / `assemble_eq178_global_sparse(...)`.
- a diferenca entre `FEM3D0` e `FEM3D1` fica reduzida explicitamente ao
  armazenamento global e ao fechamento final de `S` e `T`, nao a formulacoes
  matematicas distintas.
- as entradas públicas didáticas da Eq. `(178)` agora são
  `tp3485::build_eq178_fem3d_system_dense(...)` e
  `tp3485::build_eq178_fem3d_system_sparse(...)`, em
  `src/article/tp3485_systems.hpp`.

## 7) Regra de refatoração daqui para frente

Ao mover ou renomear qualquer rotina ligada ao artigo, este arquivo deve ser
atualizado no mesmo commit. A regra do projeto passa a ser:

1. nenhuma equação principal do artigo fica sem uma função C++ nominada;
2. nenhuma função principal fica sem ligação explícita com a numeração do
   artigo;
3. se uma equação impressa do artigo for corrigida por consistência algébrica,
   a correção deve aparecer simultaneamente:
   - na documentação teórica,
   - na documentação técnica do módulo,
   - e nesta tabela de rastreabilidade.
