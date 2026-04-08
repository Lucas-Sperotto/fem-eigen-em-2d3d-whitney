# Matriz de Casos, Executáveis e Árvore de Chamada

Este arquivo prepara a próxima etapa dos diagramas de sequência. A ideia aqui é transformar a lista dos casos numéricos do artigo em uma matriz operacional:

- qual caso do artigo cai em qual executável atual;
- qual arquivo `main_*` governa a execução;
- qual entrada pública didática, nomeada pela equação do artigo, abre o fluxo;
- qual rotina interna materializa o sistema global;
- quais arquivos sustentam a montagem, o solve e o pós-processamento;
- quantas famílias reais de diagrama existem no repositório.

Observação de arquitetura:

- os `main_*` continuam sendo os drivers reais dos executáveis;
- a entrada pública mais didática, nomeada pelas equações globais do artigo,
  agora fica em `src/article/tp3485_systems.hpp`;
- nesta matriz, as colunas de montagem distinguem explicitamente:
  - a entrada pública didática usada pelos `main_*`;
  - a rotina interna que efetivamente materializa o sistema global no solver.

## Ideia central

Os `14` casos numéricos do artigo não geram `14` árvores de chamada distintas no código atual. Eles se agrupam em `7` raízes principais de execução:

1. `HELM10` escalar 2D
2. `HELMVEC` vetorial transversal 2D
3. `HELMVEC1` misto no cutoff
4. `HELMVEC2` para `k0` com `beta` dado
5. `HELMVEC3` para `beta` com `k0` dado
6. `FEM3D0` denso
7. `FEM3D1` com montagem esparsa

Isso quer dizer que, na etapa dos diagramas, faz mais sentido desenhar primeiro por família de executável e só depois anotar quais casos do artigo entram em cada família.

## 1) Agrupamento por raiz de execução

| Raiz do diagrama | Casos do artigo | Programa histórico | Executável atual | Arquivo de entrada | Entrada didática (Eq.) | Montagem interna efetiva | Solver | Pós-processamento / comparação principal |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 1, 2, 3 | `HELM10` | `helm10_rect`, `helm10_circle`, `helm10_coax` | `src/helm10/main_helm10_*.cpp` | `tp3485::build_eq43_helm10_system(...)` | `build_helm10_scalar_system(...)` | `generalized_eigs_sym_vec(...)` | `match_*_mode_by_mass_correlation`, `helm10_post::*`, `write_vtk_unstructured_tri_scalar_vector(...)` |
| R2 | 4, 5 | `HELMVEC` | `edge_rect`, `edge_circle` | `src/helmvec/main_edge_*.cpp` | `tp3485::build_eq65_helmvec_system(...)` | `build_helm10_edge_system(...)` | `generalized_eigs_sym_vec(...)` | `match_*_edge_mode_by_mass_correlation_*`, `helmvec_post::*`, `write_vtk_unstructured_tri_cell_vector(...)` |
| R3 | 6, 7 | `HELMVEC1` | `mixed_rect`, `mixed_circle` | `src/helmvec1/main_mixed_*.cpp` | `tp3485::build_eq92_helmvec1_system_E/H(...)` | `build_system92_E(...)`, `build_system92_H(...)` | `generalized_eigs_sym_vec(...)` | `split_modes_by_block_energy(...)`, `print_rect_compare_table(...)`, `print_first_modes(...)` |
| R4 | 8 | `HELMVEC2` | `helmvec2_rect` | `src/helmvec2/main_helmvec2_rect.cpp` | `tp3485::build_eq119_helmvec2_system_E(...)` | `build_coupled_wavenumber_system_E(...)` | `generalized_eigs_real_vec(...)` | `collect_mode_candidates(...)`, `helmvec23::unique_sorted(...)`, `helmvec23::pick_closest_unused(...)` |
| R5 | 9, 10 | `HELMVEC3` | `helmvec3_fig12_rect`, `helmvec3_fig13_rect` | `src/helmvec3/main_helmvec3_fig12_rect.cpp`, `src/helmvec3/main_helmvec3_fig13_rect.cpp` | `tp3485::build_eq136_helmvec3_system_E(...)` | `build_coupled_beta_system_E(...)` | `generalized_eigs_real_vec(...)` | `beta_ratio_candidates_from_k0(...)`, `match_ratio_to_reference(...)`, `trace_ratio_branch(...)` |
| R6 | 11, 12, 13, 14 | `FEM3D0` | `fem3d0_air`, `fem3d0_half`, `fem3d0_cyl`, `fem3d0_sphere` | `src/fem3d0/main_fem3d0_air.cpp`, `src/fem3d0/main_fem3d0_half.cpp`, `src/fem3d0/main_fem3d0_cyl.cpp`, `src/fem3d0/main_fem3d0_sphere.cpp` | `tp3485::build_eq178_fem3d_system_dense(...)` | `build_helm3d_edge_system(...)` | `generalized_eigs_sym_vec(...)` | `fem3d::first_positive_k0(...)`, `fem3d::match_by_reference_with_degeneracy(...)`, `fem3d::print_table_compare(...)` |
| R7 | 11, 12, 13, 14 | `FEM3D1` | `fem3d1_air`, `fem3d1_half`, `fem3d1_cyl`, `fem3d1_sphere` | `src/fem3d1/main_fem3d1_air.cpp`, `src/fem3d1/main_fem3d1_half.cpp`, `src/fem3d1/main_fem3d1_cyl.cpp`, `src/fem3d1/main_fem3d1_sphere.cpp` | `tp3485::build_eq178_fem3d_system_sparse(...)` | `build_helm3d_edge_system_sparse(...)` | `generalized_eigs_sym_vec(...)` | `fem3d::first_positive_k0(...)`, `fem3d::match_by_reference_with_degeneracy(...)`, `fem3d::print_table_compare(...)` |

API pública didática correspondente:

- R1 -> `tp3485::build_eq43_helm10_system(...)`
- R2 -> `tp3485::build_eq65_helmvec_system(...)`
- R3 -> `tp3485::build_eq92_helmvec1_system_E/H(...)`
- R4 -> `tp3485::build_eq119_helmvec2_system_E(...)`
- R5 -> `tp3485::build_eq136_helmvec3_system_E(...)`
- R6 -> `tp3485::build_eq178_fem3d_system_dense(...)`
- R7 -> `tp3485::build_eq178_fem3d_system_sparse(...)`

## 2) Matriz caso do artigo -> execução atual

| ID | Caso do artigo | Executável atual | `main` correspondente | Parser CLI | Geração de malha / caso | Entrada didática (Eq.) | Montagem interna efetiva | Arquivos centrais da árvore |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | Retangular escalar, Tabela 1 | `helm10_rect` | `src/helm10/main_helm10_rect.cpp` | `helm10::parse_scalar_cli_options(...)` | `make_rect_mesh(...)` | `tp3485::build_eq43_helm10_system(...)` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp`, `src/core/fem_scalar_helm10.cpp`, `src/core/mode_match_rect.hpp` |
| 2 | Circular escalar, Tabela 2 | `helm10_circle` | `src/helm10/main_helm10_circle.cpp` | `helm10::parse_scalar_cli_options(...)` | `make_circle_mesh(...)` | `tp3485::build_eq43_helm10_system(...)` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp`, `src/core/mode_match_circle.hpp` |
| 3 | Coaxial escalar, Tabela 3 | `helm10_coax` | `src/helm10/main_helm10_coax.cpp` | `helm10::parse_scalar_cli_options(...)` | `make_coax_mesh(...)` | `tp3485::build_eq43_helm10_system(...)` | `build_helm10_scalar_system(...)` | `src/core/helm10_scalar_system.cpp`, `src/core/mode_match_coax.hpp` |
| 4 | Retangular edge, Tabela 4 | `edge_rect` | `src/helmvec/main_edge_rect.cpp` | `helmvec::parse_edge_cli_options(...)` | `make_rect_mesh(...)` | `tp3485::build_eq65_helmvec_system(...)` | `build_helm10_edge_system(...)` | `src/edge/edge_assembly.cpp`, `src/edge/edge_basis.cpp`, `src/edge/edge_dofs.cpp`, `src/edge/mode_match_rect_edge.hpp` |
| 5 | Circular edge, Tabela 5 | `edge_circle` | `src/helmvec/main_edge_circle.cpp` | `helmvec::parse_edge_cli_options(...)` | `make_circle_mesh(...)` | `tp3485::build_eq65_helmvec_system(...)` | `build_helm10_edge_system(...)` | `src/edge/edge_assembly.cpp`, `src/edge/mode_match_circle_edge.hpp` |
| 6 | Retangular misto, Tabela 6 | `mixed_rect` | `src/helmvec1/main_mixed_rect.cpp` | `helmvec1::parse_mixed_cli_options(...)` | `make_rect_mesh(...)` | `tp3485::build_eq92_helmvec1_system_E/H(...)` | `build_system92_E(...)`, `build_system92_H(...)` | `src/helmvec1/helmvec1_mixed_system.cpp`, `src/helmvec1/mixed_mode_utils.hpp`, `src/helmvec1/mixed_rect_reference.hpp`, `src/edge/edge_assembly.cpp`, `src/core/helm10_scalar_system.cpp` |
| 7 | Circular misto, Tabela 7 | `mixed_circle` | `src/helmvec1/main_mixed_circle.cpp` | `helmvec1::parse_mixed_cli_options(...)` | `make_circle_mesh(...)` | `tp3485::build_eq92_helmvec1_system_E/H(...)` | `build_system92_E(...)`, `build_system92_H(...)` | `src/helmvec1/helmvec1_mixed_system.cpp`, `src/helmvec1/mixed_mode_utils.hpp`, `src/edge/edge_assembly.cpp`, `src/core/helm10_scalar_system.cpp` |
| 8 | Guia quadrado parcialmente preenchido, Figura 11 / Tabela 8 | `helmvec2_rect` | `src/helmvec2/main_helmvec2_rect.cpp` | `helmvec2::parse_coupled_cli_options(...)` | `make_rect_mesh(...)` + `helmvec23::eps_step_y(...)` | `tp3485::build_eq119_helmvec2_system_E(...)` | `build_coupled_wavenumber_system_E(...)` | `src/helmvec2/helmvec2_coupled_system.cpp`, `src/helmvec2/helmvec23_shared.hpp`, `src/edge/edge_assembly.cpp`, `src/core/helm10_scalar_system.cpp` |
| 9 | Dispersão, Figura 12 / Tabela 9 | `helmvec3_fig12_rect` | `src/helmvec3/main_helmvec3_fig12_rect.cpp` | `helmvec2::parse_coupled_cli_options(...)` | `make_rect_mesh(...)` + `helmvec23::eps_step_y(...)` | `tp3485::build_eq136_helmvec3_system_E(...)` | `build_coupled_beta_system_E(...)` | `src/helmvec3/main_helmvec3_rect.cpp`, `src/helmvec3/main_helmvec3_fig12_rect.cpp`, `src/helmvec2/helmvec2_coupled_system.cpp`, `src/helmvec2/helmvec23_shared.hpp` |
| 10 | Dispersão, Figura 13 / Tabela 10 | `helmvec3_fig13_rect` | `src/helmvec3/main_helmvec3_fig13_rect.cpp` | `helmvec2::parse_coupled_cli_options(...)` | `make_rect_mesh(...)` + `helmvec23::eps_step_x(...)` | `tp3485::build_eq136_helmvec3_system_E(...)` | `build_coupled_beta_system_E(...)` | `src/helmvec3/main_helmvec3_rect.cpp`, `src/helmvec3/main_helmvec3_fig13_rect.cpp`, `src/helmvec2/helmvec2_coupled_system.cpp`, `src/helmvec2/helmvec23_shared.hpp` |
| 11 | Cavidade retangular com ar, Tabela 12 | `fem3d1_air` ou `fem3d0_air` | `src/fem3d1/main_fem3d1_air.cpp` ou `src/fem3d0/main_fem3d0_air.cpp` | `fem3d::parse_single_case_cli(...)` | `PreparedCase` | `tp3485::build_eq178_fem3d_system_sparse(...)` ou `tp3485::build_eq178_fem3d_system_dense(...)` | `build_helm3d_edge_system_sparse(...)` ou `build_helm3d_edge_system(...)` | `src/edge3d/edge3d_assembly.cpp`, `src/fem3d/fem3d_case_driver.hpp`, `src/fem3d/fem3d_reference_tables.hpp`, `src/fem3d/fem3d_compare.hpp` |
| 12 | Cavidade semi-preenchida, Tabela 13 | `fem3d1_half` ou `fem3d0_half` | `src/fem3d1/main_fem3d1_half.cpp` ou `src/fem3d0/main_fem3d0_half.cpp` | `fem3d::parse_single_case_cli(...)` | `PreparedCase` | `tp3485::build_eq178_fem3d_system_sparse(...)` ou `tp3485::build_eq178_fem3d_system_dense(...)` | `build_helm3d_edge_system_sparse(...)` ou `build_helm3d_edge_system(...)` | `src/edge3d/edge3d_assembly.cpp`, `src/fem3d/fem3d_case_driver.hpp`, `src/fem3d/fem3d_reference_tables.hpp`, `src/fem3d/fem3d_compare.hpp` |
| 13 | Cavidade cilíndrica com ar, Tabela 14 | `fem3d1_cyl` ou `fem3d0_cyl` | `src/fem3d1/main_fem3d1_cyl.cpp` ou `src/fem3d0/main_fem3d0_cyl.cpp` | `fem3d::parse_single_case_cli(...)` | `PreparedCase` | `tp3485::build_eq178_fem3d_system_sparse(...)` ou `tp3485::build_eq178_fem3d_system_dense(...)` | `build_helm3d_edge_system_sparse(...)` ou `build_helm3d_edge_system(...)` | `src/edge3d/edge3d_assembly.cpp`, `src/fem3d/fem3d_case_driver.hpp`, `src/fem3d/fem3d_reference_tables.hpp`, `src/fem3d/fem3d_compare.hpp` |
| 14 | Cavidade esférica, Tabela 15 | `fem3d1_sphere` ou `fem3d0_sphere` | `src/fem3d1/main_fem3d1_sphere.cpp` ou `src/fem3d0/main_fem3d0_sphere.cpp` | `fem3d::parse_single_case_cli(...)` | `PreparedCase` | `tp3485::build_eq178_fem3d_system_sparse(...)` ou `tp3485::build_eq178_fem3d_system_dense(...)` | `build_helm3d_edge_system_sparse(...)` ou `build_helm3d_edge_system(...)` | `src/edge3d/edge3d_assembly.cpp`, `src/fem3d/fem3d_case_driver.hpp`, `src/fem3d/fem3d_reference_tables.hpp`, `src/fem3d/fem3d_compare.hpp` |

## 3) Sequência-base de chamadas por família

### R1) Família `HELM10`

1. `helm10::parse_scalar_cli_options(...)`
2. `make_rect_mesh(...)` ou `make_circle_mesh(...)` ou `make_coax_mesh(...)`
3. `tp3485::build_eq43_helm10_system(...)`
4. `build_helm10_scalar_system(...)`
5. montagem local/global em `src/core/helm10_scalar_system.cpp` e `src/core/fem_scalar_helm10.cpp`
6. `generalized_eigs_sym_vec(...)`
7. `match_*_mode_by_mass_correlation(...)`
8. `helm10_post::extract_mode_nodal_from_Z(...)`
9. `helm10::field_reconstruction::reconstruct_transverse_fields(...)`
10. escrita de `modes.csv` e `fields_<modo>.csv`
11. `write_vtk_unstructured_tri_scalar_vector(...)`

### R2) Família `HELMVEC`

1. `helmvec::parse_edge_cli_options(...)`
2. `make_rect_mesh(...)` ou `make_circle_mesh(...)`
3. `tp3485::build_eq65_helmvec_system(...)`
4. `build_helm10_edge_system(...)`
5. montagem local/global em `src/edge/edge_assembly.cpp`, com apoio de `src/edge/edge_basis.cpp` e `src/edge/edge_dofs.cpp`
6. `generalized_eigs_sym_vec(...)`
7. `match_*_edge_mode_by_mass_correlation_*`
8. reconstrução do campo por célula em `helmvec_post::*`
9. `write_vtk_unstructured_tri_cell_vector(...)`

### R3) Família `HELMVEC1`

1. `helmvec1::parse_mixed_cli_options(...)`
2. `make_rect_mesh(...)` ou `make_circle_mesh(...)`
3. `tp3485::build_eq92_helmvec1_system_E(...)`
4. `build_system92_E(...)`
5. `generalized_eigs_sym_vec(...)`
6. `split_modes_by_block_energy(...)`
7. comparação/impressão por bloco
8. repetição do fluxo com `tp3485::build_eq92_helmvec1_system_H(...)`
9. montagem interna correspondente em `build_system92_H(...)`

Observação:
O miolo de `build_system92_E/H(...)` reutiliza dois troncos já conhecidos e,
agora, preserva explicitamente os blocos `St`, `Tt`, `Sz` e `Tz` antes de
fechar a Eq. `(92)`:

- montagem edge em `src/edge/edge_assembly.cpp`
- montagem escalar em `src/core/helm10_scalar_system.cpp`

### R4) Família `HELMVEC2`

1. `helmvec2::parse_coupled_cli_options(...)`
2. `make_rect_mesh(...)`
3. `helmvec23::eps_step_y(...)`
4. `tp3485::build_eq119_helmvec2_system_E(...)`
5. `build_coupled_wavenumber_system_E(...)`
6. montagem acoplada em `src/helmvec2/helmvec2_coupled_system.cpp`
7. `generalized_eigs_real_vec(...)`
8. `collect_mode_candidates(...)`
9. filtro físico e deduplicação em `helmvec23::*`
10. casamento final com a Tabela 8 por `helmvec23::pick_closest_unused(...)`

### R5) Família `HELMVEC3`

1. `helmvec2::parse_coupled_cli_options(...)`
2. `make_rect_mesh(...)`
3. construção do perfil material por `helmvec23::eps_step_y(...)` ou `helmvec23::eps_step_x(...)`
4. `tp3485::build_eq136_helmvec3_system_E(...)` para cada `k0` amostrado
5. `build_coupled_beta_system_E(...)`
6. `generalized_eigs_real_vec(...)`
7. `beta_ratio_candidates_from_k0(...)`
8. `match_ratio_to_reference(...)` para Tabela 9 e Tabela 10
9. `trace_ratio_branch(...)` para o ramo contínuo da Figura 13

### R6) Família `FEM3D0`

1. `fem3d::parse_cli(...)`
2. `fem3d::for_each_selected_case(...)`
3. preparação do `PreparedCase` em `src/fem3d/fem3d_case_driver.hpp`
4. `tp3485::build_eq178_fem3d_system_dense(...)`
5. `build_helm3d_edge_system(...)`
6. `build_eq178_local_tet_blocks(...)`
7. `assemble_eq178_global_generic(...)`
8. `assemble_eq178_global_dense(...)`
9. `generalized_eigs_sym_vec(...)`
10. `fem3d::first_positive_k0(...)`
11. `fem3d::match_by_reference_with_degeneracy(...)`
12. `fem3d::print_table_compare(...)`

### R7) Família `FEM3D1`

1. `fem3d::parse_cli(...)`
2. `fem3d::for_each_selected_case(...)`
3. preparação do `PreparedCase`
4. `tp3485::build_eq178_fem3d_system_sparse(...)`
5. `build_helm3d_edge_system_sparse(...)`
6. `build_eq178_local_tet_blocks(...)`
7. `assemble_eq178_global_generic(...)`
8. `assemble_eq178_global_sparse(...)`
9. conversão `SparseSymMat -> DenseMat`
10. `generalized_eigs_sym_vec(...)`
11. `fem3d::first_positive_k0(...)`
12. `fem3d::match_by_reference_with_degeneracy(...)`
13. `fem3d::print_table_compare(...)`

## 4) O que isso significa para os diagramas

Se formos produzir diagramas de sequência completos e úteis, o melhor plano é:

1. desenhar `7` diagramas-base, um para cada raiz `R1` a `R7`;
2. anexar em cada diagrama quais casos do artigo ele cobre;
3. depois, se necessário, criar subdiagramas só para trechos especiais:
   - separação TE/TM no `HELM10`
   - fluxo `E` e `H` no `HELMVEC1`
   - filtragem física em `HELMVEC2`
   - rastreamento de ramo em `HELMVEC3`
   - diferença denso vs esparso entre `FEM3D0` e `FEM3D1`

Esse caminho evita repetição e deixa a documentação mais fiel ao que o código realmente faz.
