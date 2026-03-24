# Resultados Curados

Esta pasta contem apenas os resultados leves e estaveis que valem a pena versionar no Git.
Os artefatos completos de `out/` continuam uteis localmente, mas nao entram todos no repositorio.

## Como atualizar esta pasta

Depois de rodar o pipeline completo em `out/`, esta curadoria pode ser
regerada com:

```bash
python3 scripts/export_curated_results.py
```

O script copia apenas arquivos selecionados de `out/` para `docs/results/` e
reconstrui este `README.md`.

## Campanha completa de sweep

Quando voce rodar a campanha pesada com `10` niveis por backend em
`out/sweeps/full_mesh_compare`, existe um exportador separado para curar
essa arvore sem subir todo o volume de artefatos:

```bash
python3 scripts/export_curated_sweep_results.py
```

Esse exportador copia para `docs/results/sweeps/full_mesh_compare/`:

- o `README.md` global da campanha
- os CSVs agregados `sweep_runs.csv`, `sweep_modes.csv` e `sweep_levels.csv`
- os CSVs de validacao do ultimo nivel disponivel de cada backend

Assim, o repositorio guarda uma visao comunicativa da campanha completa,
enquanto a arvore pesada continua local em `out/sweeps/...`.

## O que foi preservado

- CSVs-resumo de validacao 2D, 3D e benchmark
- resumo de modos 2D
- figuras representativas dos campos 2D
- figuras principais de validacao das secoes 2.2.2, 2.2.3 e 2.2.4

## Arquivos incluidos

- [csv/validation_2d_22_legacy.csv](csv/validation_2d_22_legacy.csv): Validacao 2D legado
- [csv/validation_2d_22_closed_form.csv](csv/validation_2d_22_closed_form.csv): Validacao 2D closed-form
- [csv/validation_3d_31_summary_legacy.csv](csv/validation_3d_31_summary_legacy.csv): Validacao 3D legado
- [csv/validation_3d_31_air_closed_form_summary.csv](csv/validation_3d_31_air_closed_form_summary.csv): Validacao 3D closed-form (caso air)
- [csv/backend_benchmark_summary.csv](csv/backend_benchmark_summary.csv): Benchmark de backends
- [csv/mode_summary.csv](csv/mode_summary.csv): Resumo de modos 2D
- [figs/2d_scalar_rect_te10.png](figs/2d_scalar_rect_te10.png): Campo escalar 2D retangular TE10
- [figs/2d_scalar_circle_te.png](figs/2d_scalar_circle_te.png): Campo escalar 2D circular TE
- [figs/2d_scalar_coax_te.png](figs/2d_scalar_coax_te.png): Campo escalar 2D coaxial TE
- [figs/2d_edge_rect_Et.png](figs/2d_edge_rect_Et.png): Campo edge 2D retangular Et
- [figs/2d_edge_circle_Et.png](figs/2d_edge_circle_Et.png): Campo edge 2D circular Et
- [figs/2d_edge_coax_Et.png](figs/2d_edge_coax_Et.png): Campo edge 2D coaxial Et
- [figs/2d_22_mixed_rect_cutoff.png](figs/2d_22_mixed_rect_cutoff.png): Validacao 2.2.2 retangular
- [figs/2d_22_mixed_circle_coax.png](figs/2d_22_mixed_circle_coax.png): Validacao 2.2.2 circular/coaxial
- [figs/2d_23_table8_k0L.png](figs/2d_23_table8_k0L.png): Validacao 2.2.3 tabela 8
- [figs/2d_23_table8_error_pct.png](figs/2d_23_table8_error_pct.png): Erro relativo 2.2.3 tabela 8
- [figs/2d_24_table9_beta_over_k0.png](figs/2d_24_table9_beta_over_k0.png): Validacao 2.2.4 tabela 9
- [figs/2d_24_table10_fem_branches.png](figs/2d_24_table10_fem_branches.png): Validacao 2.2.4 tabela 10

## Resumo 2D

### Validacao 2.2.x - legado

| Secao | Caso | Linhas | Max |err primary| (%) | Max |err secondary| (%) |
|---|---|---:|---:|---:|
| `2.2.2` | `mixed_circle_TE_edge` | 8 | - | - |
| `2.2.2` | `mixed_circle_TM_scalar` | 8 | - | - |
| `2.2.2` | `mixed_coax_TE_edge` | 8 | - | - |
| `2.2.2` | `mixed_coax_TM_scalar` | 8 | - | - |
| `2.2.2` | `mixed_rect_E_TE_table` | 10 | 1.528 | - |
| `2.2.2` | `mixed_rect_E_TM_table` | 10 | 13.547 | - |
| `2.2.3` | `Figure11_Table8` | 10 | 0.077 | 1.317 |
| `2.2.4` | `Figure12_Table9` | 5 | 9.208 | 7.277 |
| `2.2.4` | `Figure13_Table10` | 45 | 37.676 | 59.580 |

### Validacao 2.2.x - closed-form

| Secao | Caso | Linhas | Max |err primary| (%) | Max |err secondary| (%) |
|---|---|---:|---:|---:|
| `2.2.2` | `mixed_circle_TE_edge` | 8 | - | - |
| `2.2.2` | `mixed_circle_TM_scalar` | 8 | - | - |
| `2.2.2` | `mixed_coax_TE_edge` | 8 | - | - |
| `2.2.2` | `mixed_coax_TM_scalar` | 8 | - | - |
| `2.2.2` | `mixed_rect_E_TE_table` | 10 | 1.528 | - |
| `2.2.2` | `mixed_rect_E_TM_table` | 10 | 13.547 | - |
| `2.2.3` | `Figure11_Table8` | 10 | 0.077 | 1.317 |
| `2.2.4` | `Figure12_Table9` | 5 | 9.208 | 7.277 |
| `2.2.4` | `Figure13_Table10` | 45 | 37.676 | 59.580 |

### Figuras 2D escolhidas

- Escalar retangular: [figs/2d_scalar_rect_te10.png](figs/2d_scalar_rect_te10.png)
- Escalar circular: [figs/2d_scalar_circle_te.png](figs/2d_scalar_circle_te.png)
- Escalar coaxial: [figs/2d_scalar_coax_te.png](figs/2d_scalar_coax_te.png)
- Edge retangular: [figs/2d_edge_rect_Et.png](figs/2d_edge_rect_Et.png)
- Edge circular: [figs/2d_edge_circle_Et.png](figs/2d_edge_circle_Et.png)
- Edge coaxial: [figs/2d_edge_coax_Et.png](figs/2d_edge_coax_Et.png)
- Secao 2.2.2: [figs/2d_22_mixed_rect_cutoff.png](figs/2d_22_mixed_rect_cutoff.png), [figs/2d_22_mixed_circle_coax.png](figs/2d_22_mixed_circle_coax.png)
- Secao 2.2.3: [figs/2d_23_table8_k0L.png](figs/2d_23_table8_k0L.png), [figs/2d_23_table8_error_pct.png](figs/2d_23_table8_error_pct.png)
- Secao 2.2.4: [figs/2d_24_table9_beta_over_k0.png](figs/2d_24_table9_beta_over_k0.png), [figs/2d_24_table10_fem_branches.png](figs/2d_24_table10_fem_branches.png)

## Resumo 3D

### Melhor erro medio por solver/caso - legado

| Solver | Caso | Configs | Melhor mean err ana (%) | Melhor mean err ref (%) |
|---|---|---:|---:|---:|
| `fem3d0` | `air` | 3 | 0.955 | 1.531 |
| `fem3d0` | `cyl` | 3 | 2.152 | 1.777 |
| `fem3d0` | `half` | 3 | 0.312 | 0.613 |
| `fem3d0` | `sphere` | 3 | 2.328 | 1.884 |
| `fem3d1` | `air` | 3 | 0.955 | 1.531 |
| `fem3d1` | `cyl` | 3 | 2.152 | 1.777 |
| `fem3d1` | `half` | 3 | 0.312 | 0.613 |
| `fem3d1` | `sphere` | 3 | 2.328 | 1.884 |

### Closed-form 3D preservado

Nesta curadoria, foi preservado o resumo `closed-form` do caso `air`, suficiente para documentar a trilha de comparacao sem subir todos os artefatos intermediarios.

- [csv/validation_3d_31_air_closed_form_summary.csv](csv/validation_3d_31_air_closed_form_summary.csv)

## Benchmark de backends

| Caso | total_ms_mean gauss | total_ms_mean closed-form | Mais rapido |
|---|---:|---:|---|
| `edge_rect` | 1121.341 | 1428.645 | `gauss` |
| `fem3d0_air` | 1145.094 | 1074.070 | `closed-form` |
| `fem3d1_air` | 1053.529 | 1215.174 | `gauss` |
| `helm10_rect` | 504.984 | 252.696 | `closed-form` |
| `helmvec2_rect` | 30.407 | 26.853 | `closed-form` |
| `helmvec3_rect` | 2794.369 | 3822.492 | `gauss` |
| `mixed_rect` | 140.050 | 105.695 | `closed-form` |

## O que ficou de fora do Git

- pasta `out/` completa
- todos os `.vtk`
- logs integrais de execucao
- imagens repetidas por rank/modo
- CSVs detalhados demais para leitura humana

A ideia aqui e simples: o repositório guarda o que comunica bem os resultados; os artefatos pesados continuam regeneraveis localmente.
