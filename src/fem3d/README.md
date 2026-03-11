# `fem3d` - Infra compartilhada para casos 3D (Sec. 3.1)

Este diretorio centraliza componentes comuns usados por:
- `src/fem3d0` (montagem densa),
- `src/fem3d1` (montagem esparsa simetrica).

## 1) Objetivo do modulo

Evitar duplicacao de:
- definicao de casos e geometrias,
- tabelas de referencia,
- parsing/controle de CLI para selecao de casos,
- comparacao modal FEM x analitico x referencia do paper.

## 2) Arquivos e responsabilidades

- `fem3d_reference_tables.hpp`
  - geometrias de referencia:
    - Figura 15 (`1 x 0.5 x 0.75 cm`),
    - Figura 16 (`1 x 1 x 1 cm`),
    - Figura 17 (cilindro `diam=1`, `altura=0.5`),
    - Tabela 15 (esfera `raio=1`).
  - malhas padrao (`Grid3D`) por caso,
  - tabelas de referencia (12 a 15),
  - regra de `scan_limit_for_table`.

- `fem3d_compare.hpp`
  - extracao de primeiras raizes positivas (`k0`),
  - matching com agrupamento de degenerescencia analitica,
  - impressao de tabela comparativa estavel para parser.

- `fem3d_case_driver.hpp`
  - parser de CLI comum (`--air`, `--half`, `--cyl`, `--sphere`, `--all`, `--nx --ny --nz`),
  - selecao de casos,
  - construcao de malha/material por caso (`PreparedCase`),
  - fluxo `for_each_selected_case(...)`.

## 3) Estruturas principais

- `Grid3D {nx, ny, nz}`
- `RefRow {mode, analytical, ref_paper}`
- `PreparedCase`
  - `id`, `header`, `mesh`,
  - `eps_r_tet`, `mu_r_tet`,
  - `rows` (referencias modais).

## 4) Beneficios para o projeto

- `main_fem3d0_rect.cpp` e `main_fem3d1_rect.cpp` ficam curtos e focados no solver.
- Mudancas de caso/parametro de referencia ficam em um unico lugar.
- Scripts de validacao (`validate_3d_31.py`) continuam com formato de saida
  estavel mesmo apos refatoracoes internas.

## 5) Relacao com a secao 3.1 do artigo

Este modulo nao monta matrizes nem resolve EVP diretamente.
Ele organiza a camada de reproducao dos experimentos numericos da Sec. 3.1,
onde as montagens reais estao em `src/edge3d` e os solves em `fem3d0/fem3d1`.
