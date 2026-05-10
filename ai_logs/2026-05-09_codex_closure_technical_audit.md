# Auditoria técnica de fechamento — Codex

Data da auditoria: 2026-05-09  
Escopo executado: inspeção somente-leitura do repositório, comandos de diagnóstico e criação deste relatório em `ai_logs/`.  
Observação operacional: não foram executados builds, solvers, validações ou scripts de geração nesta rodada, para não sobrescrever `out/`, `docs/results/`, `scripts/`, nem artefatos científicos existentes.

## 1. Estado Git atual

- Branch atual: `main`.
- Rastreamento: `main...origin/main`.
- Arquivos modificados rastreados:
  - `.gitignore`
  - `README.md`
  - `docs/README.md`
- Arquivos não rastreados visíveis por `git status --short`:
  - `docs/INDICE.md`
  - `docs/results/`
  - `scripts/export_presentation_assets.py`
  - `scripts/generate_docs_tree.py`
- Diff stat rastreado:
  - `.gitignore`: 30 linhas alteradas.
  - `README.md`: 63 linhas alteradas.
  - `docs/README.md`: 7 linhas alteradas.
  - Total: 81 inserções, 19 remoções.
- Natureza das mudanças detectadas:
  - Código C++/solver (`src/`): sem modificações rastreadas ou não rastreadas detectadas.
  - `include/`, `apps/`, `tests/`, `data/`: diretórios ausentes no checkout atual.
  - Scripts: dois scripts novos não rastreados em `scripts/`.
  - Documentação: `README.md`, `docs/README.md`, `docs/INDICE.md` e `docs/results/` alterados/novos.
  - Resultados/logs: `out/helm10`, `out/helmvec`, `out/helmvec1`, `out/helmvec2`, `out/helmvec3`, `out/fem3d0`, `out/fem3d1` e `out/validation` existem localmente e estão ignorados pela nova regra em `.gitignore`; `out/fem_efgmi_mode_report_base/base` tem muitos artefatos rastreados já existentes, sem diff local observado.
  - Imagens de apresentação em `docs/results/img/fem/` e `docs/results/img/efgm/` existem localmente e estão ignoradas pela nova regra em `.gitignore`.

## 2. Estrutura relevante do repositório

Diretórios principais observados:

- `src/`: implementação C++ principal, com 134 arquivos.
- `src/core/`: malhas, álgebra densa/esparsa, solve LAPACK, saída CSV/VTK, métricas e utilidades comuns.
- `src/helm10/`, `src/helmvec/`, `src/helmvec1/`, `src/helmvec2/`, `src/helmvec3/`: famílias 2D do artigo.
- `src/fem3d/`, `src/fem3d0/`, `src/fem3d1/`, `src/edge3d/`: cavidades 3D e montagem edge tetraédrica.
- `src/article/`: fachada didática por equação global do artigo.
- `src/explicit/`: blocos locais explícitos/fechados.
- `src/meshfree/`: suporte EFGMI/meshfree.
- `scripts/`: build, execução, validação, diagnóstico, comparação e plotagem.
- `docs/`: documentação científica, tradução, rastreabilidade, índices, diagramas e guias de artefatos.
- `docs/results/`: páginas curadas por caso e imagens exportadas para apresentação/artigo.
- `docs/figs/`: figuras originais do artigo já rastreadas.
- `out/`: resultados locais de execução, validação e pós-processamento.
- `build/`: build CMake local.

Não há diretórios `apps/`, `include/`, `tests/` ou `data/` no estado atual.

Executáveis definidos em `CMakeLists.txt`:

- `helm10_rect`, `helm10_circle`, `helm10_coax`, `helm10_field_reconstruction_example`
- `edge_rect`, `edge_circle`, `edge_coax`
- `mixed_rect`, `mixed_circle`, `mixed_coax`
- `helmvec2_rect`
- `helmvec3_fig12_rect`, `helmvec3_fig13_rect`
- `fem3d0_air`, `fem3d0_half`, `fem3d0_cyl`, `fem3d0_sphere`
- `fem3d1_air`, `fem3d1_half`, `fem3d1_cyl`, `fem3d1_sphere`

Scripts relevantes encontrados:

- Fluxo canônico: `scripts/build_and_run_all.sh`.
- Helper legado/local: `scripts/build_and_run.sh`; ele apaga `build/` e `out/` antes de recompilar, portanto não deve ser usado em fechamento sem decisão explícita.
- Comparação/campanhas: `scripts/run_backend_compare.sh`, `scripts/benchmark_backends.py`, `scripts/run_full_mesh_sweep.py`, `scripts/run_structured_campaign.py`.
- Validação: `scripts/validate_2d_21_csv.py`, `scripts/validate_2d_22.py`, `scripts/validate_2d_221_csv.py`, `scripts/validate_2d_222_table6_csv.py`, `scripts/validate_2d_222_table7_csv.py`, `scripts/validate_2d_224_table10_csv.py`, `scripts/validate_3d_31.py`.
- Índices/veredito de validação: `scripts/generate_validation_2d_index.py`, `scripts/generate_validation_3d_index.py`, `scripts/generate_validation_master_index.py`, `scripts/generate_validation_verdict.py`.
- Plotagem/pós-processamento: `scripts/helm10.py`, `scripts/helmvec.py`, `scripts/helmvec1.py`, `scripts/helmvec2.py`, `scripts/helmvec3.py`, `scripts/fem3d.py`, `scripts/plot_vtk_quiver.py`, `scripts/plot_validation_2d_22.py`.
- Geração de relatório/curadoria: `scripts/generate_results_md.py`, `scripts/generate_fem_efgmi_mode_report.py`, `scripts/export_curated_results.py`, `scripts/export_curated_sweep_results.py`.
- Novos e não rastreados: `scripts/export_presentation_assets.py`, `scripts/generate_docs_tree.py`.
- Diagnósticos especializados, principalmente para Seção 2.2.4/Tabela 10: `scripts/diag_224_*.py`.

Documentação científica principal:

- `docs/README.md`
- `docs/traducao/*.md`
- `docs/Rastreabilidade_Equacoes_Artigo_Codigo.md`
- `docs/Casos_de_Teste_do_Artigo.md`
- `docs/Tabela_Executaveis_Entradas_Saidas.md`
- `docs/Matriz_Casos_Executaveis_Arvore_de_Chamada.md`
- `docs/diagramas_execucao/*.md`
- Guias por família de CSV/imagens: `docs/HELM10_*`, `docs/HELMVEC_*`, `docs/HELMVEC1_*`, `docs/HELMVEC2_*`, `docs/HELMVEC3_*`, `docs/FEM3D_*`, `docs/Artefatos_Espectrais_CSV_Referencia.md`.

Resultados e evidências locais:

- `out/` contém 3253 CSVs, 1987 PNGs, 1278 VTKs, 38 `run.log` e 36 `run_timing.csv`.
- `out/validation/` contém apenas `validation_2d_22.csv` nesta fotografia.
- `docs/results/` contém 16 Markdown files: índice, 14 casos e `fem_vs_efgmi.md`.
- `docs/results/img/` contém 1648 PNGs, incluindo figuras do artigo, imagens FEM e imagens EFGMI; `docs/results/img/fem/` e `docs/results/img/efgm/` estão ignorados por `.gitignore`.

## 3. Pipeline atual de build, execução e validação

Fluxo canônico documentado e implementado:

```bash
./scripts/build_and_run_all.sh
```

O script canônico:

- configura CMake em `build/`;
- compila todos os alvos;
- executa famílias 2D (`HELM10`, `HELMVEC`, `HELMVEC1`, `HELMVEC2`, `HELMVEC3`);
- executa famílias 3D (`FEM3D0`, `FEM3D1`);
- roda validações 2D e 3D, salvo `--skip-validate`;
- gera imagens pós-processadas, salvo `--skip-images`;
- gera relatório consolidado `out/RESULTS_REPORT.md`;
- escreve log padrão em `out/run_all.log`, salvo `--no-log`.

Exemplos canônicos nomeados observados em `scripts/build_and_run_all.sh`:

```bash
./scripts/build_and_run_all.sh --backend closed-form
./scripts/build_and_run_all.sh --backend closed-form --show-validation-output
./scripts/build_and_run_all.sh --case 2.1 --backend closed-form
./scripts/build_and_run_all.sh --case 2.2.4 --backend closed-form --with-validate --with-images
./scripts/build_and_run_all.sh --profile full
```

Build direto CMake:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(nproc)"
```

Execuções representativas nomeadas usadas pelo fluxo canônico:

```bash
./build/helm10_rect --ar-m 1.0 --nx 10 --ny 20 --nmodos 10 --backend closed-form
./build/helm10_circle --nr 8 --nt 15 --nmodos 10 --backend closed-form
./build/helm10_coax --nr 10 --nt 17 --nmodos 10 --backend closed-form
./build/edge_rect --nx 10 --ny 20 --nmodos 10 --backend closed-form
./build/edge_circle --nr 8 --nt 15 --nmodos 10 --backend closed-form
./build/mixed_rect --nx 10 --ny 20 --backend closed-form
./build/mixed_circle --nr 8 --nt 15 --backend closed-form
./build/helmvec2_rect --beta 10 --nx 20 --ny 20 --backend closed-form
./build/helmvec3_fig12_rect --nx 10 --ny 5 --backend closed-form
./build/helmvec3_fig13_rect --d-over-a-preview 0.20 --nx 10 --ny 5 --backend closed-form
./build/fem3d0_air --backend closed-form
./build/fem3d1_air --backend closed-form
```

Validações:

```bash
python3 scripts/validate_2d_21_csv.py --out-root out --backend closed-form --out-csv out/validation/validation_2d_21.csv
python3 scripts/validate_2d_221_csv.py --out-root out --backend closed-form --out-csv out/validation/validation_2d_221.csv
python3 scripts/validate_2d_222_table6_csv.py --out-root out --backend closed-form --out-csv out/validation/validation_2d_222_table6.csv
python3 scripts/validate_2d_222_table7_csv.py --out-root out --backend closed-form --out-csv out/validation/validation_2d_222_table7.csv
python3 scripts/validate_2d_224_table10_csv.py --out-root out --backend closed-form --out-csv out/validation/validation_2d_224_table10.csv
python3 scripts/validate_2d_22.py --build-dir build --backend closed-form --out-csv out/validation/validation_2d_22.csv
python3 scripts/validate_3d_31.py --profile quick --solver both --build-dir build --backend closed-form --out-modes out/validation/validation_3d_31_modes.csv --out-summary out/validation/validation_3d_31_summary.csv
```

Geração de índices/veredito:

```bash
python3 scripts/generate_validation_2d_index.py --validation-dir out/validation --out-csv out/validation/validation_2d_index.csv --out-md out/validation/VALIDATION_2D_INDEX.md
python3 scripts/generate_validation_3d_index.py --validation-dir out/validation --out-csv out/validation/validation_3d_index.csv --out-md out/validation/VALIDATION_3D_INDEX.md
python3 scripts/generate_validation_master_index.py --validation-dir out/validation --out-csv out/validation/validation_master_index.csv --out-md out/validation/VALIDATION_MASTER_INDEX.md
python3 scripts/generate_validation_verdict.py --master-index out/validation/validation_master_index.csv --policy-csv docs/validation_thresholds.csv --out-csv out/validation/validation_verdict.csv --out-md out/validation/VALIDATION_VERDICT.md
```

Geração de figuras/tabelas:

```bash
python3 scripts/plot_vtk_quiver.py --all-img --build-dir build --vtk-root out --out-dir out/img_all --csv out/img_all/mode_summary.csv --mode-export 10 --max-rank 10
python3 scripts/plot_validation_2d_22.py --in-csv out/validation/validation_2d_22.csv --out-dir out/img_all/validation_2d_22 --backend closed-form
python3 scripts/generate_results_md.py --out-dir out --report out/RESULTS_REPORT.md
python3 scripts/export_presentation_assets.py
python3 scripts/generate_docs_tree.py
```

Saídas esperadas, por layout atual:

- Executáveis brutos: `out/<familia>/<caso>/run.log`, `run_timing.csv`, `csv/`, `vtk/`, `linop/`.
- Validação: `out/validation/*.csv` e `out/validation/*.md`.
- Imagens: `out/<familia>/<caso>/img/`, `out/img_all/`, `docs/results/img/`.
- Relatórios: `out/RESULTS_REPORT.md`, `out/fem_efgmi_mode_report_base/FEM_EFGM_MODE_REPORT.md`, `docs/results/*.md`.

## 4. Testes e validações encontrados

Testes unitários:

- Não há diretório `tests/`.
- Não há `enable_testing()` ou `add_test()` no `CMakeLists.txt`.
- Não foi encontrado `CTestTestfile.cmake` em `build/`.
- Conclusão: não há suíte formal de testes unitários/CTest no estado atual.

Smoke tests:

- Não há smoke tests formais separados.
- O modo seletivo do fluxo canônico (`./scripts/build_and_run_all.sh --case ...`) pode funcionar como smoke operacional, mas não aparece como suíte dedicada.
- `scripts/build_and_run.sh` é helper legado e destrutivo para `build/`/`out/`; não é adequado como smoke de fechamento sem decisão explícita.

Regressão/validação numérica:

- `validate_2d_21_csv.py`: agrega regressão CSV da Seção 2.1, Tabelas 1-3.
- `validate_2d_221_csv.py`: agrega regressão CSV da Seção 2.2.1, Tabelas 4-5.
- `validate_2d_222_table6_csv.py` e `validate_2d_222_table7_csv.py`: regressões CSV para Tabelas 6-7.
- `validate_2d_224_table10_csv.py`: regressão CSV para Figura 13 / Tabela 10.
- `validate_2d_22.py`: executa/parseia saída textual para partes de 2.2.2, 2.2.3 e 2.2.4 e gera `validation_2d_22.csv`.
- `validate_3d_31.py`: executa perfis 3D quick/full para `fem3d0`, `fem3d1` ou ambos e gera CSV por modo e resumo.
- `docs/validation_thresholds.csv`: política de gate científico por seção/tabela, com limites `pass` e `warn` para métricas primárias/secundárias.
- `generate_validation_*` e `generate_validation_verdict.py`: agregação de índices e veredito científico.

Comparação com referência analítica ou artigo original:

- Presente nos CSVs de modos e nas páginas `docs/results/*.md`.
- Métricas recorrentes: `kc_fem`, `kc_ana`, `k0_fem`, `k0_analytic`, `ref_paper`, `beta_over_k0_*`, erro percentual e `rho_abs`.
- A validação 2D existente em `out/validation/validation_2d_22.csv` compara:
  - Figura 11 / Tabela 8 contra `HELMVEC2` e Hayata;
  - Figura 12 / Tabela 9 contra referência analítica e `HELMVEC3`;
  - Figura 13 / Tabela 10 contra referência analítica e `HELMVEC3`;
  - blocos de 2.2.2, com algumas linhas sem referência primária/secundária para casos snapshot.

Fragilidades observadas:

- O conjunto formal de testes está concentrado em scripts de validação, não em testes unitários/CTest.
- `validate_2d_22.py` depende de parsing textual de stdout dos executáveis; mudanças de formatação podem quebrar a validação sem mudança numérica.
- Em `out/validation/`, nesta fotografia, só existe `validation_2d_22.csv`; não há os índices consolidados, veredito, validações 2.1/2.2.1 CSV-based nem validação 3D gerados localmente.

## 5. Resultados e evidências já úteis para artigo

Casos reproduzidos/documentados:

- `docs/results/README.md` declara reprodução dos 14 casos numéricos.
- Existem páginas para:
  - Casos 01-03: `HELM10`, Tabelas 1-3.
  - Casos 04-05: `HELMVEC`, Tabelas 4-5.
  - Casos 06-07: `HELMVEC1`, Tabelas 6-7.
  - Caso 08: `HELMVEC2`, Figura 11 / Tabela 8.
  - Casos 09-10: `HELMVEC3`, Figuras 12-13 / Tabelas 9-10.
  - Casos 11-14: `FEM3D0`/`FEM3D1`, Tabelas 12-15.

Métricas disponíveis:

- Erro percentual contra referência analítica.
- Erro percentual contra referência do artigo/paper.
- Comparação contra `HELMVEC2`, `HELMVEC3` e Hayata nos casos acoplados.
- Correlação modal `rho_abs` nos casos escalares/vetoriais/mistos.
- Timing por caso: `assembly_ms`, `solve_ms`, `post_ms`, `total_ms`.
- Autopares e matrizes em `linop/` para auditoria espectral.

Tabelas/CSVs disponíveis:

- 3253 CSVs sob `out/`.
- CSVs de modos por caso em `out/<familia>/<caso>/csv/*modes.csv`.
- CSVs de campos por modo em `out/<familia>/<caso>/csv/`.
- CSVs espectrais em `out/<familia>/<caso>/linop/`.
- Relatório FEM x EFGMI rastreado em `out/fem_efgmi_mode_report_base/FEM_EFGM_MODE_REPORT.md`.
- 2713 arquivos rastreados sob `out/`, quase todos em `out/fem_efgmi_mode_report_base/base`.

Figuras disponíveis:

- 1987 PNGs em `out/`.
- 1648 PNGs em `docs/results/img/`.
- Figuras originais do artigo em `docs/figs/` e cópias em `docs/results/img/artigo/`.
- Imagens FEM e EFGMI organizadas em `docs/results/img/fem/` e `docs/results/img/efgm/`, mas ignoradas pelo `.gitignore` atual.

Evidências numéricas fortes:

- Tabelas curadas em `docs/results/caso_01` a `caso_09` e `caso_11` a `caso_14`.
- `out/validation/validation_2d_22.csv` cobre 2.2.2, 2.2.3 e 2.2.4, incluindo a Tabela 10 em CSV.
- `docs/results/fem_vs_efgmi.md` consolida timing e erro médio por caso para FEM x EFGMI.
- `docs/Rastreabilidade_Equacoes_Artigo_Codigo.md`, `docs/Tabela_Executaveis_Entradas_Saidas.md` e `docs/Matriz_Casos_Executaveis_Arvore_de_Chamada.md` dão rastreabilidade artigo -> código -> saída.

Pendências diretamente visíveis:

- `docs/results/caso_10_fig13_tab10_helmvec3.md` contém placeholders `?` na tabela da Figura 13 / Tabela 10, apesar de `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv` conter os valores.
- A causa provável, por inspeção, está em `scripts/generate_docs_tree.py`: `_read_helmvec3_fig13_table10()` procura `beta_over_k0_fem`, `beta_over_k0_ref` e `error_percent_ref`, mas o CSV atual usa `beta_over_k0_fem_matched`, `beta_over_k0_analytic`, `beta_over_k0_helmvec3`, `error_percent_analytic` e `error_percent_helmvec3`.
- `out/RESULTS_REPORT.md` não existe nesta fotografia, embora o fluxo canônico o gere ao final.
- Índices e veredito de validação (`VALIDATION_2D_INDEX.md`, `VALIDATION_3D_INDEX.md`, `VALIDATION_MASTER_INDEX.md`, `VALIDATION_VERDICT.md`) não existem em `out/validation/` nesta fotografia.

## 6. Arquivos modificados e recomendação de commit

### Pode comitar agora

- `ai_logs/2026-05-09_codex_closure_technical_audit.md`: relatório desta rodada, criado de forma isolada e sem tocar artefatos científicos.

### Revisar antes de comitar

- `.gitignore`
  - Mudança relevante de política de versionamento: passa a ignorar `out/<familia>/`, `out/validation/`, `out/img_all/`, `docs/results/img/fem/` e `docs/results/img/efgm/`.
  - Revisar se a intenção é realmente não versionar as imagens FEM/EFGMI usadas pelas páginas `docs/results/*.md`.
- `README.md`
  - Adiciona navegação rápida e aponta para `docs/results/`.
  - Ainda há exemplos posicionais legados na seção de checklist, em tensão com a política do repositório de preferir CLI nomeada quando houver aliases.
- `docs/README.md`
  - Mudança pequena, mas parte da navegação científica.
- `docs/INDICE.md`
  - Novo índice global, revisar coerência/linkagem.
- `docs/results/*.md`
  - Curadoria valiosa, mas `caso_10_fig13_tab10_helmvec3.md` tem placeholders `?`.
- `docs/results/img/artigo/*.png` e `docs/results/img/manifest.json`
  - Parecem cópias/exportações para apresentação; revisar duplicação em relação a `docs/figs/`.
- `scripts/export_presentation_assets.py`
  - Novo script copia PNGs de `out/` para `docs/results/img/` e escreve `manifest.json`; revisar política de sobrescrita e o que deve ser versionado.
- `scripts/generate_docs_tree.py`
  - Novo gerador de `docs/results/`; revisar antes de comitar por causa do problema de colunas da Tabela 10.

### Não comitar ainda

- `docs/results/img/fem/` e `docs/results/img/efgm/`, no estado atual:
  - estão ignorados;
  - são muitos arquivos;
  - precisam de decisão editorial: versionar imagens necessárias ao artigo ou manter como artefato regenerável.
- `out/helm10/`, `out/helmvec/`, `out/helmvec1/`, `out/helmvec2/`, `out/helmvec3/`, `out/fem3d0/`, `out/fem3d1/`, `out/validation/`:
  - são artefatos de execução ignorados;
  - não devem entrar sem política explícita de dataset/artefato científico.
- `scripts/__pycache__/`:
  - artefato local; não comitar.

### Gerar novamente antes de comitar

- `docs/results/caso_10_fig13_tab10_helmvec3.md`
  - Regenerar após ajustar/revisar o mapeamento de colunas da Tabela 10.
- `docs/results/*.md`
  - Se `scripts/generate_docs_tree.py` for a fonte canônica, regenerar a árvore inteira após corrigir a Tabela 10 para evitar mistura manual/gerada.
- `out/validation/*`
  - Gerar validações completas e índices/veredito antes de usar como evidência final do artigo.
- `out/RESULTS_REPORT.md`
  - Gerar via fluxo canônico se o relatório consolidado for citado.
- `docs/results/img/manifest.json`
  - Regenerar se houver decisão de versionar ou distribuir imagens de apresentação.

## 7. Bloqueios técnicos

- Tabela 10 incompleta em documentação curada:
  - `docs/results/caso_10_fig13_tab10_helmvec3.md` contém 90 placeholders `?` em blocos de Tabela 10.
  - O CSV bruto existe e contém os dados, então o bloqueio parece ser de geração/curadoria, não de solver.
- Validação final incompleta no `out/validation/` atual:
  - Presente: `validation_2d_22.csv`.
  - Ausentes: `validation_2d_21.csv`, `validation_2d_221.csv`, `validation_2d_222_table6.csv`, `validation_2d_222_table7.csv`, `validation_2d_224_table10.csv`, validação 3D, índices 2D/3D/master e veredito.
- Sem suíte formal de testes:
  - não há `tests/`, `add_test()` ou CTest.
  - regressão depende de scripts e artefatos CSV/stdout.
- Risco de inconsistência entre documentação e política de CLI:
  - `scripts/build_and_run_all.sh` e docs novas favorecem CLI nomeada.
  - `README.md` ainda contém exemplos posicionais em seção de checklist.
- Risco de política de artefatos:
  - `docs/results/*.md` referenciam imagens em `docs/results/img/fem/` e `docs/results/img/efgm/`.
  - Esses diretórios estão ignorados; em clone limpo, as páginas podem não renderizar essas imagens se elas não forem regeneradas ou publicadas em outro canal.
- Risco de sobrescrita:
  - `scripts/generate_docs_tree.py` escreve `docs/results/*.md`.
  - `scripts/export_presentation_assets.py` escreve `docs/results/img/manifest.json` e copia PNGs.
  - Rodar esses scripts sem revisão pode sobrescrever curadoria manual.
- `scripts/build_and_run.sh` é legado e apaga `build/` e `out/`:
  - documentado como helper local, mas continua perigoso para fase de fechamento.
- `out/RESULTS_REPORT.md` ausente:
  - o fluxo canônico promete esse relatório ao final; a fotografia atual não o possui.
- 3D sem validação consolidada presente:
  - existem resultados e CSVs de modos 3D em `out/fem3d0` e `out/fem3d1`, mas não há `validation_3d_31_modes.csv` e `validation_3d_31_summary.csv` no `out/validation/` atual.

## 8. Próximas tarefas pequenas recomendadas

### 1. Fechar Tabela 10 no gerador

- Arquivos envolvidos: `scripts/generate_docs_tree.py`, `docs/results/caso_10_fig13_tab10_helmvec3.md`, `out/helmvec3/fig13_rect/csv/helmvec3_fig13_rect_table10.csv`.
- Risco: baixo para solver; médio para documentação, pois o script reescreve páginas.
- Comando de teste:

```bash
python3 scripts/generate_docs_tree.py --dry-run
rg -n "\\?" docs/results/caso_10_fig13_tab10_helmvec3.md
```

- Critério de aceite: a Tabela 10 documentada não contém `?` e usa valores de `beta_over_k0_fem_matched`, referências analítica/HELMVEC3 e erros correspondentes.

### 2. Gerar pacote completo de validação

- Arquivos envolvidos: `scripts/build_and_run_all.sh`, `scripts/validate_*.py`, `scripts/generate_validation_*.py`, `docs/validation_thresholds.csv`, `out/validation/`.
- Risco: médio, por tempo de execução e sobrescrita de `out/validation/`.
- Comando de teste:

```bash
./scripts/build_and_run_all.sh --backend closed-form --skip-images --show-validation-output
```

- Critério de aceite: existem `validation_2d_index.csv`, `validation_3d_index.csv`, `validation_master_index.csv`, `validation_verdict.csv` e `VALIDATION_VERDICT.md`, sem falhas de política científica.

### 3. Smoke canônico mínimo por família

- Arquivos envolvidos: `scripts/build_and_run_all.sh`, `CMakeLists.txt`, executáveis em `build/`, saídas em `out/<familia>/<caso>/`.
- Risco: médio, pois recompila e reexecuta casos selecionados.
- Comando de teste:

```bash
./scripts/build_and_run_all.sh --case 2.1 --backend closed-form --with-validate --skip-images
./scripts/build_and_run_all.sh --case 2.2.4 --backend closed-form --with-validate --skip-images
./scripts/build_and_run_all.sh --case table12 --backend closed-form --with-validate --skip-images
```

- Critério de aceite: cada chamada termina com status zero e produz apenas os artefatos esperados no layout `out/<familia>/<caso>/` e `out/validation/`.

### 4. Revisar exemplos CLI legados

- Arquivos envolvidos: `README.md`, `docs/Casos_de_Teste_do_Artigo.md`, `docs/Tabela_Executaveis_Entradas_Saidas.md`.
- Risco: baixo; documentação operacional.
- Comando de teste:

```bash
rg -n "\\./build/[^\\n]* [0-9]" README.md docs
```

- Critério de aceite: exemplos novos usam opções nomeadas (`--nx`, `--ny`, `--backend`, etc.); posicionais ficam marcados apenas como legado/deprecated quando mencionados.

### 5. Decidir política de imagens em `docs/results`

- Arquivos envolvidos: `.gitignore`, `docs/results/img/`, `scripts/export_presentation_assets.py`, `docs/results/*.md`.
- Risco: baixo para código; médio para reprodutibilidade visual do artigo.
- Comando de teste:

```bash
git status --ignored --short -- docs/results/img out | sed -n '1,120p'
```

- Critério de aceite: decisão explícita sobre versionar, ignorar ou publicar separadamente as imagens referenciadas por `docs/results/*.md`.

### 6. Consolidar relatório final de resultados

- Arquivos envolvidos: `scripts/generate_results_md.py`, `out/RESULTS_REPORT.md`, `README.md`, `docs/results/README.md`.
- Risco: baixo/médio, pois escreve relatório em `out/`.
- Comando de teste:

```bash
python3 scripts/generate_results_md.py --out-dir out --report out/RESULTS_REPORT.md
```

- Critério de aceite: `out/RESULTS_REPORT.md` existe, referencia o layout atual e não contradiz `docs/results/` nem o fluxo canônico.

### 7. Verificar relatório FEM x EFGMI rastreado

- Arquivos envolvidos: `out/fem_efgmi_mode_report_base/FEM_EFGM_MODE_REPORT.md`, `docs/results/fem_vs_efgmi.md`, artefatos em `out/fem_efgmi_mode_report_base/base/`.
- Risco: médio, porque `out/fem_efgmi_mode_report_base/base` já é rastreado e grande.
- Comando de teste:

```bash
git diff --stat -- out/fem_efgmi_mode_report_base
```

- Critério de aceite: relatório e CSVs rastreados estão coerentes com `docs/results/fem_vs_efgmi.md` ou a divergência é documentada.

