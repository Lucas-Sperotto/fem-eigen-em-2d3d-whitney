# Plano de validação e artigo — Gemini Code Assist

**Data:** 2026-05-09  
**Auditor:** Gemini Code Assist  
**Escopo:** Análise somente-leitura do repositório para planejamento de validação e publicação.

---

## 1. Resultados existentes encontrados

A análise do repositório, com base nos relatórios de auditoria e na estrutura de arquivos, revela um conjunto substancial de artefatos:

- **CSVs de Modos e Campos:** 3253 arquivos CSV em `out/` detalhando autovalores, erros, métricas de correlação modal e campos nodais/de aresta.
- **CSVs de Validação:** `out/validation/validation_2d_22.csv` com 112 linhas cobrindo as seções 2.2.2, 2.2.3 e 2.2.4 do artigo. Outros CSVs de validação (`validation_2d_21.csv`, etc.) e os índices/veredito consolidados não estão presentes no checkout atual.
- **CSVs de Timing:** 36 arquivos `run_timing.csv` em `out/` com o tempo de `assembly`, `solve` e `post`.
- **Figuras e Imagens:** 1987 PNGs em `out/` (brutas) e 1648 PNGs em `docs/results/img/` (curadas), incluindo campos modais, figuras do artigo e comparações FEM vs. EFGMI.
- **Arquivos VTK:** 1278 arquivos em `out/` para visualização 3D dos campos modais.
- **Logs de Execução:** 38 arquivos `run.log` em `out/` com a saída textual dos executáveis.
- **Relatórios em Markdown:**
    - `docs/results/*.md`: 14 páginas de resultados curados por caso do artigo.
    - `docs/results/fem_vs_efgmi.md`: Comparativo de desempenho FEM vs. EFGMI.
    - `out/fem_efgmi_mode_report_base/FEM_EFGM_MODE_REPORT.md`: Relatório detalhado rastreado em Git comparando FEM e EFGMI para 8 casos 2D.

## 2. Classificação dos resultados

| Categoria | Exemplos | Análise |
| :--- | :--- | :--- |
| **Resultado pronto para publicação** | `docs/results/caso_01` a `09` e `11` a `14`, `docs/results/fem_vs_efgmi.md`, `out/validation/validation_2d_22.csv` (para seções 2.2.2, 2.2.3). | Evidência forte e curada, com tabelas e imagens. Pronta para ser incorporada ao artigo. |
| **Resultado bruto (útil)** | CSVs de modos e campos em `out/`, arquivos VTK, logs de execução. | Contêm todos os dados primários. Servem como fonte para gerar novas tabelas, figuras ou validações. |
| **Evidência parcial** | `docs/results/caso_10_fig13_tab10_helmvec3.md`, `out/validation/validation_2d_22.csv` (para seção 2.2.4). | A Tabela 10 está incompleta (`?`) na documentação curada, e os erros para pontos *near-cutoff* são altos (>30%), indicando um problema a ser resolvido. |
| **Evidência ausente** | Estudo de convergência de malha, validação quantitativa para `mixed-circle`/`mixed-coax`, pacote completo de validação em `out/validation/`. | São lacunas científicas. A ausência de um estudo de convergência é um bloqueio para publicação em periódicos de qualidade. |
| **Precisa ser regenerado** | `docs/results/caso_10_fig13_tab10_helmvec3.md`, índices e veredito de validação em `out/validation/`. | A Tabela 10 precisa ser corrigida no script gerador. O pacote de validação completo deve ser gerado para o fechamento. |
| **Contribuição original (não do artigo)** | Backend EFGMI e o relatório `FEM_EFGM_MODE_REPORT.md`. | É uma extensão valiosa, mas precisa ser posicionada como tal, com seu próprio estudo de convergência e justificativa. |

## 3. Experimentos mínimos para artigo

Para transformar a reprodução em um artigo publicável, os seguintes experimentos são essenciais.

### Experimento 1: Validação Cruzada de Backends
- **Objetivo:** Demonstrar que as implementações `closed-form` e `gauss` produzem resultados numericamente idênticos, validando a corretude de ambas.
- **Comando:** `./scripts/run_backend_compare.sh`
- **Entrada:** Casos de teste selecionados para cada família de solver.
- **Saída:** CSVs de modos e timing para cada backend em `out/backend_compare/`.
- **Métrica:** Diferença percentual nos autovalores entre os backends.
- **Critério de Aceite:** Diferença < 1e-9 (precisão de máquina).

### Experimento 2: Estudo de Convergência de Malha (h-refinement)
- **Objetivo:** Provar que o erro da solução FEM converge para zero com o refinamento da malha, e determinar a ordem de convergência.
- **Comando:** `python3 scripts/run_full_mesh_sweep.py --family helm10 --case rect`
- **Entrada:** Guia de onda retangular (`helm10_rect`), com malhas de `nx=6` a `nx=28`.
- **Saída:** CSV em `out/sweep/helm10_rect/` com erro vs. tamanho da malha.
- **Métrica:** Erro percentual do autovalor do modo fundamental (TE10) em relação ao valor analítico.
- **Critério de Aceite:** A ordem de convergência observada em um gráfico log-log de erro vs. `h` (tamanho do elemento) deve ser próxima de 2, como esperado para elementos lineares (P1).

### Experimento 3: Reprodução Completa das Tabelas do Artigo
- **Objetivo:** Gerar o pacote completo de validação, reproduzindo todas as tabelas numéricas do NASA TP-3485 e aplicando os gates de validação.
- **Comando:** `./scripts/build_and_run_all.sh --backend closed-form --show-validation-output`
- **Entrada:** Todos os 14 casos de teste do artigo.
- **Saída:** Artefatos em `out/` e o pacote completo em `out/validation/`, incluindo `VALIDATION_VERDICT.md`.
- **Métrica:** `err_primary_pct` e `err_secondary_pct` para cada modo, comparados com os limites em `docs/validation_thresholds.csv`.
- **Critério de Aceite:** O veredito final em `VALIDATION_VERDICT.md` deve ser `PASS` para a maioria dos casos, com justificativas para quaisquer `WARN` ou `FAIL`.

### Experimento 4: Análise de Custo Computacional (FEM vs. EFGMI)
- **Objetivo:** Comparar o custo de montagem (`assembly`) e solução (`solve`) entre o método canônico (FEM) e a extensão original (EFGMI).
- **Comando:** `python3 scripts/benchmark_backends.py` (após executar `run_backend_compare.sh` com `efgmi` habilitado).
- **Entrada:** Resultados de `out/backend_compare/`.
- **Saída:** Tabela consolidada no terminal e/ou em um arquivo Markdown.
- **Métrica:** `assembly_ms`, `solve_ms`, e a razão de tempo FEM/EFGMI.
- **Critério de Aceite:** Resultados quantitativos que demonstrem o trade-off entre os métodos (e.g., EFGMI tem montagem mais lenta mas pode ter outras vantagens).

### Experimento 5: Diagnóstico de Pontos de Instabilidade (Tabela 10)
- **Objetivo:** Investigar e resolver os altos erros nos pontos *near-cutoff* da Tabela 10.
- **Comando:** `python3 scripts/diag_224_table10_convergence.py`
- **Entrada:** Pontos problemáticos da Tabela 10 (`d/a=0.167, br/λ0=0.5`, etc.) com malhas progressivamente mais finas.
- **Saída:** CSVs e logs em `out/diag_224/`.
- **Métrica:** Erro percentual em `β/k0` em função do refinamento da malha.
- **Critério de Aceite:** O erro deve cair para dentro do gate de validação (<15%) com uma malha refinada, ou deve-se obter uma justificativa analítica para o mau condicionamento do problema nesses pontos.

## 4. Tabelas recomendadas

1.  **Tabela de Validação Principal:** Resumo dos erros primários e secundários para os 14 casos do artigo, comparando FEM, Analítico e os valores do Paper.
2.  **Tabela de Convergência de Malha:** Erro do modo TE10 (`helm10_rect`) vs. número de elementos (`N`) e ordem de convergência (`p`).
3.  **Tabela de Custo Computacional:** `Assembly Time (ms)`, `Solve Time (ms)` e `Total Time (ms)` para FEM, EFGMI e, opcionalmente, `closed-form` vs. `gauss`.
4.  **Tabela de Correção da Eq. (120):** Comparação dos autovalores `k0` obtidos com a fórmula impressa no artigo vs. a fórmula corrigida implementada no código, para o caso da Tabela 8.

## 5. Figuras recomendadas

1.  **Gráfico de Convergência:** Gráfico log-log do erro vs. `h` (inverso do número de elementos) para o Experimento 2, mostrando a inclinação que corresponde à ordem de convergência.
2.  **Gráfico de Dispersão (Figura 13):** Curva de `β/k0` vs. `br/λ0` para o caso da Tabela 10, mostrando os pontos calculados, a referência analítica e a referência do artigo.
3.  **Visualização de Modos:** Seleção de 2-3 figuras de campo modal (e.g., `helm10_rect` TE10, `fem3d0_air` TE101) para ilustrar a capacidade de pós-processamento.
4.  **Gráfico de Espectro (Espúrios):** Comparação do espectro de autovalores de `helm10` (escalar) e `helmvec` (vetorial) para uma geometria circular, mostrando a ausência de autovalores espúrios próximos de zero na formulação vetorial.

## 6. Métricas recomendadas

- **`err_primary_pct`:** Erro percentual em relação à referência analítica. Métrica principal de corretude.
- **`err_secondary_pct`:** Erro percentual em relação aos valores publicados no artigo original. Métrica de fidelidade da reprodução.
- **`rho_abs`:** Coeficiente de correlação modal. Usado para casar modos calculados com modos de referência.
- **`assembly_ms`, `solve_ms`:** Tempos de montagem e solução. Métricas de desempenho.
- **Ordem de convergência (p):** Calculada a partir do gráfico de convergência. Métrica de validação do método numérico.

## 7. Checklist de reprodutibilidade

| Item | Status | Comentário |
| :--- | :--- | :--- |
| **Dependências documentadas?** | ✔️ | O `README.md` principal e os de subdiretórios listam CMake, C++17, LAPACK, etc. |
| **Comandos de build claros?** | ✔️ | `cmake -S . -B build` e `cmake --build build` são padrão e bem documentados. |
| **Comandos de execução claros?** | ✔️ | O script `scripts/build_and_run_all.sh` é o fluxo canônico. A CLI de cada executável também é documentada. |
| **Separação de código/dados/saída?** | ✔️ | A estrutura com `src/`, `data/` (ausente mas implícito), `scripts/` e `out/` é excelente. |
| **Risco de artefatos obsoletos?** | ⚠️ | **Alto.** O repositório versiona muitos resultados em `out/fem_efgmi_mode_report_base/`. A nova política em `.gitignore` de ignorar `out/` é correta, mas exige uma decisão sobre o que fazer com os artefatos já rastreados. |
| **Scripts destrutivos?** | ⚠️ | `scripts/build_and_run.sh` apaga `build/` e `out/`, sendo perigoso. Deve ser marcado como legado ou removido. |

## 8. Riscos e inconsistências

1.  **Crítico: Tabela 10 com erros altos.** Os pontos *near-cutoff* com erro >30% invalidam a reprodução da Seção 2.2.4 se não forem resolvidos ou justificados. **Este é o principal bloqueio para publicação.**
2.  **Importante: Ausência de estudo de convergência.** Sem demonstrar que o erro diminui com o refinamento da malha, a validade dos resultados numéricos fica comprometida.
3.  **Moderado: Erro tipográfico na Eq. (120).** O código corrige um erro no artigo original, o que é uma contribuição científica valiosa. No entanto, isso precisa ser formalmente documentado no artigo, não apenas em um `README.md` de código.
4.  **Moderado: Inconsistência de artefatos versionados.** O repositório mistura a política de ignorar `out/` com o versionamento de um subdiretório (`out/fem_efgmi_mode_report_base/`). Isso pode confundir a reprodutibilidade. Recomenda-se mover os resultados EFGMI curados para `docs/results/` e limpar `out/` do controle de versão.
5.  **Menor: Tabela 10 com placeholders.** O `docs/results/caso_10_fig13_tab10_helmvec3.md` com `?` é um problema de geração de documentação que precisa ser corrigido para consistência.

## 9. Plano enxuto de execução

Para preparar o artigo, a seguinte sequência de comandos é recomendada.

```bash
# Passo 1: Corrigir o gerador da Tabela 10 (ajuste manual no script)
# vi scripts/generate_docs_tree.py

# Passo 2: Executar o fluxo completo de validação para gerar todas as evidências
# Isso irá popular 'out/' e 'out/validation/' com resultados consistentes.
./scripts/build_and_run_all.sh --backend closed-form --show-validation-output

# Passo 3: Executar o estudo de convergência de malha para o caso base
python3 scripts/run_full_mesh_sweep.py --family helm10 --case rect --backend closed-form

# Passo 4: Executar o diagnóstico específico para os pontos problemáticos da Tabela 10
python3 scripts/diag_224_table10_convergence.py

# Passo 5: Gerar o comparativo de desempenho dos backends
./scripts/run_backend_compare.sh --with-efgmi
python3 scripts/benchmark_backends.py

# Passo 6: Após a execução, analisar os novos artefatos em 'out/' e 'out/validation/'
# para escrever as seções de resultados do artigo.
ls -R out/validation/
ls out/sweep/helm10_rect/
```

Com a execução deste plano, o projeto terá todas as evidências numéricas necessárias para sustentar um artigo científico de alta qualidade, abordando as lacunas críticas de validação e convergência.

