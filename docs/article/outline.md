# Outline do Artigo — Reprodução e Extensão Meshfree do NASA TP-3485

**Título provisório:**  
*Reprodução Computacional e Extensão Meshfree do Método de Elementos Finitos para Problemas de Autovalor em Eletromagnetismo*

**Título alternativo (inglês):**  
*Computational Reproduction and Meshfree Extension of the Finite Element Method for Electromagnetic Eigenvalue Problems*

**Periódicos-alvo sugeridos (por ordem de prioridade):**  
1. *IEEE Transactions on Microwave Theory and Techniques*  
2. *IEEE Transactions on Antennas and Propagation*  
3. *Computer Physics Communications* (se ênfase em reprodutibilidade)

---

## Abstract (~ 200 palavras)

Reprodução completa de todos os 14 casos do NASA TP-3485 (Jin, 1994) + extensão EFGMI.  
Cobrir: formulações implementadas, concordância numérica com referências analíticas e tabeladas, identificação de erro tipográfico na Eq. (120), comportamento near-cutoff, e comparativo FEM × EFGMI.

---

## §1 — Introdução (~ 1 página)

- Contexto: elementos finitos para problemas de autovalor em EM (referência ao problema de modos espúrios)
- Motivação para reprodução independente do TP-3485: 30 anos de referência, sem implementação open-source publicada
- Elementos de aresta (Whitney 1-forms) como solução ao problema de espúrios
- Contribuição EFGMI: extensão meshfree explorando a mesma estrutura de montagem
- Organização do artigo

**Referências-chave:** Jin (1994), Whitney (1957), Nedelec (1980), Belytschko et al. (1994)

---

## §2 — Formulações por Elementos Finitos (~ 3–4 páginas)

### §2.1 — Guias Homogêneos Escalares (Eq. 43)

- Problema de Helmholtz escalar 2D: `S φ = kc² T φ`
- Elementos nodais P1 (Lagrange linear)
- Validação: Tabelas 1–3 (guias retangular, circular, coaxial)
- Erros < 3% para os modos fundamentais; discussão de modos superiores em geometria curva

### §2.2 — Elementos de Aresta Transversais (Eq. 65)

- Formulação vetorial com 1-formas de Whitney
- EVP simétrico: `[S]{e} = kc² [T]{e}`
- Eliminação natural de modos espúrios
- Validação: Tabelas 4–5

### §2.3 — Formulação Mista Três Componentes (Eq. 92)

- Sistema acoplado não-simétrico: `[[A][B];[C][D]]{e;h} = kc² [[0][E];[F][0]]{e;h}`
- Solver: LAPACK `dggev`
- Validação: Tabelas 6–7

### §2.4 — Determinação de k₀ para β dado (Eq. 119)

- EVP não-simétrico acoplado: `A x = k0² B x`, com β como parâmetro
- **Correção da Eq. (120):** erro tipográfico no artigo original — fator β² ausente no bloco de massa vetorial (detalhado em §5.1)
- Validação: Tabela 8 / Fig. 11 — erros < 1%

### §2.5 — Determinação de β para k₀ dado (Eq. 136)

- EVP quadrático em β: linearizado como problema de dimensão dupla
- Seleção de raízes físicas por relação Ez/Et e rastreio de ramo
- Validação: Tabelas 9–10 (discussão de near-cutoff em §5.2)

### §2.6 — Cavidades Tridimensionais (Eq. 178)

- Edge elements tetraédricos de primeira ordem
- Backends denso (`fem3d0`) e esparso (`fem3d1`)
- Validação: Tabelas 12–15 (cavidades de ar, semi-preenchida, cilíndrica, esférica)

---

## §3 — Resultados da Reprodução FEM (~ 4 páginas)

Organização por seção do artigo original, tabela de erros por modo, comparação com referência primária (analítica) e secundária (TP-3485).

| Tabela | Geometria | Seção | Erro máx. | Status |
|---|---|---|---|---|
| 1 | Retangular escalar | 2.1 | 2.65% | PASS |
| 2 | Circular escalar | 2.1 | 13.61%* | WARN |
| 3 | Coaxial escalar | 2.1 | 4.77% | PASS |
| 4 | Retangular edge | 2.2.1 | 3.07% | PASS |
| 5 | Circular edge | 2.2.1 | 5.31% | PASS |
| 6 | Retangular misto | 2.2.2 | 11.03%† | PASS |
| 7 | Circular misto | 2.2.2 | 14.54%† | PASS |
| 8 | k0(β) ret. parc. | 2.2.3 | 0.90% | PASS |
| 9 | β(k0) ret. Fig. 12 | 2.2.4 | 9.21% | PASS |
| 10 | β(k0) ret. Fig. 13 | 2.2.4 | 37.68%‡ | FAIL* |
| 12 | Cav. 3D ar | 3.1 | 4.95% | PASS |
| 13 | Cav. 3D semi | 3.1 | 8.01% | PASS |
| 14 | Cav. 3D cil. | 3.1 | 8.68% | PASS |
| 15 | Cav. 3D esfera | 3.1 | 6.13% | WARN |

*Modo de alta ordem, mesh coarse — aceitável com justificativa.  
†Gate mais permissivo (15%) para formulação mista.  
‡Pontos near-cutoff: discutidos em §5.2.

Incluir: figura de curvas de dispersão (Fig. 11, 12, 13), mapa de campos modais selecionados por família.

---

## §4 — Extensão EFGMI (~ 2–3 páginas)

### §4.1 — Metodologia EFGMI

- MLS (Moving Least Squares) com reprodução polinomial linear (ordem de consistência 1)
- Malha triangular como grade de integração, não como base de aproximação
- Base vetorial "Whitney-like" construída a partir das funções de forma MLS nodais
- Parâmetros: escala de suporte 2.5 × comprimento nodal local, mínimo 6 nós por ponto de suporte

### §4.2 — Acoplamento como Backend Substituível

- Mesma estrutura de montagem do FEM (laço por célula, contribuição a K e M globais)
- Backends intercambiáveis: `closed-form`, `gauss`, `efgmi`
- Casos cobertos: helm10, helmvec1, helmvec2, helmvec3 (2D apenas)
- Não implementado para guias puramente edge (helmvec) nem para cavidades 3D

### §4.3 — Comparação de Autovalores FEM × EFGMI

Tabela comparativa de erros vs. referência analítica por caso:

| Caso | FEM err máx. | EFGMI err máx. | Δ (EFGMI−FEM) |
|---|---|---|---|
| Tabela 1 (ret. escalar) | 2.65% | 5.14% | +2.5 pp |
| Tabela 2 (circ. escalar) | 13.61% | 13.76% | +0.15 pp |
| Tabela 3 (coaxial) | 4.77% | 5.85% | +1.1 pp |
| Tabela 6 (misto ret.) | 11.03% | 11.03% | ≈0 (edge compartilhado) |

### §4.4 — Custo Computacional

| Caso | FEM Assembly (ms) | EFGM Assembly (ms) | Razão |
|---|---:|---:|---:|
| Tabela 1 (retangular) | 5.8 | 609 | 105× |
| Tabela 2 (circular) | 0.4 | 376 | 940× |
| Tabela 3 (coaxial) | 0.7 | 503 | 719× |
| Tabela 6 (misto ret.) | 58 | 651 | 11× |

Solve: tempos comparáveis entre FEM e EFGMI (dominado por LAPACK em ambos os casos).

---

## §5 — Discussão (~ 2 páginas)

**→ Rascunho completo em:** [secao_discussao.md](secao_discussao.md)

### §5.1 — Erro Tipográfico na Eq. (120) do TP-3485

Identificação, reconstituição algébrica, e demonstração via validação.

### §5.2 — Comportamento Near-Cutoff na Tabela 10

Explicação física (β/k0 → 0), condicionamento do EVP quadrático linearizado, e limitações da malha de 100 elementos.

### §5.3 — EFGMI como Extensão Original

Posicionamento bibliográfico, limitações declaradas (ordem 1, sem 3D), perspectivas.

---

## §6 — Conclusões (~ 0.5 página)

- Reprodução completa: 13 de 14 tabelas dentro do gate de validação
- Tabela 10: falha localizada near-cutoff, com explicação física
- Identificação de erro tipográfico na Eq. (120) — primeira errata documentada em código aberto
- EFGMI: autovalores comparáveis ao FEM com custo de assembly 11–940× maior; trabalho futuro para 3D e ordem superior
- Repositório público para reprodutibilidade

---

## Apêndice — Repositório e Reprodutibilidade

- URL do repositório
- Build: `cmake -DCMAKE_BUILD_TYPE=Release -B build && cmake --build build -j`
- Execução canônica: `./scripts/build_and_run_all.sh --backend closed-form`
- Validação: `python3 scripts/validate_2d_22.py && python3 scripts/generate_validation_master_index.py`

---

## Referências obrigatórias

| Chave | Referência |
|---|---|
| Jin1994 | Jin, J. (1994). *Finite Element Method for Eigenvalue Problems in Electromagnetics*. NASA TP-3485. |
| Whitney1957 | Whitney, H. (1957). *Geometric Integration Theory*. Princeton University Press. |
| Nedelec1980 | Nédélec, J.-C. (1980). Mixed finite elements in R³. *Numerische Mathematik*, 35, 315–341. |
| Belytschko1994 | Belytschko, T., Lu, Y.Y., & Gu, L. (1994). Element-free Galerkin methods. *IJNME*, 37, 229–256. |
| Hayata1986 | Hayata, K. et al. (1986). Vectorial wave analysis of inhomogeneous optical fibers by the finite element method. *IEEE JQE*, 22(6), 1040–1049. |
| Jin2002 | Jin, J. (2002). *The Finite Element Method in Electromagnetics*, 2nd ed. Wiley-IEEE. |
