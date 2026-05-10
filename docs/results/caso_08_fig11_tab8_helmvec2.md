# Caso 08 — Figura 11 / Tabela 8 — Guia Parcialmente Preenchido, `k0` dado `β` (Sec. 2.2.3)

## Problema

- **Seção do artigo:** 2.2.3
- **Formulação:** Eq. (119) — A x = k0² B x, x = [Et; Ez]
- **Grandeza calculada:** `k0L` (número de onda × comprimento)
- **Geometria:** Guia retangular quadrado, dielétrico inferior (εr=1.5), β=10
- **Teoria:** [2.2.3_Determinação do número de onda.md](../traducao/2.2.3_Determinação do número de onda.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 11](../figs/figura11.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 441 nós, 800 triângulos (nx=20, ny=20) |
| Backend | `closed-form` |

```bash
./build/helmvec2_rect --beta 10 --nx 20 --ny 20 --backend closed-form
```

> **Nota:** A Eq. (120) impressa no artigo contém inconsistência de impressão (falta fator β² no bloco de massa vetorial). O código reconstrói A_tt a partir dos blocos elementares validados. Ver [src/helmvec2/README.md](../../src/helmvec2/README.md).

## Resultados — FEM

| Modo | `k0L` FEM | Ref. HELMVEC2 | Ref. Hayata | Erro HELMVEC2 (%) | Erro Hayata (%) |
|---|---:|---:|---:|---:|---:|
| 1 | 8.8097 | 8.8150 | 8.8093 | 0.060 | 0.004 |
| 2 | 9.3914 | 9.4430 | 9.3896 | 0.546 | 0.020 |
| 3 | 10.2813 | 10.3500 | 10.2752 | 0.664 | 0.059 |
| 4 | 11.1126 | 11.1410 | 11.1030 | 0.255 | 0.087 |
| 5 | 11.2613 | 11.2890 | 11.2677 | 0.246 | 0.057 |
| 6 | 11.4224 | 11.4246 | 11.4501 | 0.019 | 0.242 |
| 7 | 12.2547 | 12.1460 | 11.9882 | 0.895 | 2.223 |
| 8 | 12.6430 | 12.5894 | 12.6686 | 0.426 | 0.202 |
| 9 | 12.8433 | 12.8237 | 12.8092 | 0.153 | 0.266 |
| 10 | 12.9101 | 12.9987 | 12.9575 | 0.682 | 0.366 |

### Gráficos de resumo — FEM

![helmvec2_rect_candidates_k0l_by_rank](img/fem/helmvec2/rect/helmvec2_rect_candidates_k0l_by_rank.png)
![helmvec2_rect_error_by_mode](img/fem/helmvec2/rect/helmvec2_rect_error_by_mode.png)

![helmvec2_rect_table8_k0l_by_mode](img/fem/helmvec2/rect/helmvec2_rect_table8_k0l_by_mode.png)

### Campos modais — FEM

![helmvec2_rect_mode01_cand01_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode01_cand01_Et_magnitude.png)
![helmvec2_rect_mode02_cand03_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode02_cand03_Et_magnitude.png)
![helmvec2_rect_mode03_cand05_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode03_cand05_Et_magnitude.png)

![helmvec2_rect_mode04_cand08_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode04_cand08_Et_magnitude.png)
![helmvec2_rect_mode05_cand09_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode05_cand09_Et_magnitude.png)
![helmvec2_rect_mode06_cand11_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode06_cand11_Et_magnitude.png)

![helmvec2_rect_mode07_cand13_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode07_cand13_Et_magnitude.png)
![helmvec2_rect_mode08_cand16_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode08_cand16_Et_magnitude.png)
![helmvec2_rect_mode09_cand17_Et_magnitude](img/fem/helmvec2/rect/magnitude/helmvec2_rect_mode09_cand17_Et_magnitude.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 101.5 | 25653.5 | 1.6 | 27346.5 |

| EFGMI | 924.4 | 24103.6 | 1.3 | 26503.4 |


## Resultados — EFGMI

| Modo | `k0L` FEM | Ref. HELMVEC2 | Ref. Hayata | Erro HELMVEC2 (%) | Erro Hayata (%) |
|---|---:|---:|---:|---:|---:|
| 1 | 8.8100 | 8.8150 | 8.8093 | 0.056 | 0.008 |
| 2 | 9.3890 | 9.4430 | 9.3896 | 0.572 | 0.007 |
| 3 | 10.2809 | 10.3500 | 10.2752 | 0.667 | 0.056 |
| 4 | 11.1167 | 11.1410 | 11.1030 | 0.219 | 0.123 |
| 5 | 11.2556 | 11.2890 | 11.2677 | 0.296 | 0.108 |
| 6 | 11.4216 | 11.4246 | 11.4501 | 0.026 | 0.249 |
| 7 | 12.0445 | 12.1460 | 11.9882 | 0.835 | 0.470 |
| 8 | 12.6372 | 12.5894 | 12.6686 | 0.380 | 0.247 |
| 9 | 12.8596 | 12.8237 | 12.8092 | 0.280 | 0.394 |
| 10 | 12.9069 | 12.9987 | 12.9575 | 0.706 | 0.390 |

### Gráficos de resumo — EFGMI

![helmvec2_rect_candidates_k0l_by_rank](img/efgm/helmvec2/rect/helmvec2_rect_candidates_k0l_by_rank.png)
![helmvec2_rect_error_by_mode](img/efgm/helmvec2/rect/helmvec2_rect_error_by_mode.png)

![helmvec2_rect_table8_k0l_by_mode](img/efgm/helmvec2/rect/helmvec2_rect_table8_k0l_by_mode.png)

### Campos modais — EFGMI

![helmvec2_rect_mode01_cand362_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode01_cand362_Et_magnitude.png)
![helmvec2_rect_mode02_cand364_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode02_cand364_Et_magnitude.png)
![helmvec2_rect_mode03_cand366_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode03_cand366_Et_magnitude.png)

![helmvec2_rect_mode04_cand369_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode04_cand369_Et_magnitude.png)
![helmvec2_rect_mode05_cand370_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode05_cand370_Et_magnitude.png)
![helmvec2_rect_mode06_cand372_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode06_cand372_Et_magnitude.png)

![helmvec2_rect_mode07_cand373_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode07_cand373_Et_magnitude.png)
![helmvec2_rect_mode08_cand377_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode08_cand377_Et_magnitude.png)
![helmvec2_rect_mode09_cand378_Et_magnitude](img/efgm/helmvec2/rect/magnitude/helmvec2_rect_mode09_cand378_Et_magnitude.png)



---

[← Caso 07](caso_07_tab7_helmvec1_circle.md)  [Caso 09 →](caso_09_fig12_tab9_helmvec3.md)  [↑ Índice de Resultados](README.md)
