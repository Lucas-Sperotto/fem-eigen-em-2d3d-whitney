# Caso 01 — Tabela 1 — Guia Retangular Escalar (Sec. 2.1)

## Problema

- **Seção do artigo:** 2.1
- **Formulação:** Eq. (43) — S φ = kc² T φ
- **Grandeza calculada:** `kc` (cutoff escalar)
- **Geometria:** Guia retangular homogêneo, a/b = 2, PEC
- **Teoria:** [2.1_Guias de Onda Homogêneos.md](../traducao/2.1_Guias de Onda Homogêneos.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 4](../figs/figura4.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 231 nós, 400 triângulos (nx=10, ny=20) |
| Backend | `closed-form` |

```bash
./build/helm10_rect --ar-m 1.0 --nx 10 --ny 20 --nmodos 10 --backend closed-form
```

## Resultados — FEM

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_n0 | TE | 3.1416 | 3.1545 | 0.411 | 1.000000 |
| TE_m0_n1 | TE | 6.2832 | 6.2895 | 0.100 | 0.999997 |
| TE_m2_n0 | TE | 6.2832 | 6.3867 | 1.648 | 0.999971 |
| TE_m1_n1 | TE | 7.0248 | 7.0852 | 0.860 | 0.999950 |
| TE_m2_n1 | TE | 8.8858 | 9.1197 | 2.632 | 0.999630 |
| TE_m3_n0 | TE | 9.4248 | 9.7758 | 3.724 | 0.999577 |
| TE_m3_n1 | TE | 11.3272 | 11.9019 | 5.074 | 0.999243 |
| TE_m0_n2 | TE | 12.5664 | 12.6169 | 0.402 | 0.999844 |
| TE_m1_n2 | TE | 12.9531 | 13.1098 | 1.210 | 0.999749 |
| TE_m4_n0 | TE | 12.5664 | 13.3978 | 6.616 | 0.999529 |

### Gráficos de resumo — FEM

![helm10_rect_error_by_mode](img/fem/helm10/rect/helm10_rect_error_by_mode.png)
![helm10_rect_rho_by_mode](img/fem/helm10/rect/helm10_rect_rho_by_mode.png)


### Campos modais — FEM

![helm10_rect_TE_m0_n1_rank02_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m0_n1_rank02_psi.png)
![helm10_rect_TE_m0_n2_rank08_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m0_n2_rank08_psi.png)
![helm10_rect_TE_m1_n0_rank01_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m1_n0_rank01_psi.png)

![helm10_rect_TE_m1_n1_rank04_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m1_n1_rank04_psi.png)
![helm10_rect_TE_m1_n2_rank09_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m1_n2_rank09_psi.png)
![helm10_rect_TE_m2_n0_rank03_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m2_n0_rank03_psi.png)

![helm10_rect_TE_m2_n1_rank05_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m2_n1_rank05_psi.png)
![helm10_rect_TE_m3_n0_rank06_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m3_n0_rank06_psi.png)
![helm10_rect_TE_m3_n1_rank07_psi](img/fem/helm10/rect/isopotential/helm10_rect_TE_m3_n1_rank07_psi.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 1.0 | 54.1 | 138.6 | 421.9 |

| EFGMI | 609.5 | 45.7 | 1427.8 | 3225.7 |


## Resultados — EFGMI

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_n0 | TE | 3.1416 | 3.1545 | 0.409 | 0.999998 |
| TE_m0_n1 | TE | 6.2832 | 6.3541 | 1.128 | 0.999997 |
| TE_m2_n0 | TE | 6.2832 | 6.3857 | 1.632 | 0.999965 |
| TE_m1_n1 | TE | 7.0248 | 7.1233 | 1.402 | 0.999998 |
| TE_m2_n1 | TE | 8.8858 | 9.0990 | 2.399 | 0.999961 |
| TE_m3_n0 | TE | 9.4248 | 9.7683 | 3.645 | 0.999811 |
| TE_m3_n1 | TE | 11.3272 | 11.8024 | 4.195 | 0.999775 |
| TE_m0_n2 | TE | 12.5664 | 13.1796 | 4.880 | 0.999323 |
| TE_m4_n0 | TE | 12.5664 | 13.3687 | 6.385 | 0.999275 |
| TE_m1_n2 | TE | 12.9531 | 13.6182 | 5.135 | 0.999322 |

### Gráficos de resumo — EFGMI

![helm10_rect_error_by_mode](img/efgm/helm10/rect/helm10_rect_error_by_mode.png)
![helm10_rect_rho_by_mode](img/efgm/helm10/rect/helm10_rect_rho_by_mode.png)


### Campos modais — EFGMI

![helm10_rect_TE_m0_n1_rank02_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m0_n1_rank02_psi.png)
![helm10_rect_TE_m0_n2_rank08_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m0_n2_rank08_psi.png)
![helm10_rect_TE_m1_n0_rank01_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m1_n0_rank01_psi.png)

![helm10_rect_TE_m1_n1_rank04_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m1_n1_rank04_psi.png)
![helm10_rect_TE_m1_n2_rank10_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m1_n2_rank10_psi.png)
![helm10_rect_TE_m2_n0_rank03_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m2_n0_rank03_psi.png)

![helm10_rect_TE_m2_n1_rank05_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m2_n1_rank05_psi.png)
![helm10_rect_TE_m3_n0_rank06_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m3_n0_rank06_psi.png)
![helm10_rect_TE_m3_n1_rank07_psi](img/efgm/helm10/rect/isopotential/helm10_rect_TE_m3_n1_rank07_psi.png)



---

[Caso 02 →](caso_02_tab2_helm10_circle.md)  [↑ Índice de Resultados](README.md)
