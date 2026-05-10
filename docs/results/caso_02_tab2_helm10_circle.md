# Caso 02 — Tabela 2 — Guia Circular Escalar (Sec. 2.1)

## Problema

- **Seção do artigo:** 2.1
- **Formulação:** Eq. (43)
- **Grandeza calculada:** `kc`
- **Geometria:** Guia circular homogêneo, raio unitário, PEC
- **Teoria:** [2.1_Guias de Onda Homogêneos.md](../traducao/2.1_Guias de Onda Homogêneos.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 6](../figs/figura6.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 121 nós, 225 triângulos (nr=8, nt=15) |
| Backend | `closed-form` |

```bash
./build/helm10_circle --nr 8 --nt 15 --nmodos 10 --backend closed-form
```

## Resultados — FEM

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_p1 | TE | 1.8412 | 1.8776 | 1.976 | 0.898898 |
| TE_m1_p1 | TE | 1.8412 | 1.8776 | 1.976 | 0.898898 |
| TE_m2_p1 | TE | 3.0542 | 3.1533 | 3.243 | 0.723562 |
| TE_m2_p1 | TE | 3.0542 | 3.1533 | 3.243 | 0.723562 |
| TE_m0_p1 | TE | 3.8317 | 3.9412 | 2.857 | 0.999937 |
| TE_m3_p1 | TE | 4.2012 | 4.4628 | 6.226 | 0.942674 |
| TE_m3_p1 | TE | 4.2012 | 4.4628 | 6.226 | 0.942674 |
| TE_m1_p2 | TE | 5.3314 | 5.6092 | 5.210 | 0.878038 |
| TE_m1_p2 | TE | 5.3314 | 5.6092 | 5.210 | 0.878038 |
| TE_m4_p1 | TE | 5.3176 | 5.8852 | 10.675 | 0.906001 |

### Gráficos de resumo — FEM

![helm10_circle_error_by_mode](img/fem/helm10/circle/helm10_circle_error_by_mode.png)
![helm10_circle_rho_by_mode](img/fem/helm10/circle/helm10_circle_rho_by_mode.png)


### Campos modais — FEM

![helm10_circle_TE_m0_p1_rank05_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m0_p1_rank05_psi.png)
![helm10_circle_TE_m1_p1_rank01_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m1_p1_rank01_psi.png)
![helm10_circle_TE_m1_p1_rank02_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m1_p1_rank02_psi.png)

![helm10_circle_TE_m1_p2_rank08_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m1_p2_rank08_psi.png)
![helm10_circle_TE_m1_p2_rank09_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m1_p2_rank09_psi.png)
![helm10_circle_TE_m2_p1_rank03_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m2_p1_rank03_psi.png)

![helm10_circle_TE_m2_p1_rank04_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m2_p1_rank04_psi.png)
![helm10_circle_TE_m3_p1_rank06_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m3_p1_rank06_psi.png)
![helm10_circle_TE_m3_p1_rank07_psi](img/fem/helm10/circle/isopotential/helm10_circle_TE_m3_p1_rank07_psi.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 0.2 | 43.4 | 218.3 | 479.7 |

| EFGMI | 376.1 | 49.5 | 819.7 | 1725.2 |


## Resultados — EFGMI

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_p1 | TE | 1.8412 | 1.8847 | 2.363 | 0.931936 |
| TE_m1_p1 | TE | 1.8412 | 1.8847 | 2.363 | 0.931936 |
| TE_m2_p1 | TE | 3.0542 | 3.2026 | 4.858 | 0.997492 |
| TE_m2_p1 | TE | 3.0542 | 3.2026 | 4.858 | 0.997492 |
| TE_m0_p1 | TE | 3.8317 | 4.0014 | 4.429 | 0.999768 |
| TE_m3_p1 | TE | 4.2012 | 4.5498 | 8.299 | 0.872227 |
| TE_m3_p1 | TE | 4.2012 | 4.5498 | 8.299 | 0.872227 |
| TE_m1_p2 | TE | 5.3314 | 5.7384 | 7.633 | 0.891185 |
| TE_m1_p2 | TE | 5.3314 | 5.7384 | 7.633 | 0.891185 |
| TE_m4_p1 | TE | 5.3176 | 5.9428 | 11.757 | 0.822889 |

### Gráficos de resumo — EFGMI

![helm10_circle_error_by_mode](img/efgm/helm10/circle/helm10_circle_error_by_mode.png)
![helm10_circle_rho_by_mode](img/efgm/helm10/circle/helm10_circle_rho_by_mode.png)


### Campos modais — EFGMI

![helm10_circle_TE_m0_p1_rank05_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m0_p1_rank05_psi.png)
![helm10_circle_TE_m1_p1_rank01_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m1_p1_rank01_psi.png)
![helm10_circle_TE_m1_p1_rank02_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m1_p1_rank02_psi.png)

![helm10_circle_TE_m1_p2_rank08_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m1_p2_rank08_psi.png)
![helm10_circle_TE_m1_p2_rank09_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m1_p2_rank09_psi.png)
![helm10_circle_TE_m2_p1_rank03_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m2_p1_rank03_psi.png)

![helm10_circle_TE_m2_p1_rank04_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m2_p1_rank04_psi.png)
![helm10_circle_TE_m3_p1_rank06_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m3_p1_rank06_psi.png)
![helm10_circle_TE_m3_p1_rank07_psi](img/efgm/helm10/circle/isopotential/helm10_circle_TE_m3_p1_rank07_psi.png)



---

[← Caso 01](caso_01_tab1_helm10_rect.md)  [Caso 03 →](caso_03_tab3_helm10_coax.md)  [↑ Índice de Resultados](README.md)
