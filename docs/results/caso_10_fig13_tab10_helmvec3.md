# Caso 10 — Figura 13 / Tabela 10 — Dispersão, `β` dado `k0`, Exemplo 2 (Sec. 2.2.4)

## Problema

- **Seção do artigo:** 2.2.4
- **Formulação:** Eq. (136)
- **Grandeza calculada:** `β/k0` — múltiplos ramos por d/a
- **Geometria:** Guia retangular, b/a=0.45, εr=2.45, d/a variável
- **Teoria:** [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](../traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 13](../figs/figura13.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 66 nós, 100 triângulos (nx=10, ny=5) |
| Backend | `closed-form` |

```bash
./build/helmvec3_fig13_rect --d-over-a-preview 0.20 --nx 10 --ny 5 --backend closed-form
```

## Resultados — FEM

| `d/a` | `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo/HELMVEC3 | Erro analítico (%) | Erro artigo/HELMVEC3 (%) |
|---|---|---:|---:|---:|---:|---:|
| 0.000 | 0.500 | 0.0295 | 0.0300 | 0.0400 | 1.808 | 26.356 |
| 0.000 | 0.600 | 0.5619 | 0.5200 | 0.5600 | 8.058 | 0.339 |
| 0.000 | 0.700 | 0.7012 | 0.7000 | 0.7100 | 0.166 | 1.244 |
| 0.000 | 0.800 | 0.7843 | 0.7900 | 0.7800 | 0.724 | 0.549 |
| 0.000 | 0.900 | 0.8342 | 0.8300 | 0.8300 | 0.505 | 0.505 |
| 0.000 | 1.000 | 0.8681 | 0.8800 | 0.8700 | 1.348 | 0.214 |
| 0.167 | 0.500 | 0.1380 | 0.2100 | 0.1800 | 34.290 | 23.338 |
| 0.167 | 0.600 | 0.5219 | 0.6000 | 0.5900 | 13.015 | 11.540 |
| 0.167 | 0.700 | 0.7070 | 0.7200 | 0.7400 | 1.810 | 4.464 |
| 0.167 | 0.800 | 0.7659 | 0.8200 | 0.8100 | 6.601 | 5.448 |
| 0.167 | 0.900 | 0.8583 | 0.8800 | 0.8700 | 2.467 | 1.346 |
| 0.167 | 1.000 | 0.9320 | 0.9100 | 0.9000 | 2.414 | 3.552 |
| 0.286 | 0.500 | 0.7022 | 0.5100 | 0.4400 | 37.677 | 59.580 |
| 0.286 | 0.600 | 0.7070 | 0.7800 | 0.7400 | 9.363 | 4.464 |
| 0.286 | 0.700 | 0.9083 | 0.9000 | 0.8800 | 0.926 | 3.219 |
| 0.286 | 0.800 | 1.0163 | 0.9900 | 1.0500 | 2.661 | 3.205 |
| 0.286 | 0.900 | 0.9840 | 1.0300 | 1.0300 | 4.470 | 4.470 |
| 0.286 | 1.000 | 1.1038 | 1.1000 | 1.0900 | 0.345 | 1.266 |
| 0.375 | 0.500 | 0.7027 | 0.6800 | 0.6600 | 3.336 | 6.467 |
| 0.375 | 0.600 | 0.8880 | 0.9100 | 0.9000 | 2.420 | 1.335 |
| 0.375 | 0.700 | 1.1021 | 1.0500 | 1.0300 | 4.965 | 7.003 |
| 0.375 | 0.800 | 1.1066 | 1.1300 | 1.1100 | 2.074 | 0.310 |
| 0.375 | 0.900 | 1.1706 | 1.2000 | 1.1800 | 2.450 | 0.796 |
| 0.375 | 1.000 | 1.2368 | 1.2500 | 1.2300 | 1.053 | 0.556 |
| 0.500 | 0.400 | 0.5229 | 0.4000 | 0.4200 | 30.713 | 24.489 |
| 0.500 | 0.500 | 0.9339 | 0.9000 | 0.8900 | 3.765 | 4.931 |
| 0.500 | 0.600 | 1.1011 | 1.1000 | 1.0900 | 0.099 | 1.018 |
| 0.500 | 0.700 | 1.1709 | 1.2000 | 1.1900 | 2.425 | 1.606 |
| 0.500 | 0.800 | 1.2432 | 1.2500 | 1.2400 | 0.541 | 0.261 |
| 0.500 | 0.900 | 1.3013 | 1.3000 | 1.3100 | 0.100 | 0.665 |
| 0.500 | 1.000 | 1.3469 | 1.3500 | 1.3500 | 0.228 | 0.228 |
| 0.600 | 0.400 | 0.7042 | 0.7000 | 0.6700 | 0.600 | 5.104 |
| 0.600 | 0.500 | 1.0214 | 1.0200 | 1.0300 | 0.141 | 0.831 |
| 0.600 | 0.600 | 1.1675 | 1.1800 | 1.1900 | 1.061 | 1.893 |
| 0.600 | 0.700 | 1.2134 | 1.2300 | 1.2700 | 1.350 | 4.457 |
| 0.600 | 0.800 | 1.3295 | 1.3100 | 1.3300 | 1.489 | 0.037 |
| 0.600 | 0.900 | 1.3756 | 1.3800 | 1.3700 | 0.320 | 0.408 |
| 0.600 | 1.000 | 1.4093 | 1.4100 | 1.4000 | 0.048 | 0.666 |
| 0.800 | 0.400 | 0.8725 | 0.9000 | 0.9100 | 3.060 | 4.126 |
| 0.800 | 0.500 | 1.1869 | 1.1800 | 1.1800 | 0.582 | 0.582 |
| 0.800 | 0.600 | 1.3100 | 1.2900 | 1.3000 | 1.548 | 0.767 |
| 0.800 | 0.700 | 1.3808 | 1.3800 | 1.3700 | 0.060 | 0.790 |
| 0.800 | 0.800 | 1.4253 | 1.4100 | 1.4200 | 1.084 | 0.372 |
| 0.800 | 0.900 | 1.4544 | 1.4300 | 1.4400 | 1.705 | 0.999 |
| 0.800 | 1.000 | 1.4201 | 1.4400 | 1.4700 | 1.384 | 3.397 |

### Gráficos de resumo — FEM

![helmvec3_fig13_rect_preview_branch](img/fem/helmvec3/fig13_rect/helmvec3_fig13_rect_preview_branch.png)
![helmvec3_fig13_rect_table10_error_by_branch](img/fem/helmvec3/fig13_rect/helmvec3_fig13_rect_table10_error_by_branch.png)

![helmvec3_fig13_rect_table10_fem_branches](img/fem/helmvec3/fig13_rect/helmvec3_fig13_rect_table10_fem_branches.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Et_magnitude.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Et_magnitude.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Et_magnitude](img/fem/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Et_quiver.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Et_quiver.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Et_quiver](img/fem/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Ez_scalar.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Ez_scalar.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Ez_scalar.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Ez_scalar](img/fem/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Ez_scalar.png)

## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 11.3 | 2474.0 | 3102.5 | 5587.8 |

| EFGMI | 2688.4 | 2518.4 | 3205.1 | 8411.9 |


## Resultados — EFGMI

| `d/a` | `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo/HELMVEC3 | Erro analítico (%) | Erro artigo/HELMVEC3 (%) |
|---|---|---:|---:|---:|---:|---:|
| 0.000 | 0.500 | 0.0301 | 0.0300 | 0.0400 | 0.185 | 24.861 |
| 0.000 | 0.600 | 0.5619 | 0.5200 | 0.5600 | 8.057 | 0.338 |
| 0.000 | 0.700 | 0.7052 | 0.7000 | 0.7100 | 0.737 | 0.682 |
| 0.000 | 0.800 | 0.7914 | 0.7900 | 0.7800 | 0.174 | 1.458 |
| 0.000 | 0.900 | 0.8288 | 0.8300 | 0.8300 | 0.150 | 0.150 |
| 0.000 | 1.000 | 0.8681 | 0.8800 | 0.8700 | 1.347 | 0.213 |
| 0.167 | 0.500 | 0.1352 | 0.2100 | 0.1800 | 35.618 | 24.888 |
| 0.167 | 0.600 | 0.5219 | 0.6000 | 0.5900 | 13.018 | 11.543 |
| 0.167 | 0.700 | 0.7170 | 0.7200 | 0.7400 | 0.419 | 3.110 |
| 0.167 | 0.800 | 0.8209 | 0.8200 | 0.8100 | 0.109 | 1.345 |
| 0.167 | 0.900 | 0.8863 | 0.8800 | 0.8700 | 0.715 | 1.873 |
| 0.167 | 1.000 | 0.9630 | 0.9100 | 0.9000 | 5.830 | 7.006 |
| 0.286 | 0.500 | 0.7188 | 0.5100 | 0.4400 | 40.948 | 63.372 |
| 0.286 | 0.600 | 0.7797 | 0.7800 | 0.7400 | 0.039 | 5.365 |
| 0.286 | 0.700 | 0.9083 | 0.9000 | 0.8800 | 0.927 | 3.221 |
| 0.286 | 0.800 | 0.9447 | 0.9900 | 1.0500 | 4.580 | 10.033 |
| 0.286 | 0.900 | 0.9923 | 1.0300 | 1.0300 | 3.658 | 3.658 |
| 0.286 | 1.000 | 1.0878 | 1.1000 | 1.0900 | 1.113 | 0.206 |
| 0.375 | 0.500 | 0.7215 | 0.6800 | 0.6600 | 6.105 | 9.320 |
| 0.375 | 0.600 | 0.8880 | 0.9100 | 0.9000 | 2.421 | 1.337 |
| 0.375 | 0.700 | 1.0480 | 1.0500 | 1.0300 | 0.186 | 1.752 |
| 0.375 | 0.800 | 1.1447 | 1.1300 | 1.1100 | 1.297 | 3.122 |
| 0.375 | 0.900 | 1.2067 | 1.2000 | 1.1800 | 0.557 | 2.261 |
| 0.375 | 1.000 | 1.2517 | 1.2500 | 1.2300 | 0.140 | 1.768 |
| 0.500 | 0.400 | 0.5228 | 0.4000 | 0.4200 | 30.701 | 24.477 |
| 0.500 | 0.500 | 0.8561 | 0.9000 | 0.8900 | 4.882 | 3.814 |
| 0.500 | 0.600 | 1.0937 | 1.1000 | 1.0900 | 0.574 | 0.338 |
| 0.500 | 0.700 | 1.1956 | 1.2000 | 1.1900 | 0.363 | 0.474 |
| 0.500 | 0.800 | 1.2489 | 1.2500 | 1.2400 | 0.085 | 0.721 |
| 0.500 | 0.900 | 1.3013 | 1.3000 | 1.3100 | 0.099 | 0.665 |
| 0.500 | 1.000 | 1.3469 | 1.3500 | 1.3500 | 0.233 | 0.233 |
| 0.600 | 0.400 | 0.7074 | 0.7000 | 0.6700 | 1.052 | 5.577 |
| 0.600 | 0.500 | 1.0213 | 1.0200 | 1.0300 | 0.128 | 0.844 |
| 0.600 | 0.600 | 1.1828 | 1.1800 | 1.1900 | 0.236 | 0.607 |
| 0.600 | 0.700 | 1.2275 | 1.2300 | 1.2700 | 0.204 | 3.348 |
| 0.600 | 0.800 | 1.3149 | 1.3100 | 1.3300 | 0.375 | 1.135 |
| 0.600 | 0.900 | 1.3756 | 1.3800 | 1.3700 | 0.317 | 0.410 |
| 0.600 | 1.000 | 1.4093 | 1.4100 | 1.4000 | 0.047 | 0.667 |
| 0.800 | 0.400 | 0.8755 | 0.9000 | 0.9100 | 2.725 | 3.794 |
| 0.800 | 0.500 | 1.1811 | 1.1800 | 1.1800 | 0.096 | 0.096 |
| 0.800 | 0.600 | 1.2876 | 1.2900 | 1.3000 | 0.189 | 0.957 |
| 0.800 | 0.700 | 1.3810 | 1.3800 | 1.3700 | 0.074 | 0.805 |
| 0.800 | 0.800 | 1.4254 | 1.4100 | 1.4200 | 1.092 | 0.380 |
| 0.800 | 0.900 | 1.4544 | 1.4300 | 1.4400 | 1.706 | 1.000 |
| 0.800 | 1.000 | 1.4200 | 1.4400 | 1.4700 | 1.391 | 3.404 |

### Gráficos de resumo — EFGMI

![helmvec3_fig13_rect_preview_branch](img/efgm/helmvec3/fig13_rect/helmvec3_fig13_rect_preview_branch.png)
![helmvec3_fig13_rect_table10_error_by_branch](img/efgm/helmvec3/fig13_rect/helmvec3_fig13_rect_table10_error_by_branch.png)

![helmvec3_fig13_rect_table10_fem_branches](img/efgm/helmvec3/fig13_rect/helmvec3_fig13_rect_table10_fem_branches.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Et_magnitude.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Et_magnitude.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Et_magnitude.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Et_magnitude.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Et_magnitude](img/efgm/helmvec3/fig13_rect/magnitude/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Et_magnitude.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Et_quiver.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Et_quiver.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Et_quiver.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Et_quiver.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Et_quiver](img/efgm/helmvec3/fig13_rect/quiver/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Et_quiver.png)
![helmvec3_fig13_rect_preview_da0_2_br0_2_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_2_Ez_scalar.png)

![helmvec3_fig13_rect_preview_da0_2_br0_3_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_3_Ez_scalar.png)
![helmvec3_fig13_rect_preview_da0_2_br0_4_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_4_Ez_scalar.png)

![helmvec3_fig13_rect_preview_da0_2_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_preview_da0_2_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/preview_da_0_2/helmvec3_fig13_rect_preview_da0_2_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0/helmvec3_fig13_rect_table10_da0_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_167_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_167_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_167/helmvec3_fig13_rect_table10_da0_167_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_286_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_286_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_286/helmvec3_fig13_rect_table10_da0_286_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_375_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_375_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_375/helmvec3_fig13_rect_table10_da0_375_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_4_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_4_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_5_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_6_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_7_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_8_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_5_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br0_9_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_5_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_5/helmvec3_fig13_rect_table10_da0_5_br1_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_4_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_4_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_5_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_6_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_7_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_8_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_6_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br0_9_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_6_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_6/helmvec3_fig13_rect_table10_da0_6_br1_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_4_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_4_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_5_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_5_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_6_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_6_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_7_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_7_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br0_8_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_8_Ez_scalar.png)
![helmvec3_fig13_rect_table10_da0_8_br0_9_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br0_9_Ez_scalar.png)

![helmvec3_fig13_rect_table10_da0_8_br1_Ez_scalar](img/efgm/helmvec3/fig13_rect/scalar/table10_da_0_8/helmvec3_fig13_rect_table10_da0_8_br1_Ez_scalar.png)


---

[← Caso 09](caso_09_fig12_tab9_helmvec3.md)  [Caso 11 →](caso_11_tab12_fem3d_air.md)  [↑ Índice de Resultados](README.md)
