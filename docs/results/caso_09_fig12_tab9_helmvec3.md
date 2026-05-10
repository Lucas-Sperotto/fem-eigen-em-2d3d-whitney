# Caso 09 — Figura 12 / Tabela 9 — Dispersão, `β` dado `k0`, Exemplo 1 (Sec. 2.2.4)

## Problema

- **Seção do artigo:** 2.2.4
- **Formulação:** Eq. (136) — P x = β² Q x
- **Grandeza calculada:** `β/k0` (relação de dispersão)
- **Geometria:** Guia retangular, b/a=0.45, d/a=0.5, εr=2.45
- **Teoria:** [2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md](../traducao/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 12](../figs/figura12.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 66 nós, 100 triângulos (nx=10, ny=5) |
| Backend | `closed-form` |

```bash
./build/helmvec3_fig12_rect --nx 10 --ny 5 --backend closed-form
```

## Resultados — FEM

| `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo | Erro analitico (%) | Erro artigo (%) |
|---|---:|---:|---:|---:|---:|
| 0.200 | 0.4358 | 0.4800 | 0.4700 | 9.208 | 7.277 |
| 0.300 | 1.0468 | 1.0000 | 1.0100 | 4.683 | 3.646 |
| 0.400 | 1.1065 | 1.1800 | 1.1700 | 6.229 | 5.427 |
| 0.500 | 1.3329 | 1.2600 | 1.2800 | 5.784 | 4.131 |
| 0.600 | 1.3585 | 1.3000 | 1.3500 | 4.502 | 0.632 |

### Gráficos de resumo — FEM

![helmvec3_fig12_rect_table9_beta_over_k0](img/fem/helmvec3/fig12_rect/helmvec3_fig12_rect_table9_beta_over_k0.png)
![helmvec3_fig12_rect_table9_error_by_point](img/fem/helmvec3/fig12_rect/helmvec3_fig12_rect_table9_error_by_point.png)


### Campos modais — FEM

![helmvec3_fig12_rect_figure12_br0_2_Et_magnitude](img/fem/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_2_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_3_Et_magnitude](img/fem/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_3_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_4_Et_magnitude](img/fem/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_4_Et_magnitude.png)

![helmvec3_fig12_rect_figure12_br0_5_Et_magnitude](img/fem/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_5_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_6_Et_magnitude](img/fem/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_6_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_2_Et_quiver](img/fem/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_2_Et_quiver.png)

![helmvec3_fig12_rect_figure12_br0_3_Et_quiver](img/fem/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_3_Et_quiver.png)
![helmvec3_fig12_rect_figure12_br0_4_Et_quiver](img/fem/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_4_Et_quiver.png)
![helmvec3_fig12_rect_figure12_br0_5_Et_quiver](img/fem/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_5_Et_quiver.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 4.1 | 254.2 | 310.0 | 568.2 |

| EFGMI | 320.1 | 188.8 | 341.3 | 850.2 |


## Resultados — EFGMI

| `βr/λ0` | `β/k0` FEM | `β/k0` Analítico | `β/k0` Artigo | Erro analitico (%) | Erro artigo (%) |
|---|---:|---:|---:|---:|---:|
| 0.200 | 0.4401 | 0.4800 | 0.4700 | 8.322 | 6.372 |
| 0.300 | 1.0130 | 1.0000 | 1.0100 | 1.299 | 0.296 |
| 0.400 | 1.1821 | 1.1800 | 1.1700 | 0.175 | 1.031 |
| 0.500 | 1.2540 | 1.2600 | 1.2800 | 0.480 | 2.035 |
| 0.600 | 1.2929 | 1.3000 | 1.3500 | 0.545 | 4.228 |

### Gráficos de resumo — EFGMI

![helmvec3_fig12_rect_table9_beta_over_k0](img/efgm/helmvec3/fig12_rect/helmvec3_fig12_rect_table9_beta_over_k0.png)
![helmvec3_fig12_rect_table9_error_by_point](img/efgm/helmvec3/fig12_rect/helmvec3_fig12_rect_table9_error_by_point.png)


### Campos modais — EFGMI

![helmvec3_fig12_rect_figure12_br0_2_Et_magnitude](img/efgm/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_2_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_3_Et_magnitude](img/efgm/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_3_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_4_Et_magnitude](img/efgm/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_4_Et_magnitude.png)

![helmvec3_fig12_rect_figure12_br0_5_Et_magnitude](img/efgm/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_5_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_6_Et_magnitude](img/efgm/helmvec3/fig12_rect/magnitude/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_6_Et_magnitude.png)
![helmvec3_fig12_rect_figure12_br0_2_Et_quiver](img/efgm/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_2_Et_quiver.png)

![helmvec3_fig12_rect_figure12_br0_3_Et_quiver](img/efgm/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_3_Et_quiver.png)
![helmvec3_fig12_rect_figure12_br0_4_Et_quiver](img/efgm/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_4_Et_quiver.png)
![helmvec3_fig12_rect_figure12_br0_5_Et_quiver](img/efgm/helmvec3/fig12_rect/quiver/figure12_da_0_225/helmvec3_fig12_rect_figure12_br0_5_Et_quiver.png)



---

[← Caso 08](caso_08_fig11_tab8_helmvec2.md)  [Caso 10 →](caso_10_fig13_tab10_helmvec3.md)  [↑ Índice de Resultados](README.md)
