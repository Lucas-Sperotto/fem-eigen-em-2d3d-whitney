# Caso 07 — Tabela 7 — Guia Circular Misto 3 Comp. (Sec. 2.2.2)

## Problema

- **Seção do artigo:** 2.2.2
- **Formulação:** Eq. (92)
- **Grandeza calculada:** `kc` sistema misto
- **Geometria:** Guia circular homogêneo, raio unitário, formulações E e H
- **Teoria:** [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](../traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 6](../figs/figura6.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 121 nós, 225 triângulos (nr=8, nt=15) |
| Backend | `closed-form` |

```bash
./build/mixed_circle --nr 8 --nt 15 --backend closed-form
```

## Resultados — FEM

| Form. | Bloco dom. | Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---|---|---:|---:|---:|---:|
| E | edge | TE_m1_p1 | TE | 1.8412 | 1.8726 | 1.707 | 0.967477 |
| E | edge | TE_m1_p1 | TE | 1.8412 | 1.8726 | 1.707 | 0.967477 |
| E | edge | TE_m2_p1 | TE | 3.0542 | 3.1486 | 3.090 | 0.992281 |
| E | edge | TE_m2_p1 | TE | 3.0542 | 3.1486 | 3.090 | 0.992281 |
| E | edge | TE_m0_p1 | TE | 3.8317 | 3.7704 | 1.600 | 0.999622 |
| E | edge | TE_m3_p1 | TE | 4.2012 | 4.3931 | 4.569 | 0.991269 |
| E | edge | TE_m3_p1 | TE | 4.2012 | 4.3931 | 4.569 | 0.991269 |
| E | edge | TE_m1_p2 | TE | 5.3314 | 5.2175 | 2.138 | 0.878864 |
| E | edge | TE_m1_p2 | TE | 5.3314 | 5.2175 | 2.138 | 0.878864 |
| E | edge | TE_m4_p1 | TE | 5.3176 | 5.6395 | 6.054 | 0.955357 |
| E | edge | TE_m4_p1 | TE | 5.3176 | 5.6395 | 6.054 | 0.955357 |
| E | edge | TE_m2_p2 | TE | 6.7061 | 6.5733 | 1.981 | 0.713153 |

### Gráficos de resumo — FEM

![helmvec1_circle_error_by_mode](img/fem/helmvec1/circle/helmvec1_circle_error_by_mode.png)
![helmvec1_circle_rho_by_mode](img/fem/helmvec1/circle/helmvec1_circle_rho_by_mode.png)


### Campos modais — FEM

![helmvec1_circle_block_energy_by_mode](img/fem/helmvec1/circle/helmvec1_circle_block_energy_by_mode.png)
![helmvec1_circle_dominant_ratio_by_mode](img/fem/helmvec1/circle/helmvec1_circle_dominant_ratio_by_mode.png)
![helmvec1_circle_kc_by_mode](img/fem/helmvec1/circle/helmvec1_circle_kc_by_mode.png)

![helmvec1_circle_E_Et_TE_m0_p1_rank05_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m0_p1_rank05_magnitude.png)
![helmvec1_circle_E_Et_TE_m0_p2_rank14_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m0_p2_rank14_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p1_rank01_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p1_rank01_magnitude.png)

![helmvec1_circle_E_Et_TE_m1_p1_rank02_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p1_rank02_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p2_rank08_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p2_rank08_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p2_rank09_magnitude](img/fem/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p2_rank09_magnitude.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 10.562933 | 264.349257 |  | 2561.131807 |

| EFGMI | 326.437452 | 204.470579 |  | 7821.369633 |


## Resultados — EFGMI

| Form. | Bloco dom. | Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---|---|---:|---:|---:|---:|
| E | edge | TE_m1_p1 | TE | 1.8412 | 1.8726 | 1.707 | 0.967477 |
| E | edge | TE_m1_p1 | TE | 1.8412 | 1.8726 | 1.707 | 0.967477 |
| E | edge | TE_m2_p1 | TE | 3.0542 | 3.1486 | 3.090 | 0.992281 |
| E | edge | TE_m2_p1 | TE | 3.0542 | 3.1486 | 3.090 | 0.992281 |
| E | edge | TE_m0_p1 | TE | 3.8317 | 3.7704 | 1.600 | 0.999622 |
| E | edge | TE_m3_p1 | TE | 4.2012 | 4.3931 | 4.569 | 0.991269 |
| E | edge | TE_m3_p1 | TE | 4.2012 | 4.3931 | 4.569 | 0.991269 |
| E | edge | TE_m1_p2 | TE | 5.3314 | 5.2175 | 2.138 | 0.878864 |
| E | edge | TE_m1_p2 | TE | 5.3314 | 5.2175 | 2.138 | 0.878864 |
| E | edge | TE_m4_p1 | TE | 5.3176 | 5.6395 | 6.054 | 0.955357 |
| E | edge | TE_m4_p1 | TE | 5.3176 | 5.6395 | 6.054 | 0.955357 |
| E | edge | TE_m2_p2 | TE | 6.7061 | 6.5733 | 1.981 | 0.713153 |

### Gráficos de resumo — EFGMI

![helmvec1_circle_error_by_mode](img/efgm/helmvec1/circle/helmvec1_circle_error_by_mode.png)
![helmvec1_circle_rho_by_mode](img/efgm/helmvec1/circle/helmvec1_circle_rho_by_mode.png)


### Campos modais — EFGMI

![helmvec1_circle_block_energy_by_mode](img/efgm/helmvec1/circle/helmvec1_circle_block_energy_by_mode.png)
![helmvec1_circle_dominant_ratio_by_mode](img/efgm/helmvec1/circle/helmvec1_circle_dominant_ratio_by_mode.png)
![helmvec1_circle_kc_by_mode](img/efgm/helmvec1/circle/helmvec1_circle_kc_by_mode.png)

![helmvec1_circle_E_Et_TE_m0_p1_rank05_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m0_p1_rank05_magnitude.png)
![helmvec1_circle_E_Et_TE_m0_p2_rank14_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m0_p2_rank14_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p1_rank01_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p1_rank01_magnitude.png)

![helmvec1_circle_E_Et_TE_m1_p1_rank02_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p1_rank02_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p2_rank08_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p2_rank08_magnitude.png)
![helmvec1_circle_E_Et_TE_m1_p2_rank09_magnitude](img/efgm/helmvec1/circle/magnitude/helmvec1_circle_E_Et_TE_m1_p2_rank09_magnitude.png)



---

[← Caso 06](caso_06_tab6_helmvec1_rect.md)  [Caso 08 →](caso_08_fig11_tab8_helmvec2.md)  [↑ Índice de Resultados](README.md)
