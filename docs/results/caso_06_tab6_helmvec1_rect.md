# Caso 06 — Tabela 6 — Guia Retangular Misto 3 Comp. (Sec. 2.2.2)

## Problema

- **Seção do artigo:** 2.2.2
- **Formulação:** Eq. (92) — bloco diagonal [St 0; 0 Sz] e = kc² [Tt 0; 0 Tz] e
- **Grandeza calculada:** `kc` sistema misto
- **Geometria:** Guia retangular homogêneo, formulações E e H
- **Teoria:** [2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md](../traducao/2.2.2_Guias_de_Onda_Nao_Homogeneos_Tres_Componentes.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 4](../figs/figura4.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 231 nós, 400 triângulos (nx=10, ny=20) |
| Backend | `closed-form` |

```bash
./build/mixed_rect --nx 10 --ny 20 --backend closed-form
```

## Resultados — FEM

| Form. | Bloco dom. | Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---|---|---:|---:|---:|---:|
| E | edge | TE_m1_n0 | TE | 3.1416 | 3.1433 | 0.055 | 1.000000 |
| E | edge | TE_m0_n1 | TE | 6.2832 | 6.2332 | 0.795 | 0.999995 |
| E | edge | TE_m2_n0 | TE | 6.2832 | 6.2966 | 0.214 | 0.999968 |
| E | edge | TE_m1_n1 | TE | 7.0248 | 6.9971 | 0.394 | 0.999897 |
| E | edge | TE_m2_n1 | TE | 8.8858 | 8.9113 | 0.287 | 0.999779 |
| E | edge | TE_m3_n0 | TE | 9.4248 | 9.4685 | 0.463 | 0.999653 |
| E | edge | TE_m3_n1 | TE | 11.3272 | 11.4274 | 0.885 | 0.999747 |
| E | edge | TE_m0_n2 | TE | 12.5664 | 12.1809 | 3.068 | 0.999782 |
| E | edge | TE_m1_n2 | TE | 12.9531 | 12.6130 | 2.626 | 0.999657 |
| E | edge | TE_m4_n0 | TE | 12.5664 | 12.6611 | 0.754 | 0.999727 |
| E | edge | TE_m2_n2 | TE | 14.0496 | 13.8360 | 1.521 | 0.999226 |
| E | edge | TE_m4_n1 | TE | 14.0496 | 14.2466 | 1.402 | 0.999141 |

### Gráficos de resumo — FEM

![helmvec1_rect_error_by_mode](img/fem/helmvec1/rect/helmvec1_rect_error_by_mode.png)
![helmvec1_rect_rho_by_mode](img/fem/helmvec1/rect/helmvec1_rect_rho_by_mode.png)


### Campos modais — FEM

![helmvec1_rect_block_energy_by_mode](img/fem/helmvec1/rect/helmvec1_rect_block_energy_by_mode.png)
![helmvec1_rect_dominant_ratio_by_mode](img/fem/helmvec1/rect/helmvec1_rect_dominant_ratio_by_mode.png)
![helmvec1_rect_kc_by_mode](img/fem/helmvec1/rect/helmvec1_rect_kc_by_mode.png)

![helmvec1_rect_E_Et_TE_m0_n1_rank02_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n1_rank02_magnitude.png)
![helmvec1_rect_E_Et_TE_m0_n2_rank08_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n2_rank08_magnitude.png)
![helmvec1_rect_E_Et_TE_m0_n3_rank16_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n3_rank16_magnitude.png)

![helmvec1_rect_E_Et_TE_m1_n0_rank01_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n0_rank01_magnitude.png)
![helmvec1_rect_E_Et_TE_m1_n1_rank04_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n1_rank04_magnitude.png)
![helmvec1_rect_E_Et_TE_m1_n2_rank09_magnitude](img/fem/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n2_rank09_magnitude.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 46.582569 | 5369.461042 |  | 7933.143297 |

| EFGMI | 651.3185100000001 | 765.190259 |  | 16016.452565 |


## Resultados — EFGMI

| Form. | Bloco dom. | Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---|---|---:|---:|---:|---:|
| E | edge | TE_m1_n0 | TE | 3.1416 | 3.1433 | 0.055 | 1.000000 |
| E | edge | TE_m0_n1 | TE | 6.2832 | 6.2332 | 0.795 | 0.999995 |
| E | edge | TE_m2_n0 | TE | 6.2832 | 6.2966 | 0.214 | 0.999968 |
| E | edge | TE_m1_n1 | TE | 7.0248 | 6.9971 | 0.394 | 0.999897 |
| E | edge | TE_m2_n1 | TE | 8.8858 | 8.9113 | 0.287 | 0.999779 |
| E | edge | TE_m3_n0 | TE | 9.4248 | 9.4685 | 0.463 | 0.999653 |
| E | edge | TE_m3_n1 | TE | 11.3272 | 11.4274 | 0.885 | 0.999747 |
| E | edge | TE_m0_n2 | TE | 12.5664 | 12.1809 | 3.068 | 0.999782 |
| E | edge | TE_m1_n2 | TE | 12.9531 | 12.6130 | 2.626 | 0.999657 |
| E | edge | TE_m4_n0 | TE | 12.5664 | 12.6611 | 0.754 | 0.999727 |
| E | edge | TE_m2_n2 | TE | 14.0496 | 13.8360 | 1.521 | 0.999226 |
| E | edge | TE_m4_n1 | TE | 14.0496 | 14.2466 | 1.402 | 0.999141 |

### Gráficos de resumo — EFGMI

![helmvec1_rect_error_by_mode](img/efgm/helmvec1/rect/helmvec1_rect_error_by_mode.png)
![helmvec1_rect_rho_by_mode](img/efgm/helmvec1/rect/helmvec1_rect_rho_by_mode.png)


### Campos modais — EFGMI

![helmvec1_rect_block_energy_by_mode](img/efgm/helmvec1/rect/helmvec1_rect_block_energy_by_mode.png)
![helmvec1_rect_dominant_ratio_by_mode](img/efgm/helmvec1/rect/helmvec1_rect_dominant_ratio_by_mode.png)
![helmvec1_rect_kc_by_mode](img/efgm/helmvec1/rect/helmvec1_rect_kc_by_mode.png)

![helmvec1_rect_E_Et_TE_m0_n1_rank02_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n1_rank02_magnitude.png)
![helmvec1_rect_E_Et_TE_m0_n2_rank08_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n2_rank08_magnitude.png)
![helmvec1_rect_E_Et_TE_m0_n3_rank16_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m0_n3_rank16_magnitude.png)

![helmvec1_rect_E_Et_TE_m1_n0_rank01_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n0_rank01_magnitude.png)
![helmvec1_rect_E_Et_TE_m1_n1_rank04_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n1_rank04_magnitude.png)
![helmvec1_rect_E_Et_TE_m1_n2_rank09_magnitude](img/efgm/helmvec1/rect/magnitude/helmvec1_rect_E_Et_TE_m1_n2_rank09_magnitude.png)



---

[← Caso 05](caso_05_tab5_helmvec_circle.md)  [Caso 07 →](caso_07_tab7_helmvec1_circle.md)  [↑ Índice de Resultados](README.md)
