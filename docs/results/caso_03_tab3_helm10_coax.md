# Caso 03 — Tabela 3 — Linha Coaxial Escalar (Sec. 2.1)

## Problema

- **Seção do artigo:** 2.1
- **Formulação:** Eq. (43)
- **Grandeza calculada:** `kc`
- **Geometria:** Linha coaxial homogênea, r₂/r₁ = 4
- **Teoria:** [2.1_Guias de Onda Homogêneos.md](../traducao/2.1_Guias de Onda Homogêneos.md)
- **Rastreabilidade:** [Rastreabilidade_Equacoes_Artigo_Codigo.md](../Rastreabilidade_Equacoes_Artigo_Codigo.md)

## Figura do artigo

![Figura 8](../figs/figura8.png)

## Condições de cálculo

| Parâmetro | Valor |
|---|---|
| Malha | 187 nós, 340 triângulos (nr=10, nt=17) |
| Backend | `closed-form` |

```bash
./build/helm10_coax --nr 10 --nt 17 --nmodos 10 --backend closed-form
```

## Resultados — FEM

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_p1 | TE | 0.4111 | 0.4199 | 2.126 | 0.904648 |
| TE_m1_p1 | TE | 0.4111 | 0.4199 | 2.126 | 0.904648 |
| TE_m2_p1 | TE | 0.7523 | 0.7715 | 2.543 | 0.720001 |
| TE_m2_p1 | TE | 0.7523 | 0.7715 | 2.543 | 0.720001 |
| TE_m3_p1 | TE | 1.0484 | 1.0984 | 4.771 | 0.881867 |
| TE_m3_p1 | TE | 1.0484 | 1.0984 | 4.771 | 0.881867 |
| TE_m0_p1 | TE | 1.1119 | 1.1356 | 2.133 | 0.999994 |
| TE_m1_p2 | TE | 1.2511 | 1.2994 | 3.857 | 0.940877 |
| TE_m1_p2 | TE | 1.2511 | 1.2994 | 3.857 | 0.940877 |
| TE_m4_p1 | TE | 1.3291 | 1.4404 | 8.376 | 0.968396 |

### Gráficos de resumo — FEM

![helm10_coax_error_by_mode](img/fem/helm10/coax/helm10_coax_error_by_mode.png)
![helm10_coax_rho_by_mode](img/fem/helm10/coax/helm10_coax_rho_by_mode.png)


### Campos modais — FEM

![helm10_coax_TE_m0_p1_rank07_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m0_p1_rank07_psi.png)
![helm10_coax_TE_m1_p1_rank01_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m1_p1_rank01_psi.png)
![helm10_coax_TE_m1_p1_rank02_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m1_p1_rank02_psi.png)

![helm10_coax_TE_m1_p2_rank08_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m1_p2_rank08_psi.png)
![helm10_coax_TE_m1_p2_rank09_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m1_p2_rank09_psi.png)
![helm10_coax_TE_m2_p1_rank03_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m2_p1_rank03_psi.png)

![helm10_coax_TE_m2_p1_rank04_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m2_p1_rank04_psi.png)
![helm10_coax_TE_m3_p1_rank05_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m3_p1_rank05_psi.png)
![helm10_coax_TE_m3_p1_rank06_psi](img/fem/helm10/coax/isopotential/helm10_coax_TE_m3_p1_rank06_psi.png)


## Tempo de execução

| Backend | Assembly (ms) | Solve (ms) | Post (ms) | Total (ms) |
|---|---:|---:|---:|---:|
| FEM | 0.6 | 51.7 | 696.6 | 1202.7 |

| EFGMI | 503.2 | 94.6 | 2154.5 | 3834.2 |


## Resultados — EFGMI

| Modo | Familia | `kc` analitico | `kc` FEM | Erro (%) | ρ |
|---|---|---:|---:|---:|---:|
| TE_m1_p1 | TE | 0.4111 | 0.4193 | 1.992 | 0.986374 |
| TE_m1_p1 | TE | 0.4111 | 0.4193 | 1.992 | 0.986374 |
| TE_m2_p1 | TE | 0.7523 | 0.7790 | 3.542 | 0.961584 |
| TE_m2_p1 | TE | 0.7523 | 0.7790 | 3.542 | 0.961584 |
| TE_m3_p1 | TE | 1.0484 | 1.1098 | 5.854 | 0.862542 |
| TE_m3_p1 | TE | 1.0484 | 1.1098 | 5.854 | 0.862542 |
| TE_m0_p1 | TE | 1.1119 | 1.1536 | 3.748 | 0.999737 |
| TE_m1_p2 | TE | 1.2511 | 1.3103 | 4.727 | 0.711532 |
| TE_m1_p2 | TE | 1.2511 | 1.3103 | 4.727 | 0.711532 |
| TE_m4_p1 | TE | 1.3291 | 1.4430 | 8.573 | 0.990783 |

### Gráficos de resumo — EFGMI

![helm10_coax_error_by_mode](img/efgm/helm10/coax/helm10_coax_error_by_mode.png)
![helm10_coax_rho_by_mode](img/efgm/helm10/coax/helm10_coax_rho_by_mode.png)


### Campos modais — EFGMI

![helm10_coax_TE_m0_p1_rank07_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m0_p1_rank07_psi.png)
![helm10_coax_TE_m1_p1_rank01_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m1_p1_rank01_psi.png)
![helm10_coax_TE_m1_p1_rank02_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m1_p1_rank02_psi.png)

![helm10_coax_TE_m1_p2_rank08_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m1_p2_rank08_psi.png)
![helm10_coax_TE_m1_p2_rank09_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m1_p2_rank09_psi.png)
![helm10_coax_TE_m2_p1_rank03_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m2_p1_rank03_psi.png)

![helm10_coax_TE_m2_p1_rank04_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m2_p1_rank04_psi.png)
![helm10_coax_TE_m3_p1_rank05_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m3_p1_rank05_psi.png)
![helm10_coax_TE_m3_p1_rank06_psi](img/efgm/helm10/coax/isopotential/helm10_coax_TE_m3_p1_rank06_psi.png)



---

[← Caso 02](caso_02_tab2_helm10_circle.md)  [Caso 04 →](caso_04_tab4_helmvec_rect.md)  [↑ Índice de Resultados](README.md)
