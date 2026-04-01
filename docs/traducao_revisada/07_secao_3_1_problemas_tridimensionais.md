# Secao 3.1 - Problemas Tridimensionais em Cavidades

## Escopo

Esta secao estende os elementos de aresta para tetraedros 3D e formula o problema de autovalor para cavidades fechadas. No repositorio, ela corresponde aos programas `FEM3D0` e `FEM3D1`.

## Ideia central

O artigo reaproveita a filosofia da Secao 2.2:

- usar elementos de aresta para evitar modos espurios
- trabalhar com continuidade tangencial
- montar o problema vetorial em termos de graus de liberdade associados a arestas

## Revisao das equacoes

- Eq. (143): equacao vetorial de onda 3D.
- Eq. (144): forma integral inicial apos multiplicar por `T`.
- Eq. (145): identidade vetorial para `A . (curl B)`.
- Eq. (146): forma integrada por partes ainda com termo de contorno.
- Eq. (147): teorema da divergencia em 3D.
- Eq. (148): identidade de produto vetorial no contorno.
- Eq. (149): forma fraca com termo de superficie explicito.
- Eq. (150): forma fraca final em cavidade PEC.
- Eq. (151): expansao do campo em seis funcoes de aresta tetraedricas.
- Eq. (152): definicao da base de Whitney 3D.
- Eq. (153) a Eq. (156): coordenadas simplex `alpha_t1` a `alpha_t4`.
- Eq. (157): volume do tetraedro.
- Eq. (158) a Eq. (161): subvolumes usados nas coordenadas simplex.
- Eq. (162): forma afim de `alpha_ti`.
- Eq. (163): expressao explicita de `W_m` no tetraedro.
- Eq. (164): definicao de `A_xm`.
- Eq. (165): definicao de `B_xm`.
- Eq. (166): definicao de `C_xm`.
- Eq. (167): definicao de `A_ym`.
- Eq. (168): definicao de `B_ym`.
- Eq. (169): definicao de `C_ym`.
- Eq. (170): definicao de `A_zm`.
- Eq. (171): definicao de `B_zm`.
- Eq. (172): definicao de `C_zm`.
- Eq. (173): propriedade fundamental da base de aresta no contorno das arestas.
- Eq. (174): equacao elemental 3D apos substituir a expansao de `E`.
- Eq. (175): forma matricial elemental.
- Eq. (176): matriz elemental de rigidez vetorial 3D.
- Eq. (177): matriz elemental de massa vetorial 3D.
- Eq. (178): sistema global `S e = k0^2 T e`.
- Eq. (179): `curl W_m` em forma explicita.
- Eq. (180): produto `curl W_m . curl W_n`.
- Eq. (181): forma fechada de `S_el`.
- Eq. (182): forma fechada de `T_el`.

## Resultado teorico principal

O sistema global da secao e:

```text
S e = k0^2 T e
```

com base de aresta tetraedrica:

```text
W_m = L_m (alpha_i grad alpha_j - alpha_j grad alpha_i)
```

e blocos locais fechados:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(324 V^3) * K_mn
T_e(m,n) = eps_r * (L_m L_n)/(1296 V^3) * sum_{k=1}^{10} I_k
```

## Casos numericos da secao

- Figura 15 / Tabela 12:
  - cavidade retangular preenchida com ar
- Figura 16 / Tabela 13:
  - cavidade retangular meio preenchida
- Figura 17 / Tabela 14:
  - cavidade cilindrica circular preenchida com ar
- Tabela 15:
  - cavidade esferica

## Leitura para o codigo

- Eq. (178) e o alvo da montagem global em `src/edge3d`.
- Eq. (181) e Eq. (182) sao a base do backend fechado usado nos modulos 3D.
- A estrutura de coeficientes `A_xm`, `B_xm`, `C_xm`, ... , `C_zm` explica diretamente a nomenclatura dos helpers 3D do repositorio.

## Ajustes em relacao a documentacao legada

- A traducao antiga de 3D parava na Eq. (173) e nao fechava a derivacao ate as Eq. (181) e (182).
- A Eq. (152) estava quebrada por OCR como `W*m`.
- As Eq. (157) a (161) estavam com determinantes em uma unica linha, sem quebras validas de LaTeX.
- As figuras 14 a 17 nao estavam integradas como uma trilha documental coerente; a Figura 14 estava apenas como placeholder de imagem.
