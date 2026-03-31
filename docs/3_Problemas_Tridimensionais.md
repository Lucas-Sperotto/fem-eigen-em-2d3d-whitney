# 📄 3. Problemas Tridimensionais

## 3.1 Autovalores de Cavidades Tridimensionais — Formulação Vetorial

O problema de calcular frequências de ressonância em cavidades tridimensionais tem sido historicamente dificultado pela presença de **modos espúrios**.

Como mencionado na seção 2.2, esse problema foi recentemente superado utilizando **elementos de aresta tangenciais vetoriais** (refs. 15 e 16).

Nesta seção, formula-se o método de elementos finitos de Galerkin para cavidades tridimensionais e apresentam-se resultados para diversas geometrias.

---

## 3.1.1 Formulação

Essa formulação pode ser desenvolvida utilizando o campo elétrico $\mathbf{E}$ ou o campo magnético $\mathbf{H}$. Aqui, será utilizada a formulação em termos do campo elétrico.

Considere a equação vetorial de onda:

### (143)

$$
\nabla \times \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) - k_0^2 \varepsilon_r \mathbf{E} = 0
$$

onde:

$$
\mathbf{E} = E_x \hat{x} + E_y \hat{y} + E_z \hat{z}
$$

---

Introduz-se uma função teste:

$$
\mathbf{T} = T_x \hat{x} + T_y \hat{y} + T_z \hat{z}
$$

Para aplicar o método de Galerkin, multiplica-se a equação (143) por $\mathbf{T}$ e integra-se sobre o volume da cavidade:

### (144)

$$
\iiint_V \left[ \mathbf{T} \cdot \nabla \times \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) - k_0^2 \varepsilon_r \mathbf{T} \cdot \mathbf{E} \right] dv = 0
$$

---

Utilizando a identidade vetorial:

### (145)

$$
\mathbf{A} \cdot (\nabla \times \mathbf{B}) = (\nabla \times \mathbf{A}) \cdot \mathbf{B} - \nabla \cdot (\mathbf{A} \times \mathbf{B}) 
$$

a equação (144) pode ser reescrita como:

### (146)

$$
\iiint_V (\nabla \times \mathbf{T}) \cdot \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) dv = k_0^2 \varepsilon_r \iiint_V \mathbf{T} \cdot \mathbf{E} ds + \iiint_V \nabla \cdot \left[ \mathbf{T} \times \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) \right] dv
$$

---

Aplicando o **teorema da divergência**:

### (147)

$$
\iiint_V \nabla \cdot \mathbf{A} dv = \iint_S \mathbf{A} \cdot \hat{n} ds
$$

e a identidade:

### (148)

$$
(\mathbf{A} \times \mathbf{B}) \cdot \hat{n} = \mathbf{A} \cdot (\hat{n} \times \mathbf{B})
$$

onde:

* $V$: volume da cavidade
* $S$: superfície externa
* $\hat{n}$: vetor normal unitário externo

---

A equação (146) torna-se:

### (149)

$$
\iiint_V (\nabla \times \mathbf{T}) \cdot \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) dv = k_0^2 \varepsilon_r \iiint_V \mathbf{T} \cdot \mathbf{E} dv - \iint_S \mathbf{T} \cdot \left[ \hat{n} \times \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) \right] ds
$$

---

Para uma cavidade delimitada por um **condutor elétrico perfeito (PEC)**:

👉 o campo elétrico tangencial na superfície é nulo
👉 a função teste $\mathbf{T}$ também é nula na superfície

Logo, o termo de superfície desaparece.

---

## Forma Final

### (150)

$$
\iiint_V (\nabla \times \mathbf{T}) \cdot \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) dv = k_0^2 \varepsilon_r \iiint_V \mathbf{T} \cdot \mathbf{E} dv
$$

---

# 📄 3.1.2 — Discretização

O volume da cavidade é discretizado utilizando **elementos tetraédricos de primeira ordem**, como o mostrado na Figura 14.

Um tetraedro de primeira ordem possui:

* **4 nós**
* **6 arestas**

As seis arestas são formadas conforme apresentado na Tabela 11.

---

## 📌 Figura (salvar como `fig14_tetrahedron.png`)

```markdown
![Elemento tetraédrico de primeira ordem](images/fig14_tetrahedron.png)
```

Legenda:

**Figura 14. Elemento tetraédrico de primeira ordem.**

---

## 📊 Tabela 11 — Formação das Arestas do Elemento Tetraédrico

| Aresta ( m ) | Nó ( i ) | Nó ( j ) |
| ------------ | -------- | -------- |
| 1            | 1        | 2        |
| 2            | 2        | 3        |
| 3            | 1        | 3        |
| 4            | 1        | 4        |
| 5            | 2        | 4        |
| 6            | 3        | 4        |

---

O campo elétrico em um único elemento tetraédrico é representado por:

### (151)

$$
\mathbf{E} = \sum_{m=1}^{6} e_m \mathbf{W}_m
$$

onde os seis parâmetros desconhecidos associados a cada aresta são:

$$
e_1, e_2, \dots, e_6
$$

O campo total é obtido avaliando a equação (151).

---

## Elementos de Aresta Tangenciais

Os elementos de aresta vetoriais $\mathbf{W}_m$ são dados por:

### (152)

$$
\mathbf{W}*m = L_m \left( \alpha*{ti} \nabla \alpha_{tj} - \alpha_{tj} \nabla \alpha_{ti} \right)
$$

onde:

* $i$ e $j$: nós que definem a aresta $m$
* $L_m$: comprimento da aresta
* $\alpha_{ti}, \alpha_{tj}$: coordenadas simplex associadas aos nós

---

## Coordenadas Simplex

As coordenadas simplex dos nós do tetraedro são:

### (153)–(156)

$$
\alpha_{t1} = \frac{V_1}{V}, \quad \alpha_{t2} = \frac{V_2}{V}, \quad \alpha_{t3} = \frac{V_3}{V}, \quad \alpha_{t4} = \frac{V_4}{V}
$$

---

### Volume do tetraedro

### (157)

$$ 
V = \frac{1}{6} \begin{vmatrix} 1 & x_1 & y_1 & z_1 \ 1 & x_2 & y_2 & z_2 \ 1 & x_3 & y_3 & z_3 \ 1 & x_4 & y_4 & z_4 \end{vmatrix}
$$

---

### Subvolumes

### (158)

$$
V_1 = \frac{1}{6} \begin{vmatrix} 1 & x & y & z \ 1 & x_2 & y_2 & z_2 \ 1 & x_3 & y_3 & z_3 \ 1 & x_4 & y_4 & z_4 \end{vmatrix}
$$

---

### (159)

$$
V_2 = \frac{1}{6} \begin{vmatrix} 1 & x_1 & y_1 & z_1 \ 1 & x & y & z \ 1 & x_3 & y_3 & z_3 \ 1 & x_4 & y_4 & z_4 \end{vmatrix}
$$

---

### (160)

$$
V_3 = \frac{1}{6} \begin{vmatrix} 1 & x_1 & y_1 & z_1 \ 1 & x_2 & y_2 & z_2 \ 1 & x & y & z \ 1 & x_4 & y_4 & z_4 \end{vmatrix}
$$

---

### (161)

$$ 
V_4 = \frac{1}{6} \begin{vmatrix} 1 & x_1 & y_1 & z_1 \ 1 & x_2 & y_2 & z_2 \ 1 & x_3 & y_3 & z_3 \ 1 & x & y & z \end{vmatrix}
$$

Para qualquer nó $i = 1, 2, 3, 4$:

### (162)

$$
\alpha_{ti} = \frac{ a_{ti} + b_{ti} x + c_{ti} y + d_{ti} z }{6V}
$$

onde $a_{ti}, b_{ti}, c_{ti}, d_{ti}$ são coeficientes obtidos dos determinantes de $V_1, V_2, V_3, V_4$.
Substituindo (152) e (162):

### (163)

$$
\mathbf{W}*m = \frac{L_m}{36V^2} \left[ (A*{xm} + B_{xm} y + C_{xm} z)\hat{x} + (A_{ym} + B_{ym} x + C_{ym} z)\hat{y} + (A_{zm} + B_{zm} x + C_{zm} y)\hat{z} \right]
$$

---

### Coeficientes

#### (164)

$$
A_{xm} = a_{ti} b_{tj} - a_{tj} b_{ti}
$$

#### (165)

$$
B_{xm} = c_{ti} b_{tj} - c_{tj} b_{ti}
$$

#### (166)

$$
C_{xm} = d_{ti} b_{tj} - d_{tj} b_{ti}
$$

---

#### (167)

$$
A_{ym} = a_{ti} c_{tj} - a_{tj} c_{ti}
$$

#### (168)

$$
B_{ym} = b_{ti} c_{tj} - b_{tj} c_{ti} = -B_{xm}
$$

#### (169)

$$
C_{ym} = d_{ti} c_{tj} - d_{tj} c_{ti}
$$

---

#### (170)

$$
A_{zm} = a_{ti} d_{tj} - a_{tj} d_{ti}
$$

#### (171)

$$
B_{zm} = b_{ti} d_{tj} - b_{tj} d_{ti} = -C_{xm}
$$

#### (172)

$$
C_{zm} = c_{ti} d_{tj} - c_{tj} d_{ti} = -C_{ym}
$$

---

## Propriedade Fundamental dos Elementos de Aresta

Os elementos $\mathbf{W}_m$ satisfazem:

### (173)

$$
\hat{t}_m \cdot \mathbf{W}_m = \begin{cases} 1, & \text{na aresta } m \ 0, & \text{nas demais arestas} \end{cases}
$$

onde $\hat{t}_m$ é o vetor unitário ao longo da aresta.

---
