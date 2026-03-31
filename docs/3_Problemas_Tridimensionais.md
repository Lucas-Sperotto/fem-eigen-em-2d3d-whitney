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
