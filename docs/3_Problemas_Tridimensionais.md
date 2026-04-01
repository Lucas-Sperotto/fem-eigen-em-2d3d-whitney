# 📄 3. Problemas Tridimensionais

## 3.1 Autovalores de Cavidades Tridimensionais — Formulação Vetorial

O problema de calcular frequências de ressonância em cavidades tridimensionais tem sido historicamente dificultado pela presença de **modos espúrios**.

Como mencionado na seção 2.2, esse problema foi recentemente superado utilizando **elementos de aresta tangenciais vetoriais** (refs. 15 e 16).

Nesta seção, formula-se o método de elementos finitos de Galerkin para cavidades tridimensionais e apresentam-se resultados para diversas geometrias.

---

## Leitura orientadora desta seção

O leitor deve entrar nesta seção com uma expectativa correta: o salto para 3D não muda a filosofia do método, mas amplia seu cenário geométrico e algébrico. Em 2D, já havíamos aprendido que a escolha adequada dos espaços discretos é decisiva para evitar modos espúrios. Em 3D, essa mesma lição retorna com mais força, porque agora o campo vive em volume, a malha é tetraédrica e a continuidade tangencial precisa ser respeitada em todas as arestas do elemento.

Em termos de formação, esta é uma das partes mais nobres do artigo. Ela mostra que a teoria não foi construída para um caso isolado, mas para sobreviver ao aumento de dimensão sem perder a coerência física.

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

#### Passagem intermediária: a integração por partes vetorial em 3D

Esta é a etapa central da formulação variacional. Se definirmos

$$
\mathbf{B} = \frac{1}{\mu_r}\nabla \times \mathbf{E},
$$

então a identidade (145) com $\mathbf{A} = \mathbf{T}$ produz imediatamente

$$
\mathbf{T}\cdot \nabla \times \left(\frac{1}{\mu_r}\nabla \times \mathbf{E}\right) = (\nabla \times \mathbf{T})\cdot \left(\frac{1}{\mu_r}\nabla \times \mathbf{E}\right) - \nabla \cdot \left[\mathbf{T}\times \left(\frac{1}{\mu_r}\nabla \times \mathbf{E}\right)\right].
$$

Ao integrar essa identidade em todo o volume, o primeiro termo gera a parte interna da forma fraca e o segundo termo prepara o nascimento da contribuição de contorno. É exatamente a versão tridimensional da integração por partes que, em problemas escalares, costuma aparecer escondida sob o nome de fórmula de Green.

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

#### Passagem intermediária: por que a condição PEC fecha a forma fraca

Em uma cavidade PEC, a condição física essencial é

$$
\hat{n}\times \mathbf{E} = 0 \quad \text{em } S.
$$

No método de Galerkin, a função teste pertence ao mesmo espaço admissível do campo aproximado, de modo que também satisfaz

$$
\hat{n}\times \mathbf{T} = 0 \quad \text{em } S.
$$

Por isso, o termo de contorno de (149) deixa de contribuir. Em linguagem funcional, a formulação passa a viver em um subespaço de $H(\mathrm{curl},V)$ com traço tangencial homogêneo no contorno. Essa observação é importante porque não se trata de um cancelamento acidental, mas de uma consequência direta das condições físicas impostas à cavidade.

## Forma Final

### (150)

$$
\iiint_V (\nabla \times \mathbf{T}) \cdot \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) dv = k_0^2 \varepsilon_r \iiint_V \mathbf{T} \cdot \mathbf{E} dv
$$

---

## 📄 3.1.2 — Discretização

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

#### Passagem intermediária: o que representam os seis coeficientes de aresta

Os coeficientes $e_m$ não devem ser lidos como simples "valores do campo em pontos". Esse seria um pensamento nodal, adequado a problemas escalares, mas insuficiente aqui. No espírito dos elementos de aresta, cada grau de liberdade está associado à circulação tangencial do campo ao longo de uma aresta:

$$
e_m \sim \int_{l_m} \mathbf{E}\cdot \hat{t}_m\, dl.
$$

Essa é a chave geométrica do método. O campo é reconstruído a partir de quantidades tangenciais integradas ao longo das arestas, o que preserva precisamente a continuidade física que interessa em eletromagnetismo.

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

#### Passagem intermediária: por que as coordenadas simplex são tão poderosas

As coordenadas simplex organizam toda a geometria do tetraedro em uma linguagem afim muito econômica. Elas satisfazem, em qualquer ponto interno do elemento,

$$
\alpha_{t1} + \alpha_{t2} + \alpha_{t3} + \alpha_{t4} = 1,
$$

e cada $\alpha_{ti}$ vale $1$ no nó $i$ e $0$ nos outros três nós. Com isso, o tetraedro inteiro passa a ser descrito por funções lineares simples, cujos gradientes são constantes no interior do elemento. É essa linearidade que torna possível obter fórmulas fechadas para muitos blocos da formulação.

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

#### Passagem intermediária: de $\alpha_{ti}$ para seus gradientes

Como $\alpha_{ti}$ é linear em $x$, $y$ e $z$, seu gradiente é constante dentro do tetraedro:

$$
\nabla \alpha_{ti} = \frac{1}{6V} \left( b_{ti}\hat{x} + c_{ti}\hat{y} + d_{ti}\hat{z} \right).
$$

Essa observação é preciosa. Ela mostra por que os elementos de aresta de primeira ordem conseguem produzir fórmulas explícitas: os objetos básicos da interpolação já nascem com dependência afim simples, e seus gradientes não variam de ponto para ponto no interior do elemento.

Substituindo (152) e (162):

### (163)

$$
\mathbf{W}*m = \frac{L_m}{36V^2} \left[ (A*{xm} + B_{xm} y + C_{xm} z)\hat{x} + (A_{ym} + B_{ym} x + C_{ym} z)\hat{y} + (A_{zm} + B_{zm} x + C_{zm} y)\hat{z} \right]
$$

#### Passagem intermediária: como nasce a forma explícita da base de aresta

Se substituirmos a expressão linear de $\alpha_{ti}$ e o gradiente constante de $\alpha_{ti}$ em

$$
\mathbf{W}_m = L_m\left(\alpha_{ti}\nabla \alpha_{tj} - \alpha_{tj}\nabla \alpha_{ti}\right),
$$

obtemos um produto entre termos lineares e gradientes constantes. Em forma expandida, isso gera componentes afins em $x$, $y$ e $z$, todas divididas por $(6V)(6V) = 36V^2$:

$$
\mathbf{W}_m = \frac{L_m}{36V^2} \left[ \bigl(\text{termo constante} + \text{termo em } y + \text{termo em } z\bigr)\hat{x} + \bigl(\text{termo constante} + \text{termo em } x + \text{termo em } z\bigr)\hat{y} + \bigl(\text{termo constante} + \text{termo em } x + \text{termo em } y\bigr)\hat{z} \right].
$$

Os coeficientes $A_{\bullet m}$, $B_{\bullet m}$ e $C_{\bullet m}$ que aparecem em (164)–(172) são justamente o resultado desse agrupamento algébrico. Em outras palavras, a equação (163) não cai do céu: ela é a forma reorganizada da definição geométrica da base de Whitney.

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

#### Passagem intermediária: interpretação física da propriedade de aresta

A mensagem de (173) é profunda: cada função de base "reconhece" a sua própria aresta e não interfere, no sentido interpolatório, nas demais. Uma forma operacional, muito usada na literatura de elementos de aresta, é escrever essa propriedade por integrais de linha:

$$
\int_{l_n} \mathbf{W}_m \cdot \hat{t}_n\, dl = \delta_{mn}.
$$

Essa leitura ajuda a evitar uma confusão comum. O que caracteriza o grau de liberdade não é um valor pontual qualquer do campo, mas sua circulação tangencial associada à aresta. É justamente isso que alinha a discretização com a física de Maxwell.

## Fecho pedagógico da seção tridimensional inicial

Até a equação (173), o artigo constrói os alicerces 3D: a forma fraca para cavidades PEC, a interpolação por arestas tetraédricas e a maquinaria geométrica das coordenadas simplex. É uma abertura de grande valor didático, porque ensina ao leitor que o sucesso do método em 3D não depende de truques numéricos, mas de uma escolha correta de variáveis, espaços funcionais e graus de liberdade geométricos.

## Notas de revisão complementar

### Papel desta seção no desenvolvimento do artigo

O salto para 3D não abandona a lógica construída em 2D. O artigo continua trabalhando com uma formulação variacional e com graus de liberdade associados a arestas; a diferença é que agora tudo é transportado para tetraedros, com seis arestas por elemento e coordenadas simplex em volume.

### Encadeamento conceitual até a equação (173)

- As equações (143) a (150) fazem em 3D o mesmo papel estrutural que as primeiras equações de 2.2.1 fizeram em 2D: sair da forma forte e chegar à forma fraca com as condições de contorno apropriadas.
- As equações (151) e (152) introduzem a expansão do campo em bases de aresta tetraédricas.
- As equações (153) a (162) constroem a parametrização geométrica do tetraedro por coordenadas simplex.
- As equações (163) a (172) tornam explícita a forma afim das funções de base e os coeficientes geométricos que alimentam as matrizes elementares.
- A equação (173) registra a propriedade fundamental de interpolação tangencial nas arestas, que é a peça conceitual central para evitar modos espúrios.

### Passagens que merecem leitura mais cuidadosa

- A sequência (145) -> (149) -> (150) é a parte variacional mais importante desta seção, porque é nela que o termo de contorno é manipulado até restar a forma fraca compatível com PEC.
- A passagem de (152) para (163) merece atenção especial: ali a base de Whitney deixa de ser apenas uma definição abstrata e passa a ter forma explícita em função das coordenadas do tetraedro.
- As relações `B_{ym} = -B_{xm}`, `B_{zm} = -C_{xm}` e `C_{zm} = -C_{ym}` são úteis para verificar coerência algébrica entre os coeficientes auxiliares.

### Pontos de conferência e irregularidades prováveis

- Na equação (146), o termo volumétrico do lado direito termina com `ds`, embora o contexto seja uma integral em volume. A leitura mais natural, pelo desenvolvimento matemático, parece ser `dv`.
- As equações (152) e (163) apresentam marcas de OCR como `\mathbf{W}*m` e `\alpha*{ti}`. A leitura correta parece ser sem os asteriscos.
- As equações determinantes (157) a (161) foram comprimidas em uma única linha no Markdown atual. Isso preserva a informação, mas dificulta bastante a inspeção visual do sinal e da ordem das colunas.
- A equação (173) está semanticamente correta, mas o ambiente `cases` ficou sem quebra de linha explícita entre os dois ramos, o que atrapalha a renderização.
- A propriedade registrada em (173) merece leitura conceitual cuidadosa: em formulações de aresta, sua versão mais robusta costuma ser expressa por integral de linha ao longo da aresta, não apenas por leitura pontual da componente tangencial.
- A Figura 14 aparece como placeholder, então o apoio visual para o tetraedro de primeira ordem ainda está incompleto.

### Observação estrutural importante

- Este arquivo, no estado atual, encerra o desenvolvimento em (173). A continuação natural da teoria 3D, com as equações (174) a (182) e os exemplos numéricos, aparece hoje em outro arquivo da pasta. Não removi nem movi nada; estou apenas registrando essa observação para guiar a leitura.
