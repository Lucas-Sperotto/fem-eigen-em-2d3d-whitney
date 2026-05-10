# §5 — Discussão

> **Status:** rascunho científico para revisão.  
> **Nota de uso:** Este arquivo é o rascunho da Seção 5 do artigo. As subseções
> estão escritas em tom de publicação; revisar antes de submeter.

---

## 5.1 Erro Tipográfico na Equação (120) do Artigo Original

Durante a implementação do módulo de cálculo de *k*₀ para β dado
(Seção 2.2.3 do TP-3485), identificou-se uma inconsistência algébrica na
Eq. (120) do artigo impresso. A equação, que descreve o bloco local
`Sel(tt)` da submatriz de rigidez vetorial, é apresentada no artigo sob a
forma

$$S_{el}^{(tt)}(m,n) = \frac{1}{\mu_r} \frac{L_{tm} L_{tn}}{16 A^3}
\left( D_m D_n + \sum_{k=1}^{5} I_{tk} \right),
\quad [\text{Eq. (120) conforme impressa}]$$

onde *A* é a área do elemento, *L*ₜ são os comprimentos das arestas,
*D* os coeficientes geométricos de curl e *I*ₜₖ os integrais de produto
interno dos funcionais de forma de aresta.

A inconsistência torna-se evidente ao se derivar o bloco local diretamente
das Eqs. (66) e (67) do mesmo artigo, que definem separadamente a
contribuição de curl-curl e a contribuição de massa vetorial da formulação
transversal edge. A derivação elementar mostra que a forma correta é

$$S_{el}^{(tt)}(m,n) = \frac{1}{\mu_r} \frac{L_{tm} L_{tn}}{16 A^3}
\left( 4 D_m D_n + \beta^2 \sum_{k=1}^{5} I_{tk} \right),$$

revelando duas discrepâncias em relação à equação impressa: (i) o fator
numérico 4 no termo de curl-curl, que resulta da integração do produto
de curls de 1-formas de Whitney sobre triângulos com coordenadas baricêntricas;
e (ii) o fator β² multiplicando o somatório *I*ₜₖ, que representa a
contribuição do campo transversal à energia magnética longitudinal e se
anula somente em problemas de corte (β = 0).

A ausência do fator β² na Eq. (120) impressa não é detectável como erro
de consistência pelo leitor que utiliza apenas aquela equação como
referência, uma vez que ambas as formas são idênticas em β = 0, o que
corresponde ao caso de guia homogêneo tratado na Seção 2.2.1 do artigo.
O erro torna-se ativo apenas na Seção 2.2.3, onde β é fornecido como
parâmetro não nulo.

A implementação deste repositório reconstrói o bloco global `A_tt` como

$$A_{tt} = S_t + \beta^2 M_t^{(1/\mu)},$$

em que *S*ₜ é montado pela Eq. (66) e *M*ₜ pela Eq. (67), sem depender
da Eq. (120) impressa. A correção é verificada *a posteriori* pelos resultados
da Tabela 8: os dez modos reportados na Figura 11 do TP-3485 são reproduzidos
com erro máximo de 0,90% em relação à referência analítica e de 2,22% em
relação aos valores tabelados por Hayata *et al.* [REF], confirmando
que o bloco `A_tt` é algebraicamente correto.

---

## 5.2 Comportamento Near-Cutoff na Tabela 10

A Figura 13 do TP-3485 apresenta a dispersão β/*k*₀ em função de *b*r/λ₀
para diferentes valores da razão geométrica *d/a*, onde *d* é a espessura
da camada dielétrica e *a* é a largura do guia. A reprodução computacional
mostra concordância satisfatória com os valores de referência em praticamente
todos os pontos da tabela, com a exceção sistemática dos dois pontos de menor
frequência normalizada (*b*r/λ₀ = 0,5) para *d/a* = 0,167 e *d/a* = 0,286,
cujos erros atingem 34,3% e 37,7%, respectivamente.

Esses pontos correspondem a regimes de operação próximos à frequência de corte
do guia, onde β/*k*₀ → 0. A análise das demais entradas do mesmo *d/a*
mostra que, a partir de *b*r/λ₀ = 0,6, os erros recaem abaixo do gate de
validação de 15% para todos os pares (*d/a*, *b*r/λ₀) testados. O padrão é
consistente: o erro decresce monotonicamente com *b*r/λ₀ para qualquer valor
fixo de *d/a*.

A origem física da dificuldade reside no próprio estrutura do problema de
autovalor. A formulação da Seção 2.2.4 leva ao EVP quadrático em β²:

$$[K]\{e\} + \beta [C]\{e\} + \beta^2 [M]\{e\} = 0,$$

linearizado pela substituição `{f} = β{e}` para produzir o problema generalizado
de dimensão 2*N* × 2*N*

$$\begin{bmatrix} 0 & I \\ -K & -C \end{bmatrix}
\begin{Bmatrix} e \\ f \end{Bmatrix} = \beta
\begin{bmatrix} I & 0 \\ 0 & M \end{bmatrix}
\begin{Bmatrix} e \\ f \end{Bmatrix}.$$

O processo de linearização produz, para cada raiz física β, uma raiz espúria
associada que não corresponde a nenhum modo propagante. Próximo ao corte,
a raiz física β → 0 aproxima-se das raízes espúrias de pequena magnitude
geradas pela expansão do problema, dificultando a separação espectral e
tornando o autovalor sensível à discretização.

Para quantificar esse efeito, a malha utilizada na Seção 2.2.4 contém
100 elementos triangulares (10 × 5 na geometria retangular bipartida),
o que resulta em *N* = 165 graus de liberdade e um sistema linearizado
de 330 × 330. A malha foi escolhida para reproduzir as condições do TP-3485,
que não especifica o tamanho de malha exato para a Figura 13.

Conforme demonstrado na Tabela 9 (Figura 12, mesmo guia retangular, outra
configuração de *d/a*), a formulação é capaz de reproduzir curvas de dispersão
com erros de até 9,2% em *b*r/λ₀ = 0,5, sugerindo que a severidade near-cutoff
depende do valor específico de *d/a*: a camada dielétrica mais fina
(*d/a* = 0,167 e *d/a* = 0,286) concentra o modo próximo da interface e
exige maior resolução lateral para capturar a distribuição do campo.

A redução do erro nesses pontos por refinamento de malha é atualmente uma
hipótese em verificação. Os resultados disponíveis sugerem que o aumento
da razão de aspecto da malha (mais divisões na direção transversal à
interface dielétrica) produz ganho mais expressivo do que o refinamento
uniforme. Esta análise será completada em trabalho futuro.

---

## 5.3 EFGMI como Extensão Original

O método de Galerkin Sem Malha com Integração Modificada (EFGMI, do inglês
*Element-Free Galerkin Method with Modified Integration*) foi introduzido
neste repositório como extensão do arcabouço FEM para funções de base meshfree.
A motivação é explorar se a flexibilidade das funções de forma MLS (*Moving
Least Squares*) pode produzir autovalores de qualidade comparável ao FEM
com a mesma grade de pontos, especialmente em geometrias de interface irregular.

### 5.3.1 Metodologia Implementada

As funções de forma EFGMI são construídas pelo procedimento de mínimos
quadrados móvel com base polinomial linear (ordem de consistência 1).
Para um ponto de avaliação x, a função de forma φᵢ(x) associada ao nó *i*
é obtida pela minimização ponderada

$$\sum_{j \in S(x)} w(x - x_j) \left[ p(x_j)^T a - u_j \right]^2,$$

onde `S(x)` é o conjunto de suporte de x (todos os nós dentro de um raio
*h* = 2,5 × escala nodal local), `p(x)` é o vetor de monômios lineares
`[1, x, y]^T` e *w*(·) é uma função de peso cúbica spline.

A malha triangular é usada exclusivamente como grade de integração: as
quadraturas de Gauss são realizadas em cada elemento, mas os pontos de
integração utilizam os φᵢ EFGMI — e não as funções de forma de Lagrange
— para montar as matrizes globais. Essa separação entre grade de integração
e base de aproximação é a característica distintiva do método.

A base vetorial "Whitney-like" é construída a partir das funções nodais
EFGMI por analogia direta com a 1-forma de Whitney:

$$\mathbf{W}_{ij}^{\text{EFGMI}}(\mathbf{x}) =
\phi_i(\mathbf{x}) \nabla \phi_j(\mathbf{x}) -
\phi_j(\mathbf{x}) \nabla \phi_i(\mathbf{x}),$$

associada a cada aresta (*i*, *j*) da triangulação de fundo. Essa
construção preserva a propriedade de tangência contínua nas arestas,
necessária para a representação correta de campos TE.

O custo de suporte mínimo adotado é de 6 nós por ponto de avaliação;
interfaces bimaterial utilizam raio ampliado (escala 4,0 × comprimento
nodal) com mínimo de 9 nós, para garantir estabilidade da inversão da
matriz de momentos.

### 5.3.2 Acurácia Comparativa

Os autovalores obtidos pelo backend EFGMI são comparáveis aos do FEM na
mesma grade, como mostra a Tabela [N]. Para os casos escalares (Tabelas 1–3),
o erro EFGMI excede o FEM em 0,1–2,5 pontos percentuais. Para as formulações
mistas (Tabela 6), a diferença é praticamente nula nas componentes edge — que
são idênticas entre FEM e EFGMI por compartilharem a mesma base vetorial
de suporte — e ligeiramente superior nas componentes escalares nodais.

A semelhança de acurácia entre os dois backends corrobora a implementação:
a malha de fundo é suficientemente refinada para que a quadratura seja a
fonte dominante de erro em ambos os casos.

### 5.3.3 Custo de Montagem e Limitações

O custo de assembly do EFGMI é sistematicamente superior ao do FEM em razão
do cálculo das funções de forma MLS por ponto de quadratura: enquanto o FEM
avalia as funções de forma por lookup de um estêncil local de 3 nós (O(1)
por ponto), o EFGMI resolve um sistema 3 × 3 para cada ponto de Gauss
dentro do conjunto de suporte S(x), que contém tipicamente 6–12 nós.
Esse overhead resulta em razões de assembly de 11 a 940 vezes maiores do
que o FEM nas geometrias testadas (Tabela [N]).

O tempo de solução (LAPACK `dsygv`) é essencialmente idêntico entre os
dois backends, uma vez que ambos produzem matrizes densas de mesma dimensão.

As limitações atuais da extensão EFGMI são: (i) implementação restrita a
problemas 2D; (ii) ordem de consistência 1 (base de monômios lineares),
que limita a taxa de convergência dos autovalores a *O*(*h*²), a mesma dos
elementos P1; (iii) ausência de referência bibliográfica direta para a
aplicação de EFGM em eigenproblemas de guias de onda eletromagnéticos —
o método foi originalmente proposto para elasticidade linear [Belytschko
et al., 1994] e a extensão para EM é contribuição original deste trabalho.

### 5.3.4 Perspectivas

A extensão 3D do EFGMI para cavidades tetraédricas é conceitualmente
direta: a base Whitney-like tridimensional pode ser construída de modo
análogo a partir das funções MLS nodais sobre tetraedros de fundo.
O principal desafio é o custo de assembly, que em 3D escala com o volume
do suporte e cresce mais rapidamente do que em 2D. Bases de ordem superior
(quadráticas) poderiam melhorar a acurácia sem aumento proporcional de
custo, e estão previstas como trabalho futuro.

---

> **Notas para o autor:**
>
> - Substituir `[REF]` pela referência de Hayata et al. (1986) e adicionar
>   DOI quando disponível.
> - Tabela [N] na seção EFGMI deve ser numerada sequencialmente no artigo final.
> - As equações em LaTeX precisam ser convertidas para o formato do template
>   IEEE antes da submissão.
> - A frase "Esta análise será completada em trabalho futuro" na §5.2 deve
>   ser substituída pelo resultado real do experimento de convergência (T2
>   do plano de testes) antes da submissão.
