# Secao 2.2.1 - Guia Homogeneo com Campos Vetoriais Transversais

## Escopo

Esta secao introduz a formulacao vetorial 2D com elementos de aresta para eliminar modos espurios. No repositorio, ela corresponde ao programa `HELMVEC`.

## Ideia central

O artigo faz tres movimentos aqui:

1. troca a formulacao nodal vetorial por elementos de aresta
2. impõe apenas continuidade tangencial entre elementos
3. constroi o problema de autovalor diretamente para o campo transversal

## Revisao das equacoes

- Eq. (44): forma integral inicial para o problema vetorial transversal.
- Eq. (45): identidade vetorial `T_t . (curl A)`.
- Eq. (46): identidade de produto vetorial no contorno.
- Eq. (47): teorema da divergencia em 2D para o termo de contorno vetorial.
- Eq. (48): forma fraca vetorial ainda com integral de contorno.
- Eq. (49): forma fraca final em contorno PEC.
- Eq. (50): expansao de `E_t` nas funcoes de aresta do triangulo.
- Eq. (51): definicao da base de Whitney 2D `W_tm`.
- Eq. (52): forma explicita da base na aresta 1.
- Eq. (53): forma explicita da base na aresta 2.
- Eq. (54): forma explicita da base na aresta 3.
- Eq. (55): interpretacao fisica do coeficiente de aresta como componente tangencial.
- Eq. (56): propriedade `div_t W_tm = 0`.
- Eq. (57): definicao do coeficiente `A_m`.
- Eq. (58): definicao do coeficiente `B_m`.
- Eq. (59): definicao do coeficiente `C_m`.
- Eq. (60): definicao do coeficiente `D_m`.
- Eq. (61): equacao elemental antes da forma matricial.
- Eq. (62): forma matricial elemental do problema vetorial.
- Eq. (63): matriz elemental de rigidez vetorial.
- Eq. (64): matriz elemental de massa vetorial.
- Eq. (65): problema global `S e_t = kc^2 T e_t`.
- Eq. (66): forma fechada do bloco `curl-curl`.
- Eq. (67): forma fechada do bloco de massa vetorial.
- Eq. (68) a Eq. (72): integrais auxiliares do termo de massa.
- Eq. (73) a Eq. (77): reducoes geometricas dos momentos do triangulo.

## Resultado teorico principal

O ponto de chegada da secao e:

```text
S e_t = kc^2 T e_t
```

com bases de aresta:

```text
W_tm = L_tm (alpha_i grad_t alpha_j - alpha_j grad_t alpha_i)
```

e formas locais fechadas:

```text
S_e(m,n) = (1/mu_r) * (L_m L_n)/(4 A^3) * D_m D_n
T_e(m,n) = eps_r * (L_m L_n)/(16 A^3) * (It1 + It2 + It3 + It4 + It5)
```

## Por que esta secao e crucial

- Ela justifica matematicamente a troca de graus de liberdade nodais por arestas.
- Ela introduz a notacao `A_m`, `B_m`, `C_m`, `D_m`, `It1..It5`.
- Ela fornece o bloco elemental que sera reaproveitado nas secoes 2.2.2, 2.2.3 e 2.2.4.

## Figuras e tabelas desta secao

- Figura 10: configuracao dos elementos de aresta tangenciais
- Tabela 4: guia retangular
- Tabela 5: guia circular

## Leitura para o codigo

- Eq. (66) e Eq. (67) aparecem diretamente no backend `closed-form`.
- Eq. (65) explica a organizacao central de `src/helmvec`.
- As integrais auxiliares `It1..It5` sao a ponte entre o artigo e `src/explicit/tri2d_edge_explicit.hpp`.

## Ajustes em relacao a documentacao legada

- A documentacao antiga preservava bem a estrutura de 2.2.1, mas deixava a relacao entre Eq. (65), Eq. (66) e Eq. (67) pouco explicita.
- A legenda da Figura 10 estava incompleta; aqui ela foi tratada como referencia do artigo, nao como recurso grafico totalmente integrado ao Markdown legado.

