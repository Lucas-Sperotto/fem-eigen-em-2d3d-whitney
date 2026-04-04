# HELM10: Referencia dos Campos do CSV de Modos

Este arquivo documenta o significado didatico das colunas do CSV de modos
gerado pelos executaveis da familia `HELM10`:

- `helm10_rect`
- `helm10_circle`
- `helm10_coax`

Ele deve ser lido junto com:

- a traducao da Secao 2.1 do artigo em
  [docs/traducao/2.1_Guias de Onda Homogeneos.md](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md)
- a documentacao de implementacao em
  [src/helm10/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helm10/README.md)

## Cabecalho-base do caso retangular

O CSV de modos do caso retangular usa hoje o cabecalho:

```text
family,longitudinal_label,mode_label,positive_rank,eig_index,m,n,ar_m,b_m,freq_hz,eps_r,mu_r,k_medium,kc_fem,kc_ana,kc_ar_fem,kc_ar_ana,beta,ztm,below_cutoff,error_percent,rho_abs,field_status,fields_csv_file
```

Nos casos circular e coaxial, a ideia fisica e a mesma, mas os campos
geometricos mudam:

- circular: aparece `r_m` no lugar de `ar_m,b_m`
- coaxial: aparecem `r1_m,r2_m`
- circular e coaxial usam `p` no lugar de `n`

## Significado de cada coluna

- `family`
  - familia modal identificada no problema escalar.
  - no uso atual do projeto, os valores tipicos sao `TE` e `TM`.

- `longitudinal_label`
  - componente longitudinal associada a incognita escalar.
  - `TE -> Hz`
  - `TM -> Ez`

- `mode_label`
  - rotulo textual do modo ja identificado e pronto para uso em nomes de
    arquivos.
  - exemplos: `TE_m1_n0`, `TM_m2_n1`, `TE_m1_p1`.

- `positive_rank`
  - ordem didatica do modo dentro do conjunto realmente exportado.
  - conta apenas modos fisicos positivos da familia considerada.
  - exemplo: o primeiro modo TE salvo recebe `positive_rank = 1`.

- `eig_index`
  - posicao do autovalor/autovetor no espectro retornado pelo solver.
  - diferente de `positive_rank`: aqui a contagem segue o indice interno do
    problema generalizado.

- `m`
  - primeiro indice modal associado ao modo identificado por casamento com a
    referencia analitica.

- `n`
  - segundo indice modal do caso retangular.
  - no caso circular e coaxial, esse papel e exercido por `p`.

- `ar_m`
  - dimensao `a` do guia retangular, em metros.

- `b_m`
  - dimensao `b` do guia retangular, em metros.

- `r_m`
  - raio do guia circular, em metros.

- `r1_m`
  - raio interno da linha coaxial, em metros.

- `r2_m`
  - raio externo da linha coaxial, em metros.

- `freq_hz`
  - frequencia usada na reconstruicao fisica auxiliar.
  - se `--freq-hz` for informado, o valor vem do usuario.
  - se `--freq-hz` for omitido, o codigo escolhe automaticamente uma
    frequencia a partir dos `kc` exportados.
  - quando ha modos TM exportados, a politica automatica tenta usar o modo TM
    limitante como referencia e impor `Ztm = 1` para ele.

- `eps_r`
  - permissividade relativa do meio homogeneo usado na reconstrucao auxiliar.

- `mu_r`
  - permeabilidade relativa do meio homogeneo usado na reconstrucao auxiliar.

- `k_medium`
  - numero de onda do meio homogeneo.
  - no codigo, ele satisfaz `k = omega * sqrt(mu * eps)`.

- `kc_fem`
  - numero de onda de corte calculado numericamente pelo FEM.
  - como o autovalor do problema escalar e `lambda = kc^2`, este valor vem de
    `sqrt(lambda)`.

- `kc_ana`
  - numero de onda de corte analitico associado ao modo identificado.
  - serve como referencia para comparacao e erro.

- `kc_ar_fem`
  - forma adimensional do cutoff numerico no caso retangular.
  - no projeto atual, `kc_ar_fem = kc_fem * a`.

- `kc_ar_ana`
  - forma adimensional do cutoff analitico no caso retangular.
  - no projeto atual, `kc_ar_ana = kc_ana * a`.

- `beta`
  - constante de propagacao longitudinal inferida a partir de `k` e `kc`.
  - quando o modo esta abaixo do corte, o projeto registra `beta = 0`.

- `ztm`
  - impedancia modal TM auxiliar.
  - so faz sentido fisico para o ramo TM.
  - no ramo TE, fica tipicamente em `0`.
  - no ramo TM abaixo do corte, tambem fica em `0`.
  - mesmo quando `ztm` e calculado, ele nao altera o campo principal salvo em
    `Ex` e `Ey` nos CSVs de campos.

- `below_cutoff`
  - indicador logico de corte.
  - `0`: modo propagante para a frequencia escolhida.
  - `1`: modo abaixo do corte.

- `error_percent`
  - erro percentual entre `kc_fem` e `kc_ana`.
  - usado como indicador numerico de proximidade com a referencia analitica.

- `rho_abs`
  - valor absoluto da correlacao usada no casamento modal.
  - mede o quao bem o autovetor numerico combina com o modo analitico
    candidato.
  - valores mais proximos de `1` indicam casamento mais limpo.

## Como `rho_abs` e calculado e por que funciona

No `HELM10`, o casamento modal nao compara apenas autovalores. Ele compara a
forma espacial do autovetor FEM com a forma espacial de um candidato
analitico, ambos escritos sobre os mesmos graus de liberdade nodais.

No codigo, a correlacao usa a matriz global de massa `T` da formulacao
escalar:

- Eq. `(33)`: massa consistente local do triangulo linear
- Eq. `(34)`: montagem global do problema
  `S psi = kc^2 T psi`

Essa mesma matriz `T` define o produto interno discreto usado no matching:

```text
<x, y>_T = x^T T y
||x||_T = sqrt(x^T T x)
rho     = |<x, y>_T| / (||x||_T ||y||_T)
```

No repositório, isso aparece diretamente em:

- [mode_match_rect.hpp](/home/sperotto/tp3485-fem-eigen-em/src/core/mode_match_rect.hpp)
- [mode_match_circle.hpp](/home/sperotto/tp3485-fem-eigen-em/src/core/mode_match_circle.hpp)
- [mode_match_coax.hpp](/home/sperotto/tp3485-fem-eigen-em/src/core/mode_match_coax.hpp)

### Leitura matemática

- `x` e o autovetor numerico produzido pelo solve FEM.
- `y` e o modo analitico candidato, avaliado nos nos da mesma malha.
- `x^T T y` aproxima a integral do produto dos dois campos sobre a secao
  transversal.
- `||x||_T` e `||y||_T` removem a dependencia da escala dos vetores.
- o valor absoluto remove a ambiguidade de sinal global do autovetor, porque
  `x` e `-x` representam o mesmo modo fisico.

### Por que isso funciona

- A matriz de massa consistente `T` representa, no espaco discreto, o produto
  interno natural do problema escalar.
- Portanto, `rho` funciona como um "cosseno" entre dois modos no espaco
  induzido por `T`.
- Se o autovetor FEM e o candidato analitico tiverem a mesma distribuicao
  espacial, entao `rho` fica proximo de `1`.
- Se eles tiverem formas diferentes, cancelamentos no produto interno fazem
  `rho` diminuir.

Em linguagem didatica: o modo certo e aquele cuja forma "se sobrepoe melhor"
ao autovetor numerico quando a sobreposicao e medida com o peso fisicamente
correto da formulacao escalar.

### Casos circulares e coaxiais

Nos casos circular e coaxial, ha familias degeneradas em que os parceiros
`cos(m theta)` e `sin(m theta)` possuem o mesmo cutoff analitico. Nesses
casos, o codigo testa os candidatos compativeis e guarda o maior valor de
correlacao.

Por isso, nesses problemas:

- `rho_abs` continua sendo um bom indicador de identificacao modal;
- mas ele deve ser lido junto com `kc_fem`, `kc_ana` e o label modal final.

### Limites praticos

- Se a malha estiver muito grossa, `rho` pode cair mesmo para o modo correto.
- Em degenerescencias muito proximas, dois candidatos podem aparecer com
  correlacoes parecidas.
- Por isso, o projeto usa `rho_abs` como criterio principal de matching, mas
  sempre o registra ao lado do cutoff identificado e do erro correspondente.

- `field_status`
  - estado textual da reconstrucao de campos.
  - exemplos:
    - `te_sign_only_above_cutoff`
    - `tm_sign_only_above_cutoff`
    - `tm_sign_only_below_cutoff`
  - essa coluna condensa a leitura operacional mais importante:
    acima/abaixo do corte e politica usada no salvamento dos campos.

- `fields_csv_file`
  - nome do arquivo CSV de campos nodais associado a esse modo.
  - permite navegar do resumo modal para o detalhe ponto a ponto.

## Relacao com a Secao 2.1 do artigo

As colunas acima nao sao "novos simbolos" do artigo. Elas juntam, em formato
operacional:

- a identificacao modal obtida pelo solve da formulacao escalar da Secao 2.1;
- a comparacao com a referencia analitica dos casos retangular, circular e
  coaxial;
- a reconstrucao dos campos transversais a partir do potencial longitudinal;
- o estado fisico do modo para a frequencia escolhida.

Em outras palavras:

- o artigo fornece a formulacao e as equacoes da reconstrucao;
- o CSV organiza isso como registro de execucao e auditoria didatica.

## Como ler junto com o CSV de campos

O CSV de modos responde:

- qual modo foi identificado
- qual e o seu cutoff
- qual foi o erro contra a referencia
- em que estado fisico ele foi reconstruido
- qual arquivo contem os campos nodais daquele modo

O CSV de campos responde:

- qual e o valor de `psi` em cada no
- quais sao `dpsi_dx` e `dpsi_dy`
- quais sao `Ex` e `Ey`
- sempre sem multiplicacao por impedancia modal

Para o detalhe coluna por coluna do arquivo nodal, ver:

- [HELM10_CSV_Campos_Referencia.md](HELM10_CSV_Campos_Referencia.md)

## Referencias cruzadas

- teoria traduzida:
  [2.1 Guias de Onda Homogeneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md)
- implementacao:
  [src/helm10/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helm10/README.md)
- visao operacional dos executaveis:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
