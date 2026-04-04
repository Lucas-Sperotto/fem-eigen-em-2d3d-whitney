# HELM10: Referencia dos Campos do CSV de Campos

Este arquivo documenta o significado didatico das colunas do CSV de campos
nodais gerado pelos executaveis da familia `HELM10`:

- `helm10_rect`
- `helm10_circle`
- `helm10_coax`

Ele deve ser lido junto com:

- a traducao da Secao 2.1 do artigo em
  [docs/traducao/2.1_Guias de Onda Homogeneos.md](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md)
- o glossario do CSV de modos em
  [HELM10_CSV_Modos_Referencia.md](HELM10_CSV_Modos_Referencia.md)
- a documentacao de implementacao em
  [src/helm10/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helm10/README.md)

## Ideia central

Cada arquivo `fields_<modo>.csv` representa um unico modo.

Por isso, os metadados do modo nao se repetem em todas as linhas do arquivo de
campos. Eles ficam concentrados no `modes.csv`, que aponta para o arquivo
correspondente pela coluna `fields_csv_file`.

Em outras palavras:

- `modes.csv` responde "qual modo e este?"
- `fields_<modo>.csv` responde "quais sao os valores nodais deste modo?"

## Cabecalho-base

O CSV de campos usa hoje o cabecalho:

```text
node_id,x_m,y_m,psi,dpsi_dx,dpsi_dy,Ex,Ey
```

Nos casos circular e coaxial, o cabecalho nodal continua o mesmo, porque a
tabela registra valores nos nos da malha, nao parametros geometricos globais.

## Significado de cada coluna

- `node_id`
  - indice do no da malha onde a linha foi avaliada.
  - permite correlacionar a linha com a malha usada pelo solve e com o VTK
    correspondente.

- `x_m`
  - coordenada `x` do no, em metros.

- `y_m`
  - coordenada `y` do no, em metros.

- `psi`
  - valor nodal do potencial longitudinal associado ao modo.
  - no uso atual do projeto:
    - `TE -> psi` associado a `Hz`
    - `TM -> psi` associado a `Ez`

- `dpsi_dx`
  - derivada parcial nodal de `psi` em relacao a `x`.
  - e obtida a partir do gradiente por elemento, reunido de forma didatica
    para permitir reconstruir o campo transversal no no.

- `dpsi_dy`
  - derivada parcial nodal de `psi` em relacao a `y`.

- `Ex`
  - componente `x` do campo eletrico transversal exportado pelo projeto.
  - o valor salvo e sempre o campo sem multiplicacao por impedancia modal.
  - no uso atual:
    - TE: segue a convencao por gradiente com sinal do projeto
    - TM: segue a convencao sem `Ztm`

- `Ey`
  - componente `y` do campo eletrico transversal exportado pelo projeto.
  - assim como `Ex`, o valor salvo e sempre sem multiplicacao por impedancia
    modal.

## Convencao operacional atual

Para evitar ambiguidade didatica, o CSV de campos salva apenas uma versao de
campo transversal:

- TE: `Ex` e `Ey` conforme a reconstrucao transversal adotada no projeto
- TM: `Ex` e `Ey` sem multiplicacao por `Ztm`

Assim, o arquivo de campos nao mistura duas normalizacoes concorrentes.

As grandezas auxiliares que descrevem o modo como um todo, como:

- `freq_hz`
- `eps_r`
- `mu_r`
- `k_medium`
- `kc_fem`
- `kc_ana`
- `beta`
- `ztm`
- `below_cutoff`
- `error_percent`
- `rho_abs`
- `field_status`

ficam no `modes.csv`, porque sao constantes para aquele modo inteiro e nao
precisam ser repetidas linha a linha no arquivo nodal.

## Relacao com a Secao 2.1 do artigo

Este CSV e a forma operacional de registrar a etapa final da Secao 2.1:

- primeiro o problema escalar calcula o cutoff `kc`
- depois o potencial longitudinal `psi` e obtido nos nos
- por fim o gradiente transversal permite reconstruir `Ex` e `Ey`

Logo, este arquivo nao introduz nova teoria. Ele apenas organiza, em formato
tabular, o potencial e o campo transversal reconstruido.

## Como ler junto com o CSV de modos

O `modes.csv` responde:

- qual modo foi identificado
- qual e o seu cutoff
- se ele esta acima ou abaixo do corte
- qual foi o erro contra a referencia
- qual e o nome do arquivo de campos correspondente

O `fields_<modo>.csv` responde:

- qual e o valor de `psi` em cada no
- quais sao `dpsi_dx` e `dpsi_dy`
- quais sao `Ex` e `Ey` ponto a ponto

## Referencias cruzadas

- teoria traduzida:
  [2.1 Guias de Onda Homogeneos](traducao/2.1_Guias%20de%20Onda%20Homogêneos.md)
- glossario do CSV de modos:
  [HELM10_CSV_Modos_Referencia.md](HELM10_CSV_Modos_Referencia.md)
- implementacao:
  [src/helm10/README.md](/home/sperotto/tp3485-fem-eigen-em/src/helm10/README.md)
- visao operacional dos executaveis:
  [Tabela_Executaveis_Entradas_Saidas.md](Tabela_Executaveis_Entradas_Saidas.md)
