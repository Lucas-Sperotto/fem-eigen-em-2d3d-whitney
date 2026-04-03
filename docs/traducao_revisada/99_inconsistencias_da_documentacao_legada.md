# Inconsistencias da Documentacao Legada

Este arquivo registra o que foi encontrado nos arquivos antigos de `docs/` e justifica a criacao da pasta `docs/traducao_revisada/`.

Parte dessas inconsistencias ja foi corrigida na trilha principal atual em `docs/`. O objetivo deste registro e historico: preservar o diagnostico da documentacao legada e o motivo das revisoes posteriores.

## 1. Problemas estruturais graves

- `docs/2.2.3_Determinação do número de onda.md`
  - comeca com a secao `2.2.2.5` em vez de comecar em `2.2.3`
  - mistura trechos de `2.2.3` com partes inteiras de `3.1.3`, `3.1.4`, `3.1.5` e `3.1.6`
  - deixou de ser um arquivo de secao e virou um agregado incoerente
- `docs/3_Problemas_Tridimensionais.md`
  - termina na Eq. (173) e nao fecha a derivacao ate Eq. (182)
  - portanto nao documenta o problema 3D completo do artigo
- `docs/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md`
  - preserva a ordem geral da secao, mas com forte contaminacao de OCR e formulas danificadas

## 2. Problemas de numeracao e escopo

- A raiz de `docs/` nao seguia mais o sumario real do artigo.
- O sumario oficial do PDF separa:
  - 2.1
  - 2.2.1
  - 2.2.2
  - 2.2.3
  - 2.2.4
  - 2.2.5
  - 3.1.1 a 3.1.6
  - 4
  - Apendice
- A documentacao antiga fundia parte dessas secoes e deixava outras incompletas.

## 3. Problemas de OCR e formulas quebradas

Foram encontrados varios padroes repetidos:

- `W*tm`, `E*t`, `T*{tm}` em vez de `W_tm`, `E_t`, `T_t`
- matrizes com barras invertidas ausentes, por exemplo:
  - `\begin{bmatrix} ... \ ... \end{bmatrix}`
- determinantes 3D escritos em linha unica, sem formataçao valida de LaTeX
- produtos vetoriais e derivadas normais com caracteres faltando
- simbolos de media sobre o triangulo/tetraedro ora como barra, ora como texto, ora perdidos

## 4. Problemas de figuras e recursos visuais

- `docs/2.2.3_Determinação do número de onda.md`
  - trata a Figura 11 como placeholder: "coloque a imagem na pasta `/images`"
- `docs/2.2.4_Caracteristicas_de_Dispersao_de_Guias_de_Onda.md`
  - trata a Figura 13 como placeholder
- `docs/3_Problemas_Tridimensionais.md`
  - trata a Figura 14 como placeholder
- os links `images/fig11_waveguide.png`, `images/fig13_waveguide.png` e `images/fig14_tetrahedron.png` nao existem no repositorio
- em `docs/2_Problemas_Bidimensionais.md`, as legendas das Figuras 3, 4, 6 e 8 estavam como `Figura X. .`

## 5. Problemas de traducao/edicao

- `docs/0_Abstract.md`
  - contem marcador residual `:contentReference[...]`, que nao pertence ao texto final
- `docs/2_Problemas_Bidimensionais.md`
  - cita Figuras 5, 7 e 9 no corpo do texto sem integracao documental correspondente
- `docs/2.2.md`
  - a parte matematica principal estava usavel, mas o texto introdutorio ja misturava traducao do artigo com interpretacao do projeto

## 6. Impacto tecnico dessas inconsistencias

Esses problemas nao sao apenas cosmeticos. Eles prejudicam:

- a leitura sequencial do artigo
- a rastreabilidade entre equacao, secao e modulo do codigo
- a confianca na notacao usada para `St`, `Tt`, `Sz`, `Tz`, `C`, `Sel`, `Tel`
- a comparacao precisa entre paper e implementacao

## 7. O que foi feito nesta revisao

- criacao de uma nova pasta autoritativa:
  - `docs/traducao_revisada/`
- separacao por bloco logico do artigo
- restauracao do encadeamento teorico principal
- registro explicito das equacoes globais principais de cada secao
- isolamento das inconsistencias dos arquivos antigos em vez de tentar "consertar" os rascunhos misturados sem contexto

## 8. O que ainda fica pendente para uma segunda passada teorica

- redigitalizar ou redesenhar as Figuras 11, 13 e 14 se quisermos uma traducao visual completa
- revisar linha a linha a OCR de todas as equacoes auxiliares menos centrais, sobretudo nas familias de coeficientes geometricos 3D
- se desejado, produzir uma versao ainda mais literal do artigo, pagina a pagina, alem desta versao revisada por secao
