# Traducao Revisada do NASA TP-3485

Esta pasta reorganiza a teoria do artigo *Finite Element Method for Eigenvalue Problems in Electromagnetics* (NASA TP-3485) em uma trilha legivel, rastreavel e mais alinhada com o codigo do repositorio.

## Fontes usadas na revisao

- PDF local: `docs/19950011772.pdf`
- copia oficial no NTRS da NASA:
  `https://ntrs.nasa.gov/citations/19950011772`
- documentacao tecnica dos modulos em `src/*/README.md`, usada apenas para checagem de coerencia matematica com a implementacao

## Criterios desta revisao

- manter a numeracao original das equacoes do artigo
- separar a traducao por bloco logico, e nao por fragmentos misturados
- apontar quando a documentacao antiga misturava secoes, pulava numeracao ou deformava formulas
- deixar claro o que e teoria do artigo e o que ja e interpretacao voltada ao codigo

## Estrutura

- `01_abstract_e_introducao.md`
  - resumo do objetivo do artigo e mapa das secoes
- `02_secao_2_1_guias_homogeneos_escalar.md`
  - formulacao escalar 2D, Eq. (1) a Eq. (43)
- `03_secao_2_2_1_guias_homogeneos_vetorial.md`
  - formulacao vetorial transversal 2D com elementos de aresta
- `04_secao_2_2_2_guias_inhomogeneos_tres_componentes.md`
  - sistema misto no cutoff, base de `HELMVEC1`
- `05_secao_2_2_3_k0_para_beta_dado.md`
  - problema generalizado para `k0` com `beta` dado, base de `HELMVEC2`
- `06_secao_2_2_4_beta_para_k0_dado.md`
  - rearranjo para `beta` com `k0` dado, base de `HELMVEC3`
- `07_secao_3_1_problemas_tridimensionais.md`
  - formulacao vetorial 3D em tetraedros, Eq. (143) a Eq. (182)
- `08_conclusoes_apendice_referencias.md`
  - conclusoes, programas do apendice e referencias do artigo
- `99_inconsistencias_da_documentacao_legada.md`
  - auditoria do que estava inconsistente nos arquivos antigos de `docs/`

## Mapa rapido do artigo

| Parte | Tema | Equacao global principal | Programa do artigo |
| --- | --- | --- | --- |
| 2.1 | Guia homogeneo, formulacao escalar | Eq. (34) / Eq. (43) no repositorio | `HELM10` |
| 2.2.1 | Guia homogeneo, formulacao vetorial transversal | Eq. (65) | `HELMVEC` |
| 2.2.2 | Guia inhomogeneo, tres componentes no cutoff | Eq. (92) | `HELMVEC1` |
| 2.2.3 | `k0` para `beta` dado | Eq. (119) | `HELMVEC2` |
| 2.2.4 | `beta` para `k0` dado | Eq. (136) | `HELMVEC3` |
| 3.1 | Cavidades 3D com elementos de aresta | Eq. (178) | `FEM3D0` / `FEM3D1` |

## Observacoes importantes

- As figuras 11, 13 e 14 nao estao materializadas corretamente na documentacao legada; a trilha revisada registra isso explicitamente.
- Onde a OCR/traducao antiga estava ruidosa, a prioridade nesta pasta foi:
  1. preservar o papel matematico de cada equacao;
  2. restaurar o encadeamento da derivacao;
  3. alinhar a notacao com o uso do repositorio.
- Esta pasta nao substitui o artigo original. Ela serve como guia tecnico revisado para leitura e implementacao.
