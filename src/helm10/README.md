# `helm10` - Formulacao escalar 2D para `kc` (Sec. 2.1)

Este modulo implementa a etapa escalar do artigo para guias homogeneos,
usando elementos triangulares nodais `P1` e solve generalizado simetrico.

## 1) Objetivo do modulo

Resolver:

- problema TE escalar (Neumann) para cutoff `kc`,
- problema TM escalar (Dirichlet) para cutoff `kc`,
- comparacao modal com referencias analiticas de geometria padrao.

Executaveis:
- `main_helm10_rect.cpp`
- `main_helm10_circle.cpp`
- `main_helm10_coax.cpp`
- `main_field_reconstruction_example.cpp`

Utilitarios compartilhados:
- `scalar_mode_post.hpp`
- `field_reconstruction.hpp`
- `field_reconstruction.cpp`

## 2) Modelo continuo e forma fraca

No corte transversal `Gamma`:

`nabla_t^2 phi + kc^2 phi = 0`

Forma fraca:

`int_Gamma (nabla_t T_s . nabla_t phi) dA = kc^2 int_Gamma (T_s phi) dA + int_dGamma T_s (dphi/dn) dl`

No codigo, isso vira o EVP:

`S u = lambda T u`, com `lambda = kc^2`.

## 3) Condicoes de contorno

- `ScalarBC::TE_Neumann`
  - condicao natural `dphi/dn = 0`,
  - mantem nos de contorno,
  - possui modo constante (`lambda ~ 0`) removido no pos-processamento.

- `ScalarBC::TM_Dirichlet`
  - condicao essencial `phi = 0` em PEC,
  - elimina nos de contorno do espaco discreto.

## 4) Discretizacao FEM (`P1`)

Por triangulo:

`phi_h = sum_{i=1..3} phi_i N_i`

com gradiente constante:

`grad N_i = [b_i, c_i] / (2A)`

Matrizes elementares:

`S_e(i,j) = int_e (1/mu_r) grad N_i . grad N_j dA`

`T_e(i,j) = int_e eps_r N_i N_j dA`

## 4.1) Expressoes fechadas do elemento triangular

Usando:

```text
N_i(x,y) = (a_i + b_i x + c_i y) / (2A)
grad N_i = [b_i, c_i] / (2A)
```

as expressoes implementadas no backend `closed-form` ficam:

```text
S_e(i,j) = (1/mu_r) * (b_i b_j + c_i c_j) / (4A)
```

```text
T_e(i,j) = eps_r * (A/12) * [2 1 1; 1 2 1; 1 1 2]_(i,j)
```

ou seja:

```text
T_e =
eps_r * A/12 *
[ 2 1 1
  1 2 1
  1 1 2 ]
```

Na forma global:

```text
S phi = kc^2 T phi
```

com:

- TE: condicao natural de Neumann;
- TM: condicao essencial de Dirichlet.

## 4.2) Trilha de rastreabilidade

Para conectar o artigo ao codigo, a trilha principal deste modulo e:

1. Eq. `(30)` e Eq. `(33)` -> formas locais closed-form em
   `src/explicit/tri2d_scalar_explicit.hpp`
2. entrada publica didatica da Eq. `(43)` ->
   `tp3485::build_eq43_helm10_system(...)` em `src/article/tp3485_systems.hpp`
3. Eq. `(43)` -> sistema global escalar montado em
   `src/core/helm10_scalar_system.cpp`
4. funcao de montagem principal subjacente ->
   `build_helm10_scalar_system(...)`

## 5) Solver e pos-processamento

Solver:
- `generalized_eigs_sym_vec` (`LAPACKE_dsygv`).

Pos-processamento:
- filtro de autovalores fisicos,
- extracao de autovetor nodal em `scalar_mode_post.hpp`,
- reconstrucao do gradiente do potencial pelas Eq. `(36)` e `(37)` em
  `field_reconstruction.cpp`,
- reconstrucao do campo eletrico TE pelas Eq. `(38)` e `(39)`:

`Ex = -dpsi/dy`, `Ey = dpsi/dx`

- reconstrucao do campo eletrico TM pelas Eq. `(40)` e `(41)`:

`Ex = -dpsi/dx`, `Ey = -dpsi/dy`

- e o CSV tambem registra a versao fisica escalada:

`Ex_com_ztm = -Ztm * dpsi/dx`, `Ey_com_ztm = -Ztm * dpsi/dy`

- suavizacao nodal por media ponderada de area,
- normalizacao para quiver.

Observacao importante para o ramo TM:

- o problema de cutoff fornece `kc`;
- o projeto continua calculando, quando a frequencia e conhecida:
  `k = omega * sqrt(mu * eps)`,
  `beta = sqrt(k^2 - kc^2)`,
  `Ztm = beta / (omega * eps)`;
- mas os campos `Ex` e `Ey` exportados no ramo TM usam apenas gradiente e
  sinal, sem multiplicacao por `Ztm`;
- quando `k^2 - kc^2 <= 0`, o modo e marcado como abaixo do corte e um aviso
  e emitido no console, mas o padrao transversal ainda e salvo de forma
  didatica, sem escalonamento por impedancia;
- a verificacao de triangulo degenerado tambem ocorre dentro de
  `field_reconstruction.cpp`.

Observacao importante sobre a frequencia:

- se `--freq-hz` for informado, o valor do usuario e respeitado;
- se `--freq-hz` for omitido, o `helm10` escolhe automaticamente uma
  frequencia a partir dos `kc` exportados;
- quando houver modos TM exportados, a politica automatica tenta usar o
  modo TM limitante como referencia e impor `Ztm = 1` para ele;
- se necessario, a frequencia e elevada um pouco para garantir que todos os
  modos exportados permanecam propagantes;
- se nao houver TM entre os modos exportados, o codigo cai para a regra
  generica de frequencia acima do maior `kc`.

## 6) Arquivos e responsabilidades

- `src/core/helm10_scalar_system.{hpp,cpp}`: montagem escalar global.
- `src/core/fem_scalar_helm10.cpp`: blocos elementares e assembleia.
- `src/core/mode_match_rect.hpp`: matching retangular por correlacao.
- `src/core/mode_match_circle.hpp`: matching circular por correlacao.
- `src/core/mode_match_coax.hpp`: matching coaxial por correlacao.
- `src/core/io_vtk_sv.hpp`: export VTK escalar+vetorial nodal.
- `src/helm10/field_reconstruction.{hpp,cpp}`: gradiente, `Ex`, `Ey`, escolha
  automatica de frequencia, `beta`, `Ztm` informativo e verificacoes fisicas.
- `docs/HELM10_CSV_Modos_Referencia.md`: glossario didatico das colunas do
  CSV de modos exportado pelos executaveis da familia `HELM10`.
- `docs/HELM10_CSV_Campos_Referencia.md`: glossario didatico das colunas do
  CSV de campos nodais exportado pelos executaveis da familia `HELM10`.

## 7) Uso

Retangular:

```bash
./build/helm10_rect
./build/helm10_rect 1.0 14 14 20
```

Assinatura preferencial de `helm10_rect`:

```bash
./build/helm10_rect [ar_m [nx [ny [nmodos]]]]
```

Compatibilidade legada ainda aceita:

```bash
./build/helm10_rect [nx ny [nmodos]]
```

Convenio do caso retangular:

- `ar_m`: largura do guia retangular em metros; padrao `1.0`;
- `b = ar/2`, mantendo a razao didatica `a_r / b_r = 2`;
- `nmodos`: quantidade de modos TE e quantidade de modos TM exportados;
- se `nmodos` nao for informado, o padrao publico agora e `20`.

Circular:

```bash
./build/helm10_circle
./build/helm10_circle 10 48 20
```

Coaxial:

```bash
./build/helm10_coax
./build/helm10_coax 10 48 20
```

Convenio dos casos circular e coaxial:

- `nr`, `nt`: discretizacao radial/angular da malha;
- `nmodos`: quantidade de modos TE e quantidade de modos TM exportados;
- se `nmodos` nao for informado, o padrao publico agora tambem e `20`.

Exemplo didatico isolado:

```bash
./build/helm10_field_reconstruction_example
```

## 7.1) Backends e depuracao

Flags disponiveis:

- `--backend closed-form`
- `--backend gauss`
- `--freq-hz <valor>`
- `--eps-r <valor>`
- `--mu-r <valor>`
- `--debug-local-blocks`
- `--debug-candidates`
- `--debug` ou `--debug-all`

Sem `--backend`, o padrao publico do executavel e `closed-form`.
Sem `--freq-hz`, a frequencia e escolhida automaticamente acima do maior
`kc` exportado.

Exemplos:

```bash
./build/helm10_rect 1.0 14 14 8 --backend closed-form --debug-local-blocks
./build/helm10_rect 0.8 12 12 6 --freq-hz 1e9 --backend closed-form
./build/helm10_circle 10 48 8 --freq-hz 1e9 --debug-candidates
./build/helm10_coax 10 48 8 --freq-hz 1e9 --backend closed-form --debug-local-blocks --debug-candidates
./build/helm10_field_reconstruction_example
```

Interpretacao:

- `--debug-local-blocks`: imprime o primeiro triangulo com os blocos locais
  ligados as Eq. `(30)` e `(33)` e sua contribuicao para a Eq. `(43)`;
- `--debug-candidates`: imprime os primeiros `kc` positivos antes do matching.

Aviso didatico importante:

- no ramo TM, `Ex` e `Ey` continuam sendo salvos como saida principal sem
  multiplicacao por `Ztm`;
- o CSV de campos salva apenas essa versao sem impedancia, para nao misturar
  duas normalizacoes no mesmo arquivo;
- `beta` e `Ztm` permanecem nos CSVs como metadados fisicos auxiliares.

## 7.2) Geracao de imagens a partir dos CSVs

O pos-processamento grafico didatico do `HELM10` fica em:

- `python3 scripts/helm10.py`

Esse script le `modes.csv` e `fields_<modo>.csv` e gera, para cada modo:

- diagrama de isolinhas do potencial escalar `psi`
- diagrama de setas do campo transversal `(Ex, Ey)`
- grafico-resumo de erro por modo, ordenado por `kc_fem` dentro de cada familia
- grafico-resumo de `rho`, tambem ordenado por `kc_fem` dentro de cada familia

Uso basico:

```bash
python3 scripts/helm10.py
python3 scripts/helm10.py --case rect
python3 scripts/helm10.py --case circle --case coax --show-mesh
```

Pastas geradas:

- `out/helm10/rect/img/`
- `out/helm10/circle/img/`
- `out/helm10/coax/img/`

Subpastas de cada caso:

- `img/isopotential/`
- `img/quiver/`

Arquivo-resumo de erro por caso:

- `img/helm10_<caso>_error_by_mode.png`

Observacao importante:

- os CSVs sao a fonte principal dos valores;
- quando o VTK correspondente existe, o script reaproveita a conectividade
  triangular da malha para evitar artefatos geometricos, especialmente no
  caso coaxial.

Guia dedicado:

- `docs/HELM10_Imagens_Referencia.md`

## 8) Saidas tipicas

Console:

- primeiros `kc`,
- tabela FEM x analitico (`kc_ana`, `kc_fem`, `kcAr_fem`, erro %, `rho`).

Arquivos VTK gerados pelos `main`:

- retangular: `out/helm10/rect/vtk/`
- circular: `out/helm10/circle/vtk/`
- coaxial: `out/helm10/coax/vtk/`

Arquivo de log por execucao:

- retangular: `out/helm10/rect/run.log`
- circular: `out/helm10/circle/run.log`
- coaxial: `out/helm10/coax/run.log`

Esse log espelha tudo o que o executavel imprime no terminal durante a
execucao normal: informacoes de malha, tabelas de modos, avisos fisicos,
arquivos salvos e linha final de `timing`.

Arquivo estruturado de tempos por execucao:

- retangular: `out/helm10/rect/run_timing.csv`
- circular: `out/helm10/circle/run_timing.csv`
- coaxial: `out/helm10/coax/run_timing.csv`

Esse CSV registra, em uma linha, a configuracao da rodada e os tempos medidos
para:

- montagem global TE e TM;
- solve do autoproblema TE e TM;
- totais de montagem, solve, pos-processamento e execucao.

Em cada pasta `vtk/`, os nomes carregam a familia e o label modal:

- retangular: `te_rect_m*_n*_rank*_sv.vtk`, `tm_rect_m*_n*_rank*_sv.vtk`
- circular: `te_circle_m*_p*_rank*_sv.vtk`, `tm_circle_m*_p*_rank*_sv.vtk`
- coaxial: `te_coax_m*_p*_rank*_sv.vtk`, `tm_coax_m*_p*_rank*_sv.vtk`

Os aliases historicos do primeiro modo continuam sendo escritos:

- retangular: `te10_rect_sv.vtk`, `tm11_rect_sv.vtk`
- circular: `te_circle_sv.vtk`, `tm_circle_sv.vtk`
- coaxial: `te_coax_sv.vtk`, `tm_coax_sv.vtk`

Arquivos CSV adicionais do caso retangular:

- `out/helm10/rect/csv/helm10_rect_modes.csv`
  - uma linha por modo TE/TM exportado;
  - inclui `m`, `n`, `kc_fem`, `kc_ana`, `kcAr_fem`, `kcAr_ana`, `beta`,
    `Ztm`, estado abaixo/acima do corte, erro, `rho` e o nome do CSV de
    campos correspondente ao modo;
  - `beta` e `Ztm` sao informativos; eles nao escalam `Ex` e `Ey` no ramo TM.
- `out/helm10/rect/run.log`
  - copia textual da execucao completa mostrada no terminal.
- `out/helm10/rect/run_timing.csv`
  - configuracao da rodada e tempos de montagem/solve em formato tabular.
- `out/helm10/rect/csv/helm10_rect_fields_<familia>_m<m>_n<n>_rank<rr>.csv`
  - um arquivo por modo exportado;
  - uma linha por no;
  - inclui apenas `node_id`, `x_m`, `y_m`, `psi`, `dpsi_dx`, `dpsi_dy`,
    `Ex`, `Ey`;
  - todos os metadados constantes do modo ficam concentrados em
    `helm10_rect_modes.csv`;
  - no caso TM, `Ex` e `Ey` seguem a mesma politica didatica do projeto:
    saida sempre sem multiplicacao por `Ztm`.

Arquivos CSV adicionais dos casos circular e coaxial:

- `out/helm10/circle/csv/helm10_circle_modes.csv`
- `out/helm10/circle/csv/helm10_circle_fields_<familia>_m<m>_p<p>_rank<rr>.csv`
- `out/helm10/circle/run.log`
- `out/helm10/circle/run_timing.csv`
- `out/helm10/coax/csv/helm10_coax_modes.csv`
- `out/helm10/coax/csv/helm10_coax_fields_<familia>_m<m>_p<p>_rank<rr>.csv`
- `out/helm10/coax/run.log`
- `out/helm10/coax/run_timing.csv`

Todos seguem a mesma leitura didatica:

- o CSV de modos concentra os metadados do modo;
- o CSV de campos concentra apenas valores nodais;
- potencial longitudinal (`Hz` para TE, `Ez` para TM),
  gradiente transversal e `Ex`, `Ey`;
- indicacao clara, no CSV de modos, de `k`, `beta`, `Ztm` e estado acima ou
  abaixo do corte;
- aviso explicito de que o ramo TM e exportado sem multiplicacao por `Ztm`.

Arquivos de imagem gerados pelo script Python:

- `out/helm10/rect/img/isopotential/*.png`
- `out/helm10/rect/img/quiver/*.png`
- `out/helm10/rect/img/helm10_rect_error_by_mode.png`
- `out/helm10/rect/img/helm10_rect_rho_by_mode.png`
- `out/helm10/circle/img/isopotential/*.png`
- `out/helm10/circle/img/quiver/*.png`
- `out/helm10/circle/img/helm10_circle_error_by_mode.png`
- `out/helm10/circle/img/helm10_circle_rho_by_mode.png`
- `out/helm10/coax/img/isopotential/*.png`
- `out/helm10/coax/img/quiver/*.png`
- `out/helm10/coax/img/helm10_coax_error_by_mode.png`
- `out/helm10/coax/img/helm10_coax_rho_by_mode.png`

## 9) Relacao com o artigo

Corresponde ao fluxo da Sec. 2.1:
- formulacao escalar,
- discretizacao triangular nodal,
- validacao em guias retangular/circular/coaxial.

Extensoes de engenharia no repositorio:
- matching modal automatico por correlacao em massa,
- exportacao VTK pronta para pipeline de figuras.

Para a formula usada no calculo de `rho_abs` e a justificativa didatica de por
que esse matching funciona, ver:

- [docs/HELM10_CSV_Modos_Referencia.md](/home/sperotto/tp3485-fem-eigen-em/docs/HELM10_CSV_Modos_Referencia.md)
