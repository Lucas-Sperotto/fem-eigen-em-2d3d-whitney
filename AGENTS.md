# AGENTS.md

## Canonical Flow

- O fluxo canonico deste repositorio e `scripts/build_and_run_all.sh`.
- Use esse script como ponto de entrada recomendado para build, execucao, validacao e geracao de artefatos auxiliares.
- `scripts/build_and_run.sh` deve ser tratado como helper legado/local, nao como fluxo publico recomendado.

## CLI Policy

- Sempre prefira CLI nomeada quando o binario expuser aliases nomeados no `--help`.
- Posicionais podem existir por compatibilidade, mas devem ser tratados como legado/deprecated na documentacao e em exemplos novos.

## Output Layout

- Nao reintroduza `out/2d` como layout canonico.
- O layout canonico de executaveis e por familia/caso, por exemplo `out/helm10/...`, `out/helmvec/...`, `out/helmvec1/...`, `out/helmvec2/...`, `out/helmvec3/...`, `out/fem3d0/...`, `out/fem3d1/...`.
- Diferencie sempre:
  - saida de executaveis: `run.log`, `run_timing.csv`, `csv/`, `vtk/`, `linop/` sob `out/<familia>/<caso>/`
  - saida de validacao: artefatos de `scripts/validate_*.py` e relatorios agregados, tipicamente em `out/validation/`
  - saida de geracao de imagens: PNGs e resumos produzidos por scripts de pos-processamento, tipicamente em `out/<familia>/<caso>/img/` ou `out/img_all/`

## Planning

- Antes de mudancas amplas, monte e compartilhe um plano curto.
- Considere "mudanca ampla" qualquer alteracao que atravesse multiplas familias, modifique layout de saida, mude convencoes publicas de CLI ou mexa em varios scripts/READMEs de uma vez.

## Validation

- Sempre valide build apos mudancas de codigo, script ou documentacao operacional relevante.
- Quando aplicavel, rode pelo menos uma execucao representativa minima do fluxo afetado.
- Prefira uma chamada canonica nomeada nos exemplos de validacao.
- Ao documentar artefatos, diferencie o que vem do executavel bruto do que vem de pos-processamento.
