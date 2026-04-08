# FEM3D - Referência das Imagens

Este arquivo documenta o script:

- `scripts/fem3d.py`

Ele lê as saídas já exportadas por:

- `fem3d0_air`, `fem3d0_half`, `fem3d0_cyl`, `fem3d0_sphere`
- `fem3d1_air`, `fem3d1_half`, `fem3d1_cyl`, `fem3d1_sphere`

e gera imagens-resumo e projeções ortogonais dos campos 3D exportados por
modo casado.

## 1) Entradas lidas

Por caso, o script usa:

- `csv/<solver>_<caso>_modes.csv`
- `csv/<solver>_<caso>_modeXX_<modo>_E_fields.csv`
- `vtk/<solver>_<caso>_modeXX_<modo>_E.vtk`
- `run_timing.csv`

Os arquivos de campo vêm da exportação espacial dos modos casados do `FEM3D0`
e do `FEM3D1`.

## 2) Imagens geradas

Em cada pasta:

```text
out/<solver>/<caso>/img/
```

o script grava:

- `fem3d?_case_k0_by_mode.png`
  - comparação `k0_fem x k0_analytic x ref_paper`
- `fem3d?_case_error_by_mode.png`
  - erro relativo absoluto por modo
- `img/magnitude/*_magnitude_proj.png`
  - três projeções ortogonais `XY`, `XZ` e `YZ` da magnitude `|E|`
- `img/quiver/*_quiver_proj.png`
  - três projeções ortogonais `XY`, `XZ` e `YZ` das componentes vetoriais
- `img/3d_scatter/*_scatter3d.png`
  - visualização tridimensional dos centroides, colorida por `|E|`, com a casca geométrica do domínio
- `img/3d_quiver/*_quiver3d.png`
  - visualização tridimensional com setas subamostradas de `E`, também sobre a casca geométrica

## 3) Convenção geométrica das projeções

O script não faz corte volumétrico nem renderização 3D completa. Em vez disso,
ele usa os centroides dos tetraedros e produz:

- projeção `XY` com vetores `(Ex, Ey)`
- projeção `XZ` com vetores `(Ex, Ez)`
- projeção `YZ` com vetores `(Ey, Ez)`

Isso é propositalmente simples e didático:

- mantém o vínculo direto com os CSVs exportados pelo solver
- permite comparar rapidamente formas modais entre casos e malhas
- evita esconder a informação em uma visualização 3D difícil de ler

O contorno projetado do domínio é obtido a partir dos nós de fronteira do VTK
tetraédrico correspondente ao modo.

Além das projeções 2D, o script agora gera duas vistas 3D rápidas:

- `scatter3d`
  - centroides coloridos por `Emag`
- `quiver3d`
  - vetores `E` em uma subamostragem dos centroides

Nas duas vistas 3D, a geometria é desenhada como uma casca translúcida obtida
das faces de fronteira do VTK tetraédrico exportado pelo solver.

Essas figuras não substituem o uso do `VTK` em ParaView, mas ajudam bastante
na inspeção rápida dentro do próprio fluxo do repositório.

No caso `half`, o script também desenha a interface dielétrica em `z = lz/2`
nas projeções em que ela aparece:

- `XZ`
- `YZ`

e também como um plano semitransparente nas figuras 3D.

## 4) Uso

Exemplos:

```bash
python3 scripts/fem3d.py
python3 scripts/fem3d.py --solver fem3d0 --case air
python3 scripts/fem3d.py --solver fem3d1 --case half --dpi 120
python3 scripts/fem3d.py --solver fem3d0 --solver fem3d1 --case cyl --max-arrows 250
```

## 5) Leitura conjunta recomendada

Para entender de onde vêm os dados que alimentam essas figuras:

- [FEM3D_CSV_Referencia.md](FEM3D_CSV_Referencia.md)
- [Artefatos_Espectrais_CSV_Referencia.md](Artefatos_Espectrais_CSV_Referencia.md)
- [src/fem3d/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d/README.md)
- [src/fem3d0/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d0/README.md)
- [src/fem3d1/README.md](/home/sperotto/tp3485-fem-eigen-em/src/fem3d1/README.md)
