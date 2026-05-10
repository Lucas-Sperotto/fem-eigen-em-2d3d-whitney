#!/usr/bin/env python3
"""
Copia imagens curadas de out/ para docs/results/img/ e gera manifest.json.

Estrutura de saida:
    docs/results/img/
        fem/<familia>/<geometria>/*.png   <- FEM pipeline padrao
        efgm/<familia>/<geometria>/*.png  <- EFGMI (relatorio comparativo)
        artigo/*.png                      <- figuras do artigo (de docs/figs/)

    docs/results/img/manifest.json       <- mapeamento caso -> listas de imagens
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out"
DEST = ROOT / "docs" / "results" / "img"
EFGMI_BASE = OUT / "fem_efgmi_mode_report_base" / "base"


# ---------------------------------------------------------------------------
# Mapeamento: chave_caso -> (src_dir_fem, src_dir_efgm, dest_subdir)
# ---------------------------------------------------------------------------
CASE_DIRS: list[dict] = [
    # 2D Scalar (Sec. 2.1)
    {"key": "helm10_rect",    "fem": OUT / "helm10" / "rect",    "efgm": EFGMI_BASE / "efgm" / "helm10" / "rect",    "sub": "helm10/rect"},
    {"key": "helm10_circle",  "fem": OUT / "helm10" / "circle",  "efgm": EFGMI_BASE / "efgm" / "helm10" / "circle",  "sub": "helm10/circle"},
    {"key": "helm10_coax",    "fem": OUT / "helm10" / "coax",    "efgm": EFGMI_BASE / "efgm" / "helm10" / "coax",    "sub": "helm10/coax"},
    # 2D Edge (Sec. 2.2.1)
    {"key": "helmvec_rect",   "fem": OUT / "helmvec" / "rect",   "efgm": None,                                         "sub": "helmvec/rect"},
    {"key": "helmvec_circle", "fem": OUT / "helmvec" / "circle", "efgm": None,                                         "sub": "helmvec/circle"},
    # 2D Mixed (Sec. 2.2.2)
    {"key": "helmvec1_rect",   "fem": OUT / "helmvec1" / "rect",   "efgm": EFGMI_BASE / "efgm" / "helmvec1" / "rect",   "sub": "helmvec1/rect"},
    {"key": "helmvec1_circle", "fem": OUT / "helmvec1" / "circle", "efgm": EFGMI_BASE / "efgm" / "helmvec1" / "circle", "sub": "helmvec1/circle"},
    # 2D Coupled k0(beta) (Sec. 2.2.3)
    {"key": "helmvec2_rect",  "fem": OUT / "helmvec2" / "rect",  "efgm": EFGMI_BASE / "efgm" / "helmvec2" / "rect",  "sub": "helmvec2/rect"},
    # 2D Coupled beta(k0) (Sec. 2.2.4)
    {"key": "helmvec3_fig12", "fem": OUT / "helmvec3" / "fig12_rect", "efgm": EFGMI_BASE / "efgm" / "helmvec3" / "fig12_rect", "sub": "helmvec3/fig12_rect"},
    {"key": "helmvec3_fig13", "fem": OUT / "helmvec3" / "fig13_rect", "efgm": EFGMI_BASE / "efgm" / "helmvec3" / "fig13_rect", "sub": "helmvec3/fig13_rect"},
    # 3D Cavities (Sec. 3.1)
    {"key": "fem3d0_air",    "fem": OUT / "fem3d0" / "air",    "efgm": None, "sub": "fem3d0/air"},
    {"key": "fem3d0_half",   "fem": OUT / "fem3d0" / "half",   "efgm": None, "sub": "fem3d0/half"},
    {"key": "fem3d0_cyl",    "fem": OUT / "fem3d0" / "cyl",    "efgm": None, "sub": "fem3d0/cyl"},
    {"key": "fem3d0_sphere", "fem": OUT / "fem3d0" / "sphere", "efgm": None, "sub": "fem3d0/sphere"},
    {"key": "fem3d1_air",    "fem": OUT / "fem3d1" / "air",    "efgm": None, "sub": "fem3d1/air"},
    {"key": "fem3d1_half",   "fem": OUT / "fem3d1" / "half",   "efgm": None, "sub": "fem3d1/half"},
    {"key": "fem3d1_cyl",    "fem": OUT / "fem3d1" / "cyl",    "efgm": None, "sub": "fem3d1/cyl"},
    {"key": "fem3d1_sphere", "fem": OUT / "fem3d1" / "sphere", "efgm": None, "sub": "fem3d1/sphere"},
]


def _copy_img_tree(src_root: Path, dst_root: Path) -> list[str]:
    """Copia recursivamente todos os PNGs de src_root para dst_root.
    Retorna lista de caminhos relativos a DEST."""
    dst_root.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    if not src_root.exists():
        return copied
    img_root = src_root / "img"
    if not img_root.exists():
        return copied
    for png in sorted(img_root.rglob("*.png")):
        rel_to_img = png.relative_to(img_root)
        dst = dst_root / rel_to_img
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(png, dst)
        copied.append(str(dst.relative_to(DEST)))
    return copied


def _copy_article_figs() -> list[str]:
    """Copia as figuras originais do artigo para docs/results/img/artigo/."""
    src = ROOT / "docs" / "figs"
    dst = DEST / "artigo"
    dst.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    if not src.exists():
        return copied
    for png in sorted(src.glob("*.png")):
        shutil.copy2(png, dst / png.name)
        copied.append(str((dst / png.name).relative_to(DEST)))
    return copied


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, dict] = {}

    # Figuras do artigo
    artigo_imgs = _copy_article_figs()
    manifest["artigo"] = {"fem": [], "efgm": [], "artigo": artigo_imgs}
    print(f"[artigo] {len(artigo_imgs)} figuras copiadas")

    # Casos FEM + EFGMI
    for spec in CASE_DIRS:
        key = spec["key"]
        fem_src: Path = spec["fem"]
        efgm_src: Path | None = spec["efgm"]
        sub: str = spec["sub"]

        fem_dst = DEST / "fem" / sub
        efgm_dst = DEST / "efgm" / sub

        fem_imgs = _copy_img_tree(fem_src, fem_dst)
        efgm_imgs = _copy_img_tree(efgm_src, efgm_dst) if efgm_src else []

        manifest[key] = {"fem": fem_imgs, "efgm": efgm_imgs, "artigo": []}
        print(f"[{key}] fem={len(fem_imgs)} efgm={len(efgm_imgs)}")

    manifest_path = DEST / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nManifest salvo: {manifest_path.relative_to(ROOT)}")
    total_fem = sum(len(v["fem"]) for v in manifest.values())
    total_efgm = sum(len(v["efgm"]) for v in manifest.values())
    print(f"Total: {total_fem} FEM + {total_efgm} EFGM + {len(artigo_imgs)} artigo")


if __name__ == "__main__":
    main()
