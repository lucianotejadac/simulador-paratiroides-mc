"""Volúmenes para el editor de contornos de la página (grilla del fantoma, 2 mm isotrópicos).

docs/datos/editor/ct.u8.gz         CT codificado en 1 byte: código = (HU + 1000) / 10, de −1000 a 1550 HU
docs/datos/editor/etiquetas.u8.gz  regiones.npz tal como las usa la simulación
docs/datos/editor/editor.json      forma (z, y, x), tamaño de vóxel, origen LPS y nombres de las regiones
El navegador los descomprime con DecompressionStream. Las correcciones descargadas se aplican a esta misma
grilla con importar_correcciones.py.
"""
from __future__ import annotations

import gzip
import json
import os
import sys

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(RAIZ, "src"))


def main():
    from segmentar import NOMBRES   # en el renograma este módulo se llama regiones
    exportar(NOMBRES)


def exportar(nombres):
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    hu, iso, origen = f["hu"], float(f["iso"]), [float(v) for v in f["origen"]]
    reg = np.load(os.path.join(RAIZ, "salida", "regiones.npz"))["reg"].astype(np.uint8)
    d = os.path.join(RAIZ, "docs", "datos", "editor")
    os.makedirs(d, exist_ok=True)
    codigo = np.clip(np.round((hu.astype(np.float32) + 1000.0) / 10.0), 0, 255).astype(np.uint8)
    for nombre, arr in (("ct.u8.gz", codigo), ("etiquetas.u8.gz", reg)):
        with open(os.path.join(d, nombre), "wb") as fh:
            fh.write(gzip.compress(np.ascontiguousarray(arr).tobytes(), compresslevel=9, mtime=0))
    json.dump({"forma_zyx": list(reg.shape), "voxel_mm": iso, "origen_lps_mm": origen, "nombres": list(nombres),
               "codificacion_ct": "codigo = (HU + 1000) / 10"}, open(os.path.join(d, "editor.json"), "w", encoding="utf-8"), ensure_ascii=False)
    tam = {n: round(os.path.getsize(os.path.join(d, n)) / 1e6, 2) for n in ("ct.u8.gz", "etiquetas.u8.gz")}
    print(f"editor: {reg.shape} a {iso} mm; {tam} MB")


if __name__ == "__main__":
    main()
