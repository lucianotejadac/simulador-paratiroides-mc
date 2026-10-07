"""Etiquetas de las regiones (TotalSegmentator + umbral) en los cortes de la página, para dibujar contornos.

Muestrea salida/regiones.npz (grilla del fantoma, 2 mm) por vecino más cercano en los cortes de la grilla
SPECT de la página, a 2 sub-píxeles por píxel SPECT (1.65 mm), con el mismo recorte que ct_hd.bin.
Salida: docs/datos/contornos.bin (uint8, nz × alto × ancho) y docs/datos/contornos.json.
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy import ndimage

from segmentar import NOMBRES

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUB = 2
FUENTE = {"tiroides": "TotalSegmentator", "tráquea": "TotalSegmentator (laringe por umbral)", "esófago": "TotalSegmentator",
          "vasos": "TotalSegmentator (carótidas, subclavias, tronco y venas braquiocefálicas)", "hueso": "TotalSegmentator (vértebras, clavículas) + umbral",
          "pulmón/vía aérea": "TotalSegmentator + umbral", "parótidas": "aproximación por elipsoide", "submaxilares": "aproximación por elipsoide"}


def main():
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    iso = float(f["iso"])
    reg = np.load(os.path.join(RAIZ, "salida", "regiones.npz"))["reg"]
    idx = json.load(open(os.path.join(RAIZ, "docs", "datos", "indice.json"), encoding="utf-8"))
    hd = json.load(open(os.path.join(RAIZ, "docs", "datos", "ct_hd.json"), encoding="utf-8"))
    matriz, pix = 128, 3.3
    esc = iso / pix
    forma_z = [int(round(s * esc)) for s in reg.shape]
    desplaz = [(matriz - s) // 2 for s in forma_z]
    z0, nz = idx["z0"], idx["nz"]
    gz = z0 + np.arange(nz, dtype=np.float64)
    gu = (np.arange(matriz * SUB) + 0.5) / SUB - 0.5
    pz, py, px = (gz - desplaz[0]) / esc, (gu - desplaz[1]) / esc, (gu - desplaz[2]) / esc
    Z, Y, X = np.meshgrid(pz, py, px, indexing="ij")
    et = ndimage.map_coordinates(reg, [Z, Y, X], order=0, mode="constant", cval=0).astype(np.uint8)
    # mismo recorte que el CT de alta resolución (en sus unidades / 2)
    r = hd["recorte_yx"]
    k = hd["sub"] // SUB
    y0, x0, h, w = r[0] // k, r[1] // k, -(-r[2] // k), -(-r[3] // k)
    out = np.ascontiguousarray(et[:, y0:y0 + h, x0:x0 + w])
    out.tofile(os.path.join(RAIZ, "docs", "datos", "contornos.bin"))
    json.dump({"nz": int(nz), "lado": matriz * SUB, "sub": SUB, "recorte_yx": [int(y0), int(x0), int(out.shape[1]), int(out.shape[2])],
               "etiquetas": {str(i): {"nombre": n, "fuente": FUENTE.get(n, "")} for i, n in enumerate(NOMBRES)}},
              open(os.path.join(RAIZ, "docs", "datos", "contornos.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print(f"contornos.bin {out.shape} ({out.nbytes / 1e6:.1f} MB); etiquetas presentes {sorted(int(v) for v in np.unique(out))}")


if __name__ == "__main__":
    main()
