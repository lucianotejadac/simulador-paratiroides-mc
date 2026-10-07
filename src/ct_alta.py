"""CT en la resolución original del tomógrafo (0.98 mm en el plano, 3.27 mm entre cortes), para mostrar y exportar.

El fantoma de 2 mm sigue siendo el mapa de atenuación del Monte Carlo; este CT es solo para la vista:
mismo recorte del cuello, misma máscara corporal (sin camilla), mismas coordenadas LPS.

Salida:
  salida/ct_alta.npz    hu (int16, z,y,x), origen LPS del vóxel (0,0,0), espaciado (dz, dy, dx)
  docs/datos/ct_hd.bin  uint16 (nz_web, 512, 512): CT remuestreado a los cortes de la grilla SPECT de la
                        página, a 4 sub-píxeles por píxel SPECT (0.825 mm), en HU + 1024
  docs/datos/ct_hd.json forma y codificación
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy import ndimage

import fantoma

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUB = 4          # sub-píxeles por píxel de la grilla SPECT en la página


def generar():
    meta = json.load(open(os.path.join(RAIZ, "salida", "fantoma.json"), encoding="utf-8"))
    hu, z, ps, origen, _ = fantoma.leer_ct(meta["ct"])
    sub, (a, b, y0, y1, x0, x1), _ = fantoma.recortar(hu, z, ps, meta.get("abajo_mm", 120.0), meta.get("arriba_mm", 70.0))
    assert [a, b, y0, y1, x0, x1] == meta["recorte_indices"], "el recorte no coincide con el del fantoma"
    org = [origen[0] + x0 * ps[1], origen[1] + y0 * ps[0], float(z[a])]
    dz = float(np.median(np.diff(z[a:b])))
    np.savez_compressed(os.path.join(RAIZ, "salida", "ct_alta.npz"), hu=np.round(sub).astype(np.int16),
                        origen=np.array(org), espaciado=np.array([dz, ps[0], ps[1]]))
    return np.round(sub).astype(np.int16), org, (dz, ps[0], ps[1]), meta


def grilla():
    """Matriz y píxel de la adquisición SPECT (los de las proyecciones del caso normal)."""
    d = np.load(os.path.join(RAIZ, "salida", "casos", "normal", "proyecciones_precoz.npz"))
    return int(d["matriz"]), float(d["pixel_mm"])


def para_web(sub, org, esp, meta):
    """Muestrea el CT original en los cortes de la grilla SPECT de la página (misma correspondencia que exportar_web)."""
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    iso = float(f["iso"])
    forma_f = f["hu"].shape
    matriz, pix = grilla()
    esc = iso / pix
    forma_z = [int(round(s * esc)) for s in forma_f]
    desplaz = [(matriz - s) // 2 for s in forma_z]
    idx = json.load(open(os.path.join(RAIZ, "docs", "datos", "indice.json"), encoding="utf-8"))
    z0, nz = idx["z0"], idx["nz"]
    o_f = meta["origen_mm"]
    # coordenadas de cada sub-píxel: grilla g -> fantoma p = (g - desplaz)/esc -> LPS = o_f + p*iso -> índice del CT original
    gz = z0 + np.arange(nz, dtype=np.float64)
    gu = (np.arange(matriz * SUB) + 0.5) / SUB - 0.5
    def a_lps(g, eje):
        return o_f[eje] + (g - desplaz[2 - eje]) / esc * iso
    zl, yl, xl = a_lps(gz, 2), a_lps(gu, 1), a_lps(gu, 0)
    iz = (zl - org[2]) / esp[0]
    iy = (yl - org[1]) / esp[1]
    ix = (xl - org[0]) / esp[2]
    Z, Y, X = np.meshgrid(iz, iy, ix, indexing="ij")
    vol = ndimage.map_coordinates(sub.astype(np.float32), [Z, Y, X], order=1, mode="constant", cval=-1000.0)
    cuerpo = (vol > -900).any(axis=0)
    ys, xs = np.nonzero(cuerpo)
    y0, y1, x0, x1 = max(0, ys.min() - 4), min(vol.shape[1], ys.max() + 5), max(0, xs.min() - 4), min(vol.shape[2], xs.max() + 5)
    out = np.clip(np.round(vol[:, y0:y1, x0:x1] + 1024.0), 0, 65535).astype(np.uint16)
    out.tofile(os.path.join(RAIZ, "docs", "datos", "ct_hd.bin"))
    json.dump({"nz": int(nz), "lado": matriz * SUB, "sub": SUB, "recorte_yx": [int(y0), int(x0), int(y1 - y0), int(x1 - x0)],
               "codificacion": "uint16, HU + 1024", "pixel_mm": pix / SUB,
               "fuente": "CT original del tomógrafo (0.98 mm en el plano, 3.27 mm entre cortes), interpolación lineal"},
              open(os.path.join(RAIZ, "docs", "datos", "ct_hd.json"), "w", encoding="utf-8"), ensure_ascii=False)
    print(f"web: ct_hd.bin {out.shape} ({out.nbytes / 1e6:.1f} MB)")


def main():
    sub, org, esp, meta = generar()
    print(f"CT alta: {sub.shape}, espaciado {esp}, origen {org}")
    para_web(sub, org, esp, meta)


if __name__ == "__main__":
    main()
