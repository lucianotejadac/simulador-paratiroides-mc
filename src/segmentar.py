"""Paso 2: regiones del fantoma para asignar captación de sestamibi.

Etiquetas (salida/regiones.npz, uint8):
    0 aire, 1 tejido blando, 2 hueso, 3 pulmón/vía aérea, 4 tiroides, 5 parótidas, 6 submaxilares, 7 tráquea
La tiroides se segmenta por densidad (es hiperdensa en CT sin contraste por el yodo: 70–200 HU) a ambos
lados de la tráquea, entre 3 y 9 cm bajo el mínimo del cuello. Las glándulas salivales se aproximan con
elipsoides a partir de puntos de referencia (rama mandibular y tráquea), porque en CT sin contraste no se
separan bien del tejido vecino; importan solo como captación de fondo.
"""
from __future__ import annotations

import json
import os

import numpy as np
from scipy import ndimage

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AIRE, BLANDO, HUESO, PULMON, TIROIDES, PAROTIDA, SUBMAXILAR, TRAQUEA = range(8)
NOMBRES = ["aire", "tejido blando", "hueso", "pulmón/vía aérea", "tiroides", "parótidas", "submaxilares", "tráquea"]


def traquea_por_corte(hu, cuerpo):
    """Centro (y, x) del lumen traqueal por corte: componente de aire dentro del cuerpo, cerca de la línea media, anterior."""
    nz, ny, nx = hu.shape
    centros = {}
    for k in range(nz):
        aire = (hu[k] < -500) & cuerpo[k]
        if not cuerpo[k].any():
            continue
        xm = np.nonzero(cuerpo[k].any(axis=0))[0].mean()     # línea media del cuerpo en este corte
        et, n = ndimage.label(aire)
        mejor = None
        for i in range(1, n + 1):
            m = et == i
            v = m.sum()
            if v < 20 or v > 1500:
                continue
            ys, xs = np.nonzero(m)
            cy, cx = ys.mean(), xs.mean()
            if abs(cx - xm) > 25 / 2.0:            # a menos de 25 mm de la línea media (vóxel 2 mm)
                continue
            puntaje = v - 3 * abs(cx - xm)
            if mejor is None or puntaje > mejor[0]:
                mejor = (puntaje, cy, cx, v)
        if mejor:
            centros[k] = (mejor[1], mejor[2])
    return centros


def elipsoide(forma, centro_zyx, radios_zyx):
    z, y, x = np.ogrid[:forma[0], :forma[1], :forma[2]]
    return (((z - centro_zyx[0]) / radios_zyx[0]) ** 2 + ((y - centro_zyx[1]) / radios_zyx[1]) ** 2
            + ((x - centro_zyx[2]) / radios_zyx[2]) ** 2) <= 1.0


def segmentar(hu: np.ndarray, iso: float, z_cuello_idx: int):
    nz, ny, nx = hu.shape
    cuerpo = np.stack([ndimage.binary_fill_holes(c) for c in (hu > -300)])   # corte a corte: la tráquea y los pulmones son huecos en 2D, no en 3D
    reg = np.full(hu.shape, AIRE, dtype=np.uint8)
    reg[cuerpo] = BLANDO
    reg[cuerpo & (hu < -300)] = PULMON
    reg[cuerpo & (hu > 200)] = HUESO
    centros = traquea_por_corte(hu, cuerpo)
    for k, (cy, cx) in centros.items():
        aire = (hu[k] < -500) & cuerpo[k]
        et, _ = ndimage.label(aire)
        reg[k][et == et[int(round(cy)), int(round(cx))]] = TRAQUEA
    # tiroides: 3 a 9 cm bajo el mínimo del cuello, densa, a ≤ 32 mm del centro traqueal, no detrás de él
    k0, k1 = z_cuello_idx - int(90 / iso), z_cuello_idx - int(54 / iso)   # por debajo del cartílago tiroides, cuyas láminas calcificadas también miden 70-200 HU
    cand = np.zeros(hu.shape, dtype=bool)
    for k in range(max(0, k0), min(nz, k1)):
        if k not in centros:
            continue
        cy, cx = centros[k]
        yy, xx = np.ogrid[:ny, :nx]
        cerca = ((yy - cy) ** 2 + (xx - cx) ** 2) <= (32 / iso) ** 2
        adelante = yy <= cy + 12 / iso            # anterior o lateral a la tráquea (y crece hacia posterior)
        cand[k] = cerca & adelante & (hu[k] >= 70) & (hu[k] <= 200) & (reg[k] == BLANDO)
    cand = ndimage.binary_opening(cand, structure=np.ones((1, 2, 2)))
    et, n = ndimage.label(cand)
    if n:
        tam = ndimage.sum(cand, et, range(1, n + 1))
        # lóbulos: los componentes de más de 0.8 mL (100 vóxeles de 8 mm³)
        for i, v in enumerate(tam, start=1):
            if v >= 100:
                reg[et == i] = TIROIDES
    reg[ndimage.binary_closing(reg == TIROIDES, structure=np.ones((3, 3, 3))) & (reg == BLANDO)] = TIROIDES
    # glándulas salivales: elipsoides referidos a la tráquea en el corte del cuello (mandíbula ~ 1-4 cm sobre el mínimo)
    if centros:
        ks = sorted(centros)
        kref = min(ks, key=lambda k: abs(k - z_cuello_idx))
        cy, cx = centros[kref]
        z_par = z_cuello_idx + int(25 / iso)
        for lado in (-1, 1):
            par = elipsoide(hu.shape, (z_par, cy + 5 / iso, cx + lado * 52 / iso), (22 / iso, 12 / iso, 10 / iso))
            reg[par & (reg == BLANDO)] = PAROTIDA
            sub = elipsoide(hu.shape, (z_cuello_idx + int(5 / iso), cy - 18 / iso, cx + lado * 28 / iso), (12 / iso, 10 / iso, 12 / iso))
            reg[sub & (reg == BLANDO)] = SUBMAXILAR
    return reg, centros


def montaje(hu, reg, ruta, k_tiroides):
    from PIL import Image

    def ventana(a, c=40, w=400):
        return np.clip((a - (c - w / 2)) / w * 255, 0, 255).astype(np.uint8)

    colores = {TIROIDES: (255, 80, 80), PAROTIDA: (80, 200, 255), SUBMAXILAR: (120, 255, 120), TRAQUEA: (255, 255, 0), PULMON: (60, 60, 160), HUESO: (200, 200, 200)}

    def pintar(h, r):
        g = ventana(h)
        rgb = np.stack([g, g, g], -1).astype(np.float32)
        for et, col in colores.items():
            m = r == et
            rgb[m] = 0.55 * rgb[m] + 0.45 * np.array(col)
        return Image.fromarray(rgb.astype(np.uint8))

    nz, ny, nx = hu.shape
    cor = pintar(hu[:, ny // 2, :][::-1], reg[:, ny // 2, :][::-1])
    ks = [max(0, k_tiroides - 8), k_tiroides - 4, k_tiroides, k_tiroides + 4, k_tiroides + 8, min(nz - 1, k_tiroides + 30)]
    axiales = [pintar(hu[k], reg[k]) for k in ks]
    W = max(cor.width, len(axiales) * (nx + 4))
    lienzo = Image.new("RGB", (W, cor.height + ny + 10), (0, 0, 0))
    lienzo.paste(cor, (0, 0))
    for i, im in enumerate(axiales):
        lienzo.paste(im, (i * (nx + 4), cor.height + 10))
    lienzo.save(ruta)
    return ks


def main():
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    meta = json.load(open(os.path.join(RAIZ, "salida", "fantoma.json"), encoding="utf-8"))
    hu, iso = f["hu"].astype(np.float32), float(f["iso"])
    i_a = meta["recorte_indices"][0]
    dz_ct = meta["ct_espaciado_mm"][0]
    # índice del corte del cuello dentro del fantoma remuestreado
    z_cuello_idx = int(round((meta["z_cuello_mm"] - meta["origen_mm"][2]) / iso))
    reg, centros = segmentar(hu, iso, z_cuello_idx)
    vol_ml = {NOMBRES[i]: round(float((reg == i).sum()) * iso ** 3 / 1000.0, 1) for i in range(8)}
    zt = np.nonzero((reg == TIROIDES).any(axis=(1, 2)))[0]
    k_t = int(zt.mean()) if len(zt) else z_cuello_idx - int(60 / iso)
    ks = montaje(hu, reg, os.path.join(RAIZ, "salida", "regiones.png"), k_t)
    np.savez_compressed(os.path.join(RAIZ, "salida", "regiones.npz"), reg=reg)
    json.dump({"volumen_ml": vol_ml, "nombres": NOMBRES, "z_cuello_idx": z_cuello_idx, "cortes_montaje": [int(k) for k in ks],
               "tiroides_cortes": [int(zt.min()), int(zt.max())] if len(zt) else None,
               "traquea_cortes": [min(centros), max(centros)] if centros else None},
              open(os.path.join(RAIZ, "salida", "regiones.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print("volúmenes (mL):", vol_ml)
    print("tiroides en cortes", zt.min() if len(zt) else None, "-", zt.max() if len(zt) else None, "| cuello idx", z_cuello_idx, "| tráquea en", min(centros) if centros else None, "-", max(centros) if centros else None)


if __name__ == "__main__":
    main()
