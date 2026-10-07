"""Paso 3: mapas de actividad de sestamibi (kBq/mL) por caso, fase precoz (15 min) y tardía (2 h).

Concentraciones de fondo tomadas de la biodistribución clínica del Tc-99m MIBI (orden de magnitud,
740 MBq inyectados, cuello): la tiroides lava ~60 % entre fases y el adenoma retiene ~80 %.
El adenoma se coloca por referencia anatómica a partir de la segmentación: polo inferior de un
lóbulo (posterior), retroesofágico o mediastínico (borde inferior del fantoma). La "verdad" de
cada caso (centro, diámetro, relación de captación) se guarda en verdad.json.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from scipy import ndimage

from segmentar import AIRE, BLANDO, HUESO, PULMON, TIROIDES, PAROTIDA, SUBMAXILAR, TRAQUEA, ESOFAGO, VASOS

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# kBq/mL por región: (precoz, tardía)
CONCENTRACION = {AIRE: (0.0, 0.0), BLANDO: (4.0, 2.5), HUESO: (3.0, 2.0), PULMON: (1.0, 0.6), TIROIDES: (40.0, 16.0),
                 PAROTIDA: (30.0, 18.0), SUBMAXILAR: (30.0, 18.0), TRAQUEA: (0.0, 0.0),
                 ESOFAGO: (4.0, 2.5), VASOS: (6.0, 3.0)}      # vasos: pool sanguíneo
RETENCION_ADENOMA = 0.8

CASOS = {
    "normal": None,
    "inferior-derecho-12": {"sitio": "inferior-derecho", "diametro_mm": 12, "relacion": 2.5},
    "inferior-izquierdo-8": {"sitio": "inferior-izquierdo", "diametro_mm": 8, "relacion": 2.0},
    "retroesofagico-10": {"sitio": "retroesofagico", "diametro_mm": 10, "relacion": 2.5},
    "mediastinico-15": {"sitio": "mediastinico", "diametro_mm": 15, "relacion": 3.0},
}


def _objetivo(reg: np.ndarray, iso: float, sitio: str, r: float):
    """Punto anatómico al que debe acercarse el centro del adenoma (z, y, x en vóxeles)."""
    tir = reg == TIROIDES
    zs, ys, xs = np.nonzero(tir)
    kz = int(np.median(zs))
    tys, txs = np.nonzero(reg[kz] == TRAQUEA)
    tx = txs.mean()
    if sitio in ("inferior-derecho", "inferior-izquierdo"):
        lado = xs < tx if sitio == "inferior-derecho" else xs > tx
        lz, ly, lx = zs[lado], ys[lado], xs[lado]
        polo = lz <= lz.min() + 10 / iso                     # últimos 10 mm del lóbulo
        return (lz.min() + 1, ly[polo].max(), lx[polo].mean())   # cara posterior del polo inferior
    if sitio == "retroesofagico":
        ez, ey, ex = np.nonzero(reg == ESOFAGO)
        cerca = np.abs(ez - kz) <= 2
        return (kz, ey[cerca].max() + r, ex[cerca].mean())
    if sitio == "mediastinico":
        k = 3 + int(round(8 / iso))
        tys2, txs2 = np.nonzero(reg[k] == TRAQUEA)
        return (k, tys2.mean() - 8 / iso, txs2.mean() - 8 / iso)
    raise ValueError(sitio)


def centro_adenoma(reg: np.ndarray, iso: float, sitio: str, diametro_mm: float, hu: np.ndarray):
    """Centro (z, y, x) del vóxel más cercano al objetivo anatómico donde la esfera entera cae en tejido blando
    libre (−200 a 150 HU, sin tiroides, esófago, tráquea, vasos, hueso ni pulmón). Erosión del tejido libre
    por la esfera = centros admisibles."""
    r = diametro_mm / 2.0 / iso
    libre = (reg == BLANDO) & (hu > -200) & (hu < 150)
    n = int(np.ceil(r))
    zz, yy, xx = np.ogrid[-n:n + 1, -n:n + 1, -n:n + 1]
    bola = (zz ** 2 + yy ** 2 + xx ** 2) <= r * r
    admisible = ndimage.binary_erosion(libre, structure=bola)
    obj = _objetivo(reg, iso, sitio, r)
    cz, cy, cx = np.nonzero(admisible)
    d = np.sqrt((cz - obj[0]) ** 2 + (cy - obj[1]) ** 2 + (cx - obj[2]) ** 2)
    i = int(np.argmin(d))
    return (float(cz[i]), float(cy[i]), float(cx[i])), float(d[i] * iso), [float(v) for v in obj]


def mapas(reg: np.ndarray, iso: float, caso: dict | None, hu: np.ndarray):
    pre = np.zeros(reg.shape, np.float32)
    tar = np.zeros(reg.shape, np.float32)
    for et, (a, b) in CONCENTRACION.items():
        pre[reg == et] = a
        tar[reg == et] = b
    verdad = {"adenoma": None}
    if caso:
        c, desvio_mm, obj = centro_adenoma(reg, iso, caso["sitio"], caso["diametro_mm"], hu)
        r = caso["diametro_mm"] / 2.0 / iso
        z, y, x = np.ogrid[:reg.shape[0], :reg.shape[1], :reg.shape[2]]
        esfera = ((z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2) <= r * r
        pre[esfera] = CONCENTRACION[TIROIDES][0] * caso["relacion"]
        tar[esfera] = CONCENTRACION[TIROIDES][0] * caso["relacion"] * RETENCION_ADENOMA
        verdad["adenoma"] = {**caso, "centro_voxel_zyx": [float(v) for v in c], "centro_mm_zyx": [float(v * iso) for v in c],
                             "volumen_ml": float(esfera.sum() * iso ** 3 / 1000.0), "objetivo_voxel_zyx": obj, "desvio_del_objetivo_mm": desvio_mm,
                             "concentracion_kbq_ml": [CONCENTRACION[TIROIDES][0] * caso["relacion"], CONCENTRACION[TIROIDES][0] * caso["relacion"] * RETENCION_ADENOMA]}
    return pre, tar, verdad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--casos", default=",".join(CASOS))
    a = ap.parse_args()
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    reg = np.load(os.path.join(RAIZ, "salida", "regiones.npz"))["reg"]
    iso = float(f["iso"])
    vox_ml = iso ** 3 / 1000.0
    for nombre in a.casos.split(","):
        pre, tar, verdad = mapas(reg, iso, CASOS[nombre], f["hu"])
        carpeta = os.path.join(RAIZ, "salida", "casos", nombre)
        os.makedirs(carpeta, exist_ok=True)
        np.save(os.path.join(carpeta, "actividad_precoz.npy"), pre)
        np.save(os.path.join(carpeta, "actividad_tardia.npy"), tar)
        verdad.update({"caso": nombre, "actividad_total_MBq": [float(pre.sum() * vox_ml / 1000.0), float(tar.sum() * vox_ml / 1000.0)],
                       "concentraciones_kbq_ml": {str(k): v for k, v in CONCENTRACION.items()}, "retencion_adenoma": RETENCION_ADENOMA})
        json.dump(verdad, open(os.path.join(carpeta, "verdad.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
        ad = verdad["adenoma"]
        print(f"{nombre}: actividad en el campo {verdad['actividad_total_MBq'][0]:.1f} / {verdad['actividad_total_MBq'][1]:.1f} MBq"
              + (f"; adenoma {ad['diametro_mm']} mm ({ad['volumen_ml']:.2f} mL) en z,y,x = {[round(v) for v in ad['centro_mm_zyx']]} mm, a {ad['desvio_del_objetivo_mm']:.1f} mm del objetivo, relación {ad['relacion']}" if ad else ""))


if __name__ == "__main__":
    main()
