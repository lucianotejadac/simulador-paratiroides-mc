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

from segmentar import AIRE, BLANDO, HUESO, PULMON, TIROIDES, PAROTIDA, SUBMAXILAR, TRAQUEA

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# kBq/mL por región: (precoz, tardía)
CONCENTRACION = {AIRE: (0.0, 0.0), BLANDO: (4.0, 2.5), HUESO: (3.0, 2.0), PULMON: (1.0, 0.6), TIROIDES: (40.0, 16.0),
                 PAROTIDA: (30.0, 18.0), SUBMAXILAR: (30.0, 18.0), TRAQUEA: (0.0, 0.0)}
RETENCION_ADENOMA = 0.8

CASOS = {
    "normal": None,
    "inferior-derecho-12": {"sitio": "inferior-derecho", "diametro_mm": 12, "relacion": 2.5},
    "inferior-izquierdo-8": {"sitio": "inferior-izquierdo", "diametro_mm": 8, "relacion": 2.0},
    "retroesofagico-10": {"sitio": "retroesofagico", "diametro_mm": 10, "relacion": 2.5},
    "mediastinico-15": {"sitio": "mediastinico", "diametro_mm": 15, "relacion": 3.0},
}


def centro_adenoma(reg: np.ndarray, iso: float, sitio: str):
    """Centro (z, y, x) en vóxeles. z crece hacia la cabeza; y hacia posterior; x hacia la izquierda del paciente."""
    tir = reg == TIROIDES
    zs, ys, xs = np.nonzero(tir)
    tr = reg == TRAQUEA
    kz = int(np.median(zs))
    tys, txs = np.nonzero(tr[kz])
    ty, tx = tys.mean(), txs.mean()
    if sitio in ("inferior-derecho", "inferior-izquierdo"):
        lado = xs < tx if sitio == "inferior-derecho" else xs > tx
        lz, ly, lx = zs[lado], ys[lado], xs[lado]
        z0 = lz.min()
        polo = lz <= z0 + 4 / iso                      # últimos 8 mm del lóbulo
        return (z0 - 2 / iso, ly[polo].mean() + 7 / iso, lx[polo].mean())
    if sitio == "retroesofagico":
        return (kz - 6 / iso, ty + 26 / iso, tx + 3 / iso)
    if sitio == "mediastinico":
        return (4 + 8 / iso, ty - 8 / iso, tx - 8 / iso)
    raise ValueError(sitio)


def mapas(reg: np.ndarray, iso: float, caso: dict | None):
    pre = np.zeros(reg.shape, np.float32)
    tar = np.zeros(reg.shape, np.float32)
    for et, (a, b) in CONCENTRACION.items():
        pre[reg == et] = a
        tar[reg == et] = b
    verdad = {"adenoma": None}
    if caso:
        c = centro_adenoma(reg, iso, caso["sitio"])
        r = caso["diametro_mm"] / 2.0 / iso
        z, y, x = np.ogrid[:reg.shape[0], :reg.shape[1], :reg.shape[2]]
        esfera = ((z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2) <= r * r
        esfera &= reg != AIRE
        pre[esfera] = CONCENTRACION[TIROIDES][0] * caso["relacion"]
        tar[esfera] = CONCENTRACION[TIROIDES][0] * caso["relacion"] * RETENCION_ADENOMA
        verdad["adenoma"] = {**caso, "centro_voxel_zyx": [float(v) for v in c], "centro_mm_zyx": [float(v * iso) for v in c],
                             "volumen_ml": float(esfera.sum() * iso ** 3 / 1000.0),
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
        pre, tar, verdad = mapas(reg, iso, CASOS[nombre])
        carpeta = os.path.join(RAIZ, "salida", "casos", nombre)
        os.makedirs(carpeta, exist_ok=True)
        np.save(os.path.join(carpeta, "actividad_precoz.npy"), pre)
        np.save(os.path.join(carpeta, "actividad_tardia.npy"), tar)
        verdad.update({"caso": nombre, "actividad_total_MBq": [float(pre.sum() * vox_ml / 1000.0), float(tar.sum() * vox_ml / 1000.0)],
                       "concentraciones_kbq_ml": {str(k): v for k, v in CONCENTRACION.items()}, "retencion_adenoma": RETENCION_ADENOMA})
        json.dump(verdad, open(os.path.join(carpeta, "verdad.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
        ad = verdad["adenoma"]
        print(f"{nombre}: actividad en el campo {verdad['actividad_total_MBq'][0]:.1f} / {verdad['actividad_total_MBq'][1]:.1f} MBq"
              + (f"; adenoma {ad['diametro_mm']} mm ({ad['volumen_ml']:.2f} mL) en z,y,x = {[round(v) for v in ad['centro_mm_zyx']]} mm, relación {ad['relacion']}" if ad else ""))


if __name__ == "__main__":
    main()
