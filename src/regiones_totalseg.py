"""Paso 2b: reemplaza las regiones por umbral con las estructuras de TotalSegmentator (Wasserthal et al., 2023).

De TotalSegmentator se toman tiroides, tráquea, esófago, vértebras, clavículas, pulmones y vasos
(carótidas comunes, subclavias, tronco braquiocefálico, venas braquiocefálicas). Lo que TotalSegmentator
no cubre (vía aérea alta y laringe, glándulas salivales, el resto del hueso) se conserva de la
segmentación por umbral. La tiroides por umbral queda descartada: incluía cartílago laríngeo y estaba
corrida 16 mm hacia craneal (BITACORA 0002).

Entrada: salida/regiones.npz (umbral) y salida/totalseg/etiquetas.npz. Salida: salida/regiones.npz
(se guarda la anterior como regiones_umbral.npz) y montaje salida/regiones_totalseg.png.
"""
from __future__ import annotations

import json
import os
import shutil

import numpy as np

from segmentar import AIRE, BLANDO, HUESO, PULMON, TIROIDES, TRAQUEA, ESOFAGO, VASOS, NOMBRES, montaje

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    from totalsegmentator.map_to_binary import class_map
    nombres = {v: int(k) for k, v in class_map["total"].items()}
    s = os.path.join(RAIZ, "salida")
    if not os.path.exists(os.path.join(s, "regiones_umbral.npz")):
        shutil.copy(os.path.join(s, "regiones.npz"), os.path.join(s, "regiones_umbral.npz"))
    reg = np.load(os.path.join(s, "regiones_umbral.npz"))["reg"].copy()
    et = np.load(os.path.join(s, "totalseg", "etiquetas.npz"))["et"]
    f = np.load(os.path.join(s, "fantoma.npz"))
    hu, iso = f["hu"], float(f["iso"])

    def m(*ns):
        return np.isin(et, [nombres[n] for n in ns if n in nombres])

    reg[reg == TIROIDES] = BLANDO                                  # descartar la tiroides por umbral
    reg[m("thyroid_gland")] = TIROIDES
    reg[m("esophagus")] = ESOFAGO
    reg[m("common_carotid_artery_left", "common_carotid_artery_right", "subclavian_artery_left", "subclavian_artery_right",
          "brachiocephalic_trunk", "brachiocephalic_vein_left", "brachiocephalic_vein_right")] = VASOS
    reg[m("trachea") & (hu < -300)] = TRAQUEA
    reg[m("lung_upper_lobe_left", "lung_upper_lobe_right") & (hu < -300)] = PULMON
    hueso_ts = np.isin(et, [v for k, v in nombres.items() if k.startswith("vertebrae_") or k.startswith("clavicula")])
    reg[hueso_ts & (hu > 100)] = HUESO
    np.savez_compressed(os.path.join(s, "regiones.npz"), reg=reg)
    vol = {NOMBRES[i]: round(float((reg == i).sum()) * iso ** 3 / 1000.0, 1) for i in range(len(NOMBRES))}
    zt = np.nonzero((reg == TIROIDES).any(axis=(1, 2)))[0]
    montaje(hu.astype(np.float32), reg, os.path.join(s, "regiones_totalseg.png"), int(zt.mean()))
    json.dump({"volumen_ml": vol, "nombres": NOMBRES, "fuente": "TotalSegmentator 2 (total, roi_subset) + umbral para lo no cubierto",
               "tiroides_cortes": [int(zt.min()), int(zt.max())]}, open(os.path.join(s, "regiones.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print("volúmenes (mL):", vol, "| tiroides en cortes", zt.min(), "-", zt.max())


if __name__ == "__main__":
    main()
