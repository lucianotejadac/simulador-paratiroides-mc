"""Compara la segmentación de TotalSegmentator con la propia y revisa dónde cae cada adenoma.

Lee salida/totalseg/total.nii.gz (multietiqueta, misma grilla que el fantoma), lo lleva a (z, y, x),
y reporta: volumen y Dice de la tiroides, posición del esófago, estructuras que ocupa cada adenoma y
su distancia a tiroides, esófago, tráquea y vértebras. Escribe salida/totalseg/etiquetas.npz y un
montaje de control salida/totalseg/comparacion.png.
"""
from __future__ import annotations

import json
import os

import nibabel as nib
import numpy as np
from scipy import ndimage

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cargar():
    from totalsegmentator.map_to_binary import class_map
    img = nib.load(os.path.join(RAIZ, "salida", "totalseg", "total.nii.gz"))
    et = np.transpose(np.asarray(img.dataobj).astype(np.uint8), (2, 1, 0))      # (x,y,z) -> (z,y,x)
    nombres = class_map["total"]
    return et, {int(k): v for k, v in nombres.items()}


def dice(a, b):
    return 2.0 * (a & b).sum() / max(1, a.sum() + b.sum())


def main():
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    hu, iso = f["hu"], float(f["iso"])
    reg = np.load(os.path.join(RAIZ, "salida", "regiones.npz"))["reg"]
    et, nombres = cargar()
    assert et.shape == hu.shape, (et.shape, hu.shape)
    ml = iso ** 3 / 1000.0
    presentes = {nombres[i]: round(float((et == i).sum()) * ml, 1) for i in np.unique(et) if i}
    print("TotalSegmentator, volúmenes (mL):", presentes)
    inv = {v: k for k, v in nombres.items()}
    tir_ts = et == inv["thyroid_gland"]
    tir_mio = reg == 4
    print(f"tiroides: TotalSegmentator {tir_ts.sum() * ml:.1f} mL, propia {tir_mio.sum() * ml:.1f} mL, Dice {dice(tir_ts, tir_mio):.2f}; "
          f"cortes TS {np.nonzero(tir_ts.any(axis=(1, 2)))[0][[0, -1]]}, propios {np.nonzero(tir_mio.any(axis=(1, 2)))[0][[0, -1]]}")
    tr_ts, tr_mio = et == inv["trachea"], reg == 7
    print(f"tráquea: Dice {dice(tr_ts, tr_mio):.2f}")
    eso = et == inv["esophagus"]
    vert = np.isin(et, [inv[n] for n in nombres.values() if n.startswith("vertebrae_")])
    dist = {n: ndimage.distance_transform_edt(~m) * iso for n, m in (("tiroides", tir_ts), ("esófago", eso), ("tráquea", tr_ts), ("vértebra", vert))}
    resumen = {}
    for caso in ["inferior-derecho-12", "inferior-izquierdo-8", "retroesofagico-10", "mediastinico-15"]:
        v = json.load(open(os.path.join(RAIZ, "salida", "casos", caso, "verdad.json"), encoding="utf-8"))["adenoma"]
        c = v["centro_voxel_zyx"]
        r = v["diametro_mm"] / 2 / iso
        z, y, x = np.ogrid[:hu.shape[0], :hu.shape[1], :hu.shape[2]]
        esf = ((z - c[0]) ** 2 + (y - c[1]) ** 2 + (x - c[2]) ** 2) <= r * r
        ocupa = {nombres[i]: int((et[esf] == i).sum()) for i in np.unique(et[esf]) if i}
        ci = tuple(int(round(t)) for t in c)
        d = {n: round(float(m[ci]), 1) for n, m in dist.items()}
        resumen[caso] = {"vóxeles": int(esf.sum()), "ocupa": ocupa, "distancia_centro_mm": d}
        print(f"{caso}: {int(esf.sum())} vóxeles; ocupa {ocupa or 'solo tejido no etiquetado'}; centro a {d} mm")
    np.savez_compressed(os.path.join(RAIZ, "salida", "totalseg", "etiquetas.npz"), et=et)
    json.dump({"volumenes_ml": presentes, "adenomas": resumen}, open(os.path.join(RAIZ, "salida", "totalseg", "comparacion.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    montaje(hu, et, inv, tir_mio, os.path.join(RAIZ, "salida", "totalseg", "comparacion.png"))


def montaje(hu, et, inv, tir_mio, ruta):
    from PIL import Image
    colores = {"thyroid_gland": (255, 70, 70), "esophagus": (80, 220, 120), "trachea": (255, 230, 0), "common_carotid_artery_left": (230, 60, 230),
               "common_carotid_artery_right": (230, 60, 230)}
    zs = np.nonzero((et == inv["thyroid_gland"]).any(axis=(1, 2)))[0]
    ks = np.linspace(max(0, zs.min() - 6), zs.max() + 6, 8).astype(int)
    tr = np.nonzero(et[int(zs.mean())] == inv["trachea"])
    cy, cx = int(tr[0].mean()), int(tr[1].mean())
    ims = []
    for k in ks:
        sl = (slice(cy - 30, cy + 30), slice(cx - 40, cx + 40))
        g = np.clip((hu[k][sl] + 160) / 400 * 255, 0, 255).astype(np.float32)
        rgb = np.stack([g, g, g], -1)
        for n, col in colores.items():
            if n in inv:
                m = et[k][sl] == inv[n]
                rgb[m] = 0.5 * rgb[m] + 0.5 * np.array(col)
        borde = tir_mio[k][sl] & ~ndimage.binary_erosion(tir_mio[k][sl])
        rgb[borde] = (60, 160, 255)
        ims.append(Image.fromarray(rgb.astype(np.uint8)).resize((240, 180), Image.NEAREST))
    L = Image.new("RGB", (4 * 244, 2 * 184), (0, 0, 0))
    for i, im in enumerate(ims):
        L.paste(im, ((i % 4) * 244, (i // 4) * 184))
    L.save(ruta)
    print("montaje", ruta, "cortes", list(ks))


if __name__ == "__main__":
    main()
