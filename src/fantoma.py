"""Paso 1: del CT de TCIA al fantoma de atenuación.

Lee la serie CT, recorta el cuello (desde la mandíbula hasta ~12 cm bajo el mínimo del
cuello, que cubre clavículas y mediastino alto), remuestrea a vóxeles isotrópicos y
convierte unidades Hounsfield a coeficiente de atenuación lineal a 140 keV con la
conversión bilineal (agua 0.1537 /cm, hueso cortical ~0.28 /cm a 140 keV).

Salida: salida/fantoma.npz con hu (int16), mu (float32, 1/cm), espaciado, origen y z0 del
CT original, para poder exportar después con la misma identidad geométrica.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pydicom
from scipy import ndimage

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CT_POR_DEFECTO = r"C:\Users\lucia\Downloads\PET CT\ENTREGA PET\PET Magdalena\Caso 2\CT"
MU_AGUA_140 = 0.1537     # 1/cm, NIST, 140 keV
MU_HUESO_140 = 0.2860    # 1/cm, hueso cortical ICRU a 140 keV (densidad 1.92)
HU_HUESO_REF = 1000.0    # HU que se asigna a MU_HUESO_140 en la rama de hueso


def hu_a_mu(hu: np.ndarray) -> np.ndarray:
    """Conversión bilineal HU -> mu(140 keV) en 1/cm. Aire -1000 -> 0, agua 0 -> mu_agua, hueso lineal aparte."""
    hu = hu.astype(np.float32)
    mu = np.where(hu <= 0.0, MU_AGUA_140 * (1.0 + hu / 1000.0),
                  MU_AGUA_140 + (MU_HUESO_140 - MU_AGUA_140) * hu / HU_HUESO_REF)
    return np.clip(mu, 0.0, 0.6).astype(np.float32)


def leer_ct(carpeta: str):
    cortes = []
    for a in os.listdir(carpeta):
        d = pydicom.dcmread(os.path.join(carpeta, a), force=True)
        if getattr(d, "Modality", "") == "CT" and hasattr(d, "PixelData"):
            cortes.append(d)
    cortes.sort(key=lambda d: float(d.ImagePositionPatient[2]))
    hu = np.stack([d.pixel_array.astype(np.float32) * float(d.RescaleSlope) + float(d.RescaleIntercept) for d in cortes])
    z = np.array([float(d.ImagePositionPatient[2]) for d in cortes])
    d0 = cortes[0]
    ps = [float(x) for x in d0.PixelSpacing]          # [fila(y), columna(x)]
    origen = [float(v) for v in d0.ImagePositionPatient]
    return hu, z, ps, origen, cortes


def recortar(hu, z, ps, abajo_mm=120.0, arriba_mm=70.0):
    """Localiza el cuello (mínimo de área corporal en el 40 % superior) y recorta en z y en el plano."""
    area = (hu > -300).sum(axis=(1, 2))
    n = len(area)
    i0 = int(n * 0.6)
    i_cuello = i0 + int(np.argmin(area[i0:]))
    dz = float(np.median(np.diff(z)))
    a = max(0, i_cuello - int(round(abajo_mm / dz)))
    b = min(n, i_cuello + int(round(arriba_mm / dz)) + 1)
    sub = hu[a:b]
    cuerpo = ndimage.binary_opening(sub > -300, iterations=2)
    # quedarse con el componente conexo mayor (saca la camilla y el soporte de cabeza)
    et, k = ndimage.label(cuerpo)
    if k > 1:
        tam = ndimage.sum(cuerpo, et, range(1, k + 1))
        cuerpo = et == (1 + int(np.argmax(tam)))
    ys, xs = np.nonzero(cuerpo.any(axis=0))
    margen = int(round(15.0 / ps[0]))
    y0, y1 = max(0, ys.min() - margen), min(sub.shape[1], ys.max() + margen + 1)
    x0, x1 = max(0, xs.min() - margen), min(sub.shape[2], xs.max() + margen + 1)
    sub = np.where(cuerpo, sub, -1000.0)[:, y0:y1, x0:x1]
    return sub, (a, b, y0, y1, x0, x1), i_cuello


def remuestrear(vol, dz, ps, iso=2.0):
    factores = (dz / iso, ps[0] / iso, ps[1] / iso)
    return ndimage.zoom(vol, factores, order=1).astype(np.float32)


def montaje(hu, iso, ruta):
    from PIL import Image

    def ventana(a, c=40, w=400):
        return np.clip((a - (c - w / 2)) / w * 255, 0, 255).astype(np.uint8)

    nz, ny, nx = hu.shape
    cor = ventana(hu[:, ny // 2, :])[::-1]
    sag = ventana(hu[:, :, nx // 2])[::-1]
    axiales = [ventana(hu[k]) for k in np.linspace(5, nz - 6, 6).astype(int)]
    W = cor.shape[1] + sag.shape[1] + 10
    H = max(cor.shape[0], sag.shape[0]) + ny + 10
    lienzo = Image.new("L", (max(W, 6 * (nx + 4)), H), 0)
    lienzo.paste(Image.fromarray(cor), (0, 0))
    lienzo.paste(Image.fromarray(sag), (cor.shape[1] + 10, 0))
    for k, a in enumerate(axiales):
        lienzo.paste(Image.fromarray(a), (k * (nx + 4), max(cor.shape[0], sag.shape[0]) + 10))
    lienzo.save(ruta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ct", default=CT_POR_DEFECTO)
    ap.add_argument("--iso", type=float, default=2.0, help="mm por vóxel, isotrópico")
    ap.add_argument("--abajo", type=float, default=120.0, help="mm bajo el mínimo del cuello")
    ap.add_argument("--arriba", type=float, default=70.0, help="mm sobre el mínimo del cuello")
    ap.add_argument("--salida", default=os.path.join(RAIZ, "salida"))
    a = ap.parse_args()
    os.makedirs(a.salida, exist_ok=True)
    hu, z, ps, origen, cortes = leer_ct(a.ct)
    dz = float(np.median(np.diff(z)))
    sub, (i_a, i_b, y0, y1, x0, x1), i_cuello = recortar(hu, z, ps, a.abajo, a.arriba)
    iso = remuestrear(sub, dz, ps, a.iso)
    mu = hu_a_mu(iso)
    origen_rec = [origen[0] + x0 * ps[1], origen[1] + y0 * ps[0], z[i_a]]
    meta = {"ct": a.ct, "paciente": str(getattr(cortes[0], "PatientID", "")), "iso_mm": a.iso, "forma_zyx": list(iso.shape),
            "origen_mm": origen_rec, "recorte_indices": [int(v) for v in (i_a, i_b, y0, y1, x0, x1)], "z_cuello_mm": float(z[i_cuello]),
            "ct_espaciado_mm": [dz, ps[0], ps[1]], "ct_cortes": len(z), "mu_agua_140": MU_AGUA_140,
            "study_uid": str(getattr(cortes[0], "StudyInstanceUID", "")), "frame_uid": str(getattr(cortes[0], "FrameOfReferenceUID", "")),
            "sexo": str(getattr(cortes[0], "PatientSex", "")), "edad": str(getattr(cortes[0], "PatientAge", ""))}
    np.savez_compressed(os.path.join(a.salida, "fantoma.npz"), hu=np.round(iso).astype(np.int16), mu=mu,
                        iso=a.iso, origen=np.array(origen_rec))
    json.dump(meta, open(os.path.join(a.salida, "fantoma.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    montaje(iso, a.iso, os.path.join(a.salida, "fantoma.png"))
    print(f"fantoma {iso.shape} a {a.iso} mm ({iso.shape[0] * a.iso / 10:.0f} x {iso.shape[1] * a.iso / 10:.0f} x {iso.shape[2] * a.iso / 10:.0f} cm); "
          f"mu: aire {mu.min():.3f}, mediana cuerpo {np.median(mu[iso > -300]):.3f}, máx {mu.max():.3f} /cm; cuello en z={z[i_cuello]:.0f}")


if __name__ == "__main__":
    main()
