"""Paso 6: reconstrucción OSEM de las proyecciones simuladas (para revisar y para exportar un volumen NM).

Proyector: rotación del volumen con interpolación bilineal (scipy.ndimage.rotate) y suma a lo largo del
eje perpendicular al detector; opcionalmente con corrección de atenuación (mapa mu del fantoma,
remuestreado a la grilla de reconstrucción) y con una PSF gaussiana fija (resolución al radio medio).
No es el proyector del Monte Carlo: es el que un equipo clínico usaría, con sus simplificaciones.
Grilla de reconstrucción: matriz × matriz × matriz con el píxel de adquisición.
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
from scipy import ndimage

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _rotar(vol, grados):
    return ndimage.rotate(vol, grados, axes=(2, 1), reshape=False, order=1, mode="constant", cval=0.0, prefilter=False)


def proyectar(vol, angulos_rad, mu=None, pix_cm=0.33, fwhm_pix=0.0):
    """vol (z, y, x) -> proyecciones (ang, z, u). El detector en el ángulo phi mira desde la dirección (cos phi, sin phi):
    se rota el volumen -phi para que esa dirección quede en +x y se integra en x."""
    proy = np.zeros((len(angulos_rad), vol.shape[0], vol.shape[2]), np.float32)
    for ia, phi in enumerate(angulos_rad):
        g = np.degrees(phi)           # el sentido de giro de ndimage.rotate con axes=(2,1) es el contrario al de la cámara
        v = _rotar(vol, g)
        if mu is not None:
            m = _rotar(mu, g)
            # atenuación acumulada desde cada vóxel hacia +x (hacia el detector)
            tau = np.cumsum(m[:, :, ::-1], axis=2)[:, :, ::-1] * pix_cm
            v = v * np.exp(-(tau - 0.5 * m * pix_cm))
        p = v.sum(axis=2)          # (z, y) ; el eje y del volumen rotado es el eje u del detector
        if fwhm_pix > 0:
            p = ndimage.gaussian_filter(p, fwhm_pix / 2.3548)
        proy[ia] = p
    return proy


def retroproyectar(proy, angulos_rad, forma, mu=None, pix_cm=0.33, fwhm_pix=0.0):
    vol = np.zeros(forma, np.float32)
    for ia, phi in enumerate(angulos_rad):
        p = proy[ia]
        if fwhm_pix > 0:
            p = ndimage.gaussian_filter(p, fwhm_pix / 2.3548)
        v = np.repeat(p[:, :, None], forma[2], axis=2).astype(np.float32)
        if mu is not None:
            m = _rotar(mu, np.degrees(phi))
            tau = np.cumsum(m[:, :, ::-1], axis=2)[:, :, ::-1] * pix_cm
            v = v * np.exp(-(tau - 0.5 * m * pix_cm))
        vol += _rotar(v, -np.degrees(phi))
    return vol


def osem(proy, angulos_rad, iteraciones=4, subconjuntos=6, mu=None, pix_cm=0.33, fwhm_pix=0.0):
    n_ang, nz, nu = proy.shape
    forma = (nz, nu, nu)
    x = np.ones(forma, np.float32)
    orden = np.arange(n_ang)
    for it in range(iteraciones):
        for s in range(subconjuntos):
            idx = orden[s::subconjuntos]
            ang = angulos_rad[idx]
            est = proyectar(x, ang, mu, pix_cm, fwhm_pix)
            razon = proy[idx] / np.maximum(est, 1e-3)
            razon[est < 1e-3] = 0.0
            num = retroproyectar(razon, ang, forma, mu, pix_cm, fwhm_pix)
            den = retroproyectar(np.ones_like(razon), ang, forma, mu, pix_cm, fwhm_pix)
            x = x * num / np.maximum(den, 1e-3)
            x[den < 1e-3] = 0.0
    return x


def mu_en_grilla(mu_fantoma, iso_mm, matriz, pix_mm):
    """Lleva el mapa mu (z,y,x a iso_mm) a la grilla de reconstrucción (matriz³ a pix_mm), centrado."""
    nz, ny, nx = mu_fantoma.shape
    esc = iso_mm / pix_mm
    z = ndimage.zoom(mu_fantoma, esc, order=1)
    out = np.zeros((matriz, matriz, matriz), np.float32)
    cz, cy, cx = [(m - s) // 2 for m, s in zip((matriz, matriz, matriz), z.shape)]
    sz = [slice(max(0, c), max(0, c) + min(s, m)) for c, s, m in zip((cz, cy, cx), z.shape, (matriz,) * 3)]
    so = [slice(max(0, -c), max(0, -c) + (sl.stop - sl.start)) for c, sl in zip((cz, cy, cx), sz)]
    out[sz[0], sz[1], sz[2]] = z[so[0], so[1], so[2]]
    return out


def montaje(vol, ruta, filas=3, cols=8):
    from PIL import Image
    nz = vol.shape[0]
    idx = np.linspace(nz * 0.2, nz * 0.8, filas * cols).astype(int)
    tope = np.percentile(vol, 99.9) or 1.0
    ims = [Image.fromarray(np.clip(vol[i] / tope * 255, 0, 255).astype(np.uint8)) for i in idx]
    w, h = ims[0].size
    L = Image.new("L", (cols * (w + 2), filas * (h + 2)), 0)
    for k, im in enumerate(ims):
        L.paste(im, ((k % cols) * (w + 2), (k // cols) * (h + 2)))
    L.save(ruta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="inferior-derecho-12")
    ap.add_argument("--fase", default="precoz")
    ap.add_argument("--iteraciones", type=int, default=4)
    ap.add_argument("--subconjuntos", type=int, default=6)
    ap.add_argument("--sin-ac", action="store_true", help="sin corrección de atenuación")
    ap.add_argument("--psf", type=float, default=0.0, help="FWHM de la PSF del proyector en mm (0 = sin PSF)")
    ap.add_argument("--filtro", type=float, default=6.0, help="FWHM en mm del filtro gaussiano posterior (0 = ninguno)")
    a = ap.parse_args()
    carpeta = os.path.join(RAIZ, "salida", "casos", a.caso)
    d = np.load(os.path.join(carpeta, f"proyecciones_{a.fase}.npz"))
    proy = d["ruido"].astype(np.float32)
    ang = d["angulos_rad"]
    pix_mm = float(d["pixel_mm"])
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    mu = None if a.sin_ac else mu_en_grilla(f["mu"], float(f["iso"]), int(d["matriz"]), pix_mm)
    vol = osem(proy, ang, a.iteraciones, a.subconjuntos, mu, pix_mm / 10.0, a.psf / pix_mm)
    if a.filtro > 0:
        vol = ndimage.gaussian_filter(vol, a.filtro / 2.3548 / pix_mm)
    sufijo = f"{a.fase}{'_noac' if a.sin_ac else '_ac'}"
    np.save(os.path.join(carpeta, f"recon_{sufijo}.npy"), vol)
    montaje(vol, os.path.join(carpeta, f"recon_{sufijo}.png"))
    json.dump({"caso": a.caso, "fase": a.fase, "iteraciones": a.iteraciones, "subconjuntos": a.subconjuntos, "ac": not a.sin_ac, "psf_mm": a.psf, "filtro_mm": a.filtro,
               "pixel_mm": pix_mm, "max": float(vol.max())}, open(os.path.join(carpeta, f"recon_{sufijo}.json"), "w", encoding="utf-8"), indent=2)
    print(f"recon {a.caso} {sufijo}: {vol.shape}, máx {vol.max():.1f}")


if __name__ == "__main__":
    main()
