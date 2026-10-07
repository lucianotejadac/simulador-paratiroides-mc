"""Paso 5: adquisición SPECT de un caso (fase precoz y tardía) con el motor Monte Carlo.

Salida en salida/casos/<caso>/: proyecciones_<fase>.npz (esperado, ruido, ángulos, geometría),
adquisicion.json (parámetros y estadísticas) y proyecciones_<fase>.png (montaje de ángulos).
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

import montecarlo as mc

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def montaje(ruido: np.ndarray, ruta: str, cuantos: int = 8):
    from PIL import Image
    n = ruido.shape[0]
    idx = np.linspace(0, n - 1, cuantos).astype(int)
    tope = np.percentile(ruido, 99.8)
    ims = [Image.fromarray(np.clip(ruido[i][::-1] / tope * 255, 0, 255).astype(np.uint8)).resize((ruido.shape[2] * 2, ruido.shape[1] * 2)) for i in idx]
    W = sum(i.width + 4 for i in ims)
    L = Image.new("L", (W, ims[0].height), 0)
    x = 0
    for im in ims:
        L.paste(im, (x, 0))
        x += im.width + 4
    L.save(ruta)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="inferior-derecho-12")
    ap.add_argument("--fases", default="precoz,tardia")
    ap.add_argument("--historias", type=int, default=2_000_000)
    ap.add_argument("--semilla", type=int, default=1)
    ap.add_argument("--angulos", type=int, default=60)
    ap.add_argument("--arco", type=float, default=360.0)
    ap.add_argument("--radio", type=float, default=20.0, help="cm, órbita circular")
    ap.add_argument("--matriz", type=int, default=128)
    ap.add_argument("--pixel", type=float, default=3.3, help="mm")
    ap.add_argument("--tiempo", type=float, default=25.0, help="s por proyección")
    ap.add_argument("--paso", type=float, default=0.1, help="cm, paso de integración de la transmisión")
    a = ap.parse_args()
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    mu, hu, iso = f["mu"], f["hu"], float(f["iso"])
    carpeta = os.path.join(RAIZ, "salida", "casos", a.caso)
    cam = mc.Camara(a.angulos, a.arco, 0.0, a.radio, a.matriz, a.pixel, a.tiempo)
    registro = {"caso": a.caso, "camara": {"angulos": a.angulos, "arco_grados": a.arco, "radio_cm": a.radio, "matriz": a.matriz, "pixel_mm": a.pixel,
                                           "tiempo_s": a.tiempo, "colimador": "LEHR", "ventana_keV": list(mc.VENTANA), "resolucion_energia": mc.RES_ENERGIA_140,
                                           "fwhm_intrinseca_mm": mc.FWHM_INTRINSECA, "g_geom": mc.G_GEOM},
                "historias": a.historias, "semilla": a.semilla, "fases": {}}
    for fase in a.fases.split(","):
        act = np.load(os.path.join(carpeta, f"actividad_{fase}.npy"))
        t0 = time.time()
        esperado, ruido, est = mc.simular(act, mu, hu, iso / 10.0, cam, n_hist=a.historias, semilla=a.semilla + (0 if fase == "precoz" else 1000), paso_cm=a.paso)
        est["segundos"] = round(time.time() - t0, 1)
        registro["fases"][fase] = est
        np.savez_compressed(os.path.join(carpeta, f"proyecciones_{fase}.npz"), esperado=esperado.astype(np.float32), ruido=ruido,
                            angulos_rad=cam.angulos, radio_cm=a.radio, pixel_mm=a.pixel, matriz=a.matriz)
        montaje(ruido, os.path.join(carpeta, f"proyecciones_{fase}.png"))
        print(f"{a.caso} {fase}: {est['segundos']} s, {est['historias']} historias, {est['compton']} Compton, {est['fotoelectricos']} fotoeléctricos; "
              f"{est['actividad_MBq']:.1f} MBq en el campo, {est['cuentas_medias_por_proyeccion']:.0f} cuentas por proyección (media), "
              f"máximo por píxel {ruido.max()}", flush=True)
    json.dump(registro, open(os.path.join(carpeta, "adquisicion.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
