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
    ap.add_argument("--radio", type=float, default=0.0, help="cm, órbita circular (0 = contorno del cuerpo + 2 cm)")
    ap.add_argument("--matriz", type=int, default=128)
    ap.add_argument("--pixel", type=float, default=4.8, help="mm (128 sin zoom: 61 cm de campo)")
    ap.add_argument("--tiempo", type=float, default=25.0, help="s por proyección")
    ap.add_argument("--paso", type=float, default=0.1, help="cm, paso de integración de la transmisión")
    ap.add_argument("--planar", type=int, default=4_000_000, help="historias de la planar anterior por fase (0 = sin planar)")
    ap.add_argument("--tiempo-planar", type=float, default=300.0, help="s")
    a = ap.parse_args()
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    mu, hu, iso = f["mu"], f["hu"], float(f["iso"])
    carpeta = os.path.join(RAIZ, "salida", "casos", a.caso)
    vox_cm = iso / 10.0
    cuerpo = mu > 0.02
    yy, xx = np.nonzero(cuerpo.any(axis=0))
    cy, cx = mu.shape[1] * vox_cm / 2, mu.shape[2] * vox_cm / 2
    r_cuerpo = float(np.sqrt(((yy + 0.5) * vox_cm - cy) ** 2 + ((xx + 0.5) * vox_cm - cx) ** 2).max())
    if a.radio <= 0:
        a.radio = round(r_cuerpo + 2.0, 1)
    anterior_cm = float(cy - (yy.min() * vox_cm))                  # del eje al punto más anterior del cuerpo
    cam = mc.Camara(a.angulos, a.arco, 0.0, a.radio, a.matriz, a.pixel, a.tiempo)
    cam_planar = mc.Camara(1, 360.0, 270.0, round(anterior_cm + 1.0, 1), 256, 2.4, a.tiempo_planar)   # 270°: detector anterior
    registro = {"caso": a.caso, "radio_auto_cm": a.radio, "camara": {"angulos": a.angulos, "arco_grados": a.arco, "radio_cm": a.radio, "matriz": a.matriz, "pixel_mm": a.pixel,
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
        if a.planar > 0:
            t1 = time.time()
            esp_p, ruido_p, est_p = mc.simular(act, mu, hu, iso / 10.0, cam_planar, n_hist=a.planar, semilla=a.semilla + (500 if fase == "precoz" else 1500), paso_cm=a.paso)
            est_p["segundos"] = round(time.time() - t1, 1)
            registro["fases"][fase]["planar"] = {**est_p, "radio_cm": cam_planar.radio_cm, "matriz": 256, "pixel_mm": 2.4, "tiempo_s": a.tiempo_planar}
            np.savez_compressed(os.path.join(carpeta, f"planar_{fase}.npz"), esperado=esp_p[0].astype(np.float32), ruido=ruido_p[0], pixel_mm=2.4,
                                distancia_cm=cam_planar.radio_cm, tiempo_s=a.tiempo_planar)
            print(f"   planar anterior {fase}: {est_p['segundos']} s, {est_p['cuentas_medias_por_proyeccion']:.0f} cuentas", flush=True)
        print(f"{a.caso} {fase}: {est['segundos']} s, {est['historias']} historias, {est['compton']} Compton, {est['fotoelectricos']} fotoeléctricos; "
              f"{est['actividad_MBq']:.1f} MBq en el campo, {est['cuentas_medias_por_proyeccion']:.0f} cuentas por proyección (media), "
              f"máximo por píxel {ruido.max()}", flush=True)
    json.dump(registro, open(os.path.join(carpeta, "adquisicion.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
