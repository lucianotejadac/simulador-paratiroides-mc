"""Aplica un archivo de correcciones de contornos descargado de la página a las regiones del fantoma.

En la página los contornos viven en los cortes de la grilla SPECT (píxel de 4.8 mm, cortes de 4.8 mm) a 2
sub-píxeles por píxel (2.4 mm), recortados al cuerpo (ver contornos_web.py). El archivo trae, por corte, solo
los píxeles cambiados como corridas [inicio, largo, etiqueta]. Para cada vóxel del fantoma (2 mm) se busca el
sub-píxel de la página que lo contiene, con la misma correspondencia que usó contornos_web.py, y si ese
sub-píxel fue corregido el vóxel toma la etiqueta nueva. Lo no corregido no se toca.

Uso:
    python src/importar_correcciones.py contornos-paratiroides-v2.json            # aplica y resume
    python src/importar_correcciones.py contornos-paratiroides-v2.json --rehacer  # y rehace actividad, simulación y página
La versión anterior queda en salida/regiones_antes_de_<fecha>.npz. Después de importar no se debe volver a correr
regiones_totalseg.py, que regeneraría las regiones desde TotalSegmentator y borraría las correcciones.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import subprocess
import sys

import numpy as np

from contornos_web import grilla
from segmentar import NOMBRES

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archivo")
    ap.add_argument("--rehacer", action="store_true")
    a = ap.parse_args()
    d = json.load(open(a.archivo, encoding="utf-8"))
    if d.get("formato") != "contornos-editados-v1" or d.get("meta", {}).get("pagina") != "simulador-paratiroides-mc":
        sys.exit("el archivo no es de correcciones de paratiroides")
    m = d["meta"]
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    iso = float(f["iso"])
    ruta = os.path.join(RAIZ, "salida", "regiones.npz")
    reg = np.load(ruta)["reg"]
    idx = json.load(open(os.path.join(RAIZ, "docs", "datos", "indice.json"), encoding="utf-8"))
    matriz, pix = grilla()
    esc = iso / pix
    desplaz = [(matriz - int(round(s * esc))) // 2 for s in reg.shape]
    sub = int(m["sub"])
    y0, x0, h, w = m["recorte_yx"]
    z0 = idx["z0"]
    # corrección por corte de la página: arreglo (h, w) con -1 donde no hubo cambio
    corr = {}
    for k, corridas in d["cortes"].items():
        plano = np.full(h * w, -1, np.int16)
        for ini, largo, et in corridas:
            plano[ini:ini + largo] = et
        corr[int(k)] = plano.reshape(h, w)
    # correspondencia vóxel del fantoma -> (corte, fila, columna) de la página
    gz = np.arange(reg.shape[0]) * esc + desplaz[0]
    gy = np.arange(reg.shape[1]) * esc + desplaz[1]
    gx = np.arange(reg.shape[2]) * esc + desplaz[2]
    kz = np.rint(gz).astype(int) - z0
    fy = np.rint((gy + 0.5) * sub - 0.5).astype(int) - y0
    fx = np.rint((gx + 0.5) * sub - 0.5).astype(int) - x0
    vy = (fy >= 0) & (fy < h)
    vx = (fx >= 0) & (fx < w)
    nuevo = reg.copy()
    tocados = 0
    for pz in range(reg.shape[0]):
        c = corr.get(int(kz[pz]))
        if c is None:
            continue
        sub_c = np.full(reg.shape[1:], -1, np.int16)
        sub_c[np.ix_(vy, vx)] = c[np.ix_(fy[vy], fx[vx])]
        cambia = sub_c >= 0
        nuevo[pz][cambia] = sub_c[cambia]
        tocados += int(cambia.sum())
    copia = os.path.join(RAIZ, "salida", f"regiones_antes_de_{dt.datetime.now():%Y%m%d_%H%M%S}.npz")
    shutil.copy(ruta, copia)
    np.savez_compressed(ruta, reg=nuevo)
    ml = iso ** 3 / 1000.0
    print(f"{len(corr)} cortes de la página corregidos -> {tocados} vóxeles del fantoma; copia anterior en {copia}")
    for i, n in enumerate(NOMBRES):
        dv = (int((nuevo == i).sum()) - int((reg == i).sum())) * ml
        if abs(dv) > 0.001:
            print(f"   {n:24s} {int((reg == i).sum()) * ml:8.1f} mL -> {int((nuevo == i).sum()) * ml:8.1f} mL ({dv:+.1f})")
    if a.rehacer:
        subprocess.run([sys.executable, os.path.join(RAIZ, "src", "actividad.py")], check=True)
        subprocess.run([sys.executable, os.path.join(RAIZ, "scratch", "rehacer.py")], check=True, cwd=RAIZ)


if __name__ == "__main__":
    main()
