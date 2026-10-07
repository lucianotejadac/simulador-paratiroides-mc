"""Paso 9: datos compactos para la página de GitHub Pages (docs/).

Por caso escribe en docs/datos/<caso>/:
  nm_precoz.bin     uint8 (nz, 128, 128): reconstrucción OSEM filtrada, escalada a su máximo
  nm_tardia.bin     idem
  proy_precoz.bin   uint8 (60, 128, 128): proyecciones con ruido, escaladas al percentil 99.8
  proy_tardia.bin   idem
  meta.json         formas, espaciado, ángulos, verdad (centro del adenoma en vóxeles de la grilla) y estadísticas
y docs/dicom/<caso>.zip con los DICOM del caso.
"""
from __future__ import annotations

import json
import os
import zipfile

import numpy as np
from scipy import ndimage

from reconstruir import mu_en_grilla

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASOS = ["normal", "inferior-derecho-12", "inferior-izquierdo-8", "retroesofagico-10", "mediastinico-15"]


def grilla():
    """Matriz y píxel de la adquisición SPECT (los de las proyecciones del caso normal)."""
    d = np.load(os.path.join(RAIZ, "salida", "casos", "normal", "proyecciones_precoz.npz"))
    return int(d["matriz"]), float(d["pixel_mm"])


def main():
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    meta_f = json.load(open(os.path.join(RAIZ, "salida", "fantoma.json"), encoding="utf-8"))
    hu, iso = f["hu"].astype(np.float32), float(f["iso"])
    matriz, pix_mm = grilla()
    # CT en la grilla del SPECT (misma convención que la corrección de atenuación y que el DICOM NM)
    ct_g = mu_en_grilla(hu + 1000.0, iso, matriz, pix_mm) - 1000.0     # el relleno fuera del fantoma queda en -1000 HU
    ocupado = np.nonzero((ct_g > -900).any(axis=(1, 2)))[0]
    z0, z1 = int(ocupado.min()), int(ocupado.max()) + 1
    ct8 = np.clip((ct_g[z0:z1] - (-160.0)) / 400.0 * 255.0, 0, 255).astype(np.uint8)
    nz = z1 - z0
    # transformación fantoma(vóxel 2 mm) -> grilla SPECT (vóxel 3.3 mm): el fantoma zoomeado se centra en 128³
    esc = iso / pix_mm
    forma_z = [int(round(s * esc)) for s in hu.shape]
    desplaz = [(matriz - s) // 2 for s in forma_z]       # (z, y, x) índice de inicio en la grilla

    def a_grilla(c_zyx):
        return [c_zyx[0] * esc + desplaz[0] - z0, c_zyx[1] * esc + desplaz[1], c_zyx[2] * esc + desplaz[2]]

    indice = {"casos": [], "nz": nz, "z0": z0, "matriz": matriz, "pixel_mm": pix_mm, "paciente_ct": "TCIA, caso 2 de la entrega docente PET/CT (anónimo)"}
    os.makedirs(os.path.join(RAIZ, "docs", "dicom"), exist_ok=True)
    for caso in CASOS:
        carpeta = os.path.join(RAIZ, "salida", "casos", caso)
        destino = os.path.join(RAIZ, "docs", "datos", caso)
        os.makedirs(destino, exist_ok=True)
        verdad = json.load(open(os.path.join(carpeta, "verdad.json"), encoding="utf-8"))
        adq = json.load(open(os.path.join(carpeta, "adquisicion.json"), encoding="utf-8"))
        meta = {"caso": caso, "nz": nz, "matriz": matriz, "pixel_mm": pix_mm, "fases": {}, "adenoma": None,
                "camara": adq["camara"], "historias": adq["historias"]}
        for fase in ("precoz", "tardia"):
            vol = np.load(os.path.join(carpeta, f"recon_{fase}_ac.npy"))[z0:z1]
            tope = float(np.percentile(vol, 99.95)) or 1.0
            np.clip(vol / tope * 255.0, 0, 255).astype(np.uint8).tofile(os.path.join(destino, f"nm_{fase}.bin"))
            d = np.load(os.path.join(carpeta, f"proyecciones_{fase}.npz"))
            pr = d["ruido"].astype(np.float32)
            tp = float(np.percentile(pr, 99.8)) or 1.0
            np.clip(pr / tp * 255.0, 0, 255).astype(np.uint8).tofile(os.path.join(destino, f"proy_{fase}.bin"))
            ruta_p = os.path.join(carpeta, f"planar_{fase}.npz")
            planar = None
            if os.path.exists(ruta_p):
                pl = np.load(ruta_p)["ruido"].astype(np.float32)
                np.clip(pl / (float(np.percentile(pl, 99.7)) or 1.0) * 255.0, 0, 255).astype(np.uint8).tofile(os.path.join(destino, f"planar_{fase}.bin"))
                planar = {"matriz": int(pl.shape[0]), "pixel_mm": float(np.load(ruta_p)["pixel_mm"]), "cuentas": int(pl.sum()),
                          "tiempo_s": float(np.load(ruta_p)["tiempo_s"]), "max_pixel": int(pl.max())}
            meta["fases"][fase] = {"planar": planar, "n_proy": int(pr.shape[0]), "angulos_grados": [round(float(np.degrees(a)), 1) for a in d["angulos_rad"]],
                                   "cuentas_por_proyeccion": int(round(pr.sum(axis=(1, 2)).mean())), "max_pixel": int(pr.max()),
                                   "actividad_MBq": round(adq["fases"][fase]["actividad_MBq"], 1), "segundos_mc": adq["fases"][fase]["segundos"]}
        if verdad["adenoma"]:
            a = verdad["adenoma"]
            meta["adenoma"] = {"sitio": a["sitio"], "diametro_mm": a["diametro_mm"], "relacion": a["relacion"], "volumen_ml": round(a["volumen_ml"], 2),
                               "centro_grilla_zyx": [round(v, 2) for v in a_grilla(a["centro_voxel_zyx"])]}
            # posición en la planar anterior (misma geometría que el Monte Carlo: u = x - centro, v = z - centro)
            c = a["centro_voxel_zyx"]
            pl_pix = 2.4
            u = ((c[2] + 0.5) * iso - hu.shape[2] * iso / 2.0) / pl_pix + 128
            v = ((c[0] + 0.5) * iso - hu.shape[0] * iso / 2.0) / pl_pix + 128
            meta["adenoma"]["planar_uv"] = [round(u, 1), round(v, 1)]
        json.dump(meta, open(os.path.join(destino, "meta.json"), "w", encoding="utf-8"), ensure_ascii=False)
        # zip de los DICOM
        dic = os.path.join(RAIZ, "salida", "dicom", caso)
        with zipfile.ZipFile(os.path.join(RAIZ, "docs", "dicom", f"{caso}.zip"), "w", zipfile.ZIP_DEFLATED) as z:
            for raiz, _, archivos in os.walk(dic):
                for arch in archivos:
                    p = os.path.join(raiz, arch)
                    z.write(p, os.path.relpath(p, dic))
        indice["casos"].append({"caso": caso, "adenoma": meta["adenoma"] is not None, "zip_mb": round(os.path.getsize(os.path.join(RAIZ, "docs", "dicom", f"{caso}.zip")) / 1e6, 1)})
        print(f"{caso}: nz {nz}, zip {indice['casos'][-1]['zip_mb']} MB, adenoma {meta['adenoma']}")
    json.dump(indice, open(os.path.join(RAIZ, "docs", "datos", "indice.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
