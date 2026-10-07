"""Paso 2a: TotalSegmentator sobre el fantoma (CPU) y etiquetas en la grilla del fantoma.

Escribe el fantoma como NIfTI con su geometría exacta (la salida cae en la misma grilla), corre
TotalSegmentator con el subconjunto de estructuras de cuello y tórax, y guarda las etiquetas en
(z, y, x) en salida/totalseg/etiquetas.npz.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

import nibabel as nib
import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ESTRUCTURAS = ["thyroid_gland", "esophagus", "trachea", "heart", "aorta", "pulmonary_vein", "superior_vena_cava",
               "liver", "spleen", "stomach",
               "lung_upper_lobe_left", "lung_upper_lobe_right", "lung_middle_lobe_right", "lung_lower_lobe_left", "lung_lower_lobe_right",
               "common_carotid_artery_left", "common_carotid_artery_right", "brachiocephalic_trunk", "subclavian_artery_left",
               "subclavian_artery_right", "brachiocephalic_vein_left", "brachiocephalic_vein_right",
               "sternum", "clavicula_left", "clavicula_right", "scapula_left", "scapula_right", "humerus_left", "humerus_right"] + \
              [f"vertebrae_C{i}" for i in range(2, 8)] + [f"vertebrae_T{i}" for i in range(1, 13)] + \
              [f"rib_{l}_{i}" for l in ("left", "right") for i in range(1, 13)]


def main():
    d = os.path.join(RAIZ, "salida", "totalseg")
    os.makedirs(d, exist_ok=True)
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    m = json.load(open(os.path.join(RAIZ, "salida", "fantoma.json"), encoding="utf-8"))
    hu, iso, o = f["hu"].astype(np.int16), float(f["iso"]), m["origen_mm"]
    aff = np.array([[-iso, 0, 0, -o[0]], [0, -iso, 0, -o[1]], [0, 0, iso, o[2]], [0, 0, 0, 1]], float)
    nib.save(nib.Nifti1Image(np.transpose(hu, (2, 1, 0)), aff), os.path.join(d, "cuello_ct.nii.gz"))
    exe = os.path.join(os.path.dirname(sys.executable), "TotalSegmentator")
    r = subprocess.run([exe, "-i", os.path.join(d, "cuello_ct.nii.gz"), "-o", os.path.join(d, "total.nii.gz"), "--ml", "-d", "cpu",
                        "--roi_subset"] + ESTRUCTURAS, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
        sys.exit(1)
    et = np.transpose(np.asarray(nib.load(os.path.join(d, "total.nii.gz")).dataobj).astype(np.uint8), (2, 1, 0))
    assert et.shape == hu.shape
    np.savez_compressed(os.path.join(d, "etiquetas.npz"), et=et)
    from totalsegmentator.map_to_binary import class_map
    nombres = class_map["total"]
    vol = {nombres[int(i)]: round(float((et == i).sum()) * iso ** 3 / 1000.0, 1) for i in np.unique(et) if i}
    print("TotalSegmentator (mL):", {k: v for k, v in vol.items() if not k.startswith(("rib_", "vertebrae_"))})


if __name__ == "__main__":
    main()
