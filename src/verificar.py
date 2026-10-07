"""Paso 8: verificación de los DICOM exportados contra las reglas de core.js del visor_dicom y de la fusión.

Uso: python src/verificar.py --caso inferior-derecho-12
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pydicom

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fallas = []


def check(ok, texto, detalle=""):
    print(f"  [{'ok' if ok else 'FALLA'}] {texto}{('  -> ' + str(detalle)) if detalle else ''}")
    if not ok:
        fallas.append(texto)
    return ok


def verificar_nm(ruta):
    ds = pydicom.dcmread(ruta)
    print(f"\n== NM {os.path.basename(ruta)}")
    n = int(ds.NumberOfFrames)
    check(str(ds.file_meta.TransferSyntaxUID) in ("1.2.840.10008.1.2", "1.2.840.10008.1.2.1", "1.2.840.10008.1.2.2"), "transfer syntax sin comprimir")
    check(ds.Modality == "NM", "modalidad NM")
    check(list(ds.ImageType)[2] == "RECON TOMO", "ImageType[2] = RECON TOMO", "\\".join(ds.ImageType))
    check(int(ds.SamplesPerPixel) == 1 and ds.PhotometricInterpretation == "MONOCHROME2", "MONOCHROME2, 1 muestra")
    check(int(ds.BitsAllocated) == 16 and int(ds.BitsStored) == 16 and int(ds.HighBit) == 15 and int(ds.PixelRepresentation) == 0, "16 bits sin signo")
    det = ds.DetectorInformationSequence
    check(len(det) == 1, "un único detector en (0054,0022)")
    ipp = [float(v) for v in det[0].ImagePositionPatient]
    iop = [float(v) for v in det[0].ImageOrientationPatient]
    check(len(ipp) == 3 and len(iop) == 6 and all(abs(a - b) <= 1e-4 for a, b in zip(iop, [1, 0, 0, 0, 1, 0])), "geometría axial LPS en el detector")
    check(float(ds.RescaleSlope) != 0, "rescale válido", f"slope {ds.RescaleSlope}")
    fip = ds["FrameIncrementPointer"]
    check(pydicom.tag.Tag(fip.value) == pydicom.tag.Tag(0x0054, 0x0080), "FrameIncrementPointer = SliceVector")
    check(int(ds.NumberOfEnergyWindows) == 1 and int(ds.NumberOfDetectors) == 1, "una ventana y un detector")
    sv = list(ds.SliceVector)
    check(len(sv) == n and int(ds.NumberOfSlices) == n and all(v == i + 1 for i, v in enumerate(sv)), "SliceVector consecutivo desde 1")
    check(abs(float(ds.SpacingBetweenSlices)) >= 1e-4, "SpacingBetweenSlices válido", ds.SpacingBetweenSlices)
    check(len(ds.PixelData) >= int(ds.Rows) * int(ds.Columns) * n * 2, "píxeles completos")
    check(ds.Units == "CNTS", "Units = CNTS")
    ps = [float(v) for v in ds.PixelSpacing]
    return ds, ipp, [ps[1], ps[0], float(ds.SpacingBetweenSlices)], (n, int(ds.Rows), int(ds.Columns))


def verificar_ct(carpeta):
    filas = []
    for nombre in sorted(os.listdir(carpeta)):
        ds = pydicom.dcmread(os.path.join(carpeta, nombre))
        filas.append((float(ds.ImagePositionPatient[2]), ds))
    filas.sort(key=lambda r: r[0])
    f0 = filas[0][1]
    print(f"\n== CT ({len(filas)} cortes)")
    ok_geo = all(ds.SeriesInstanceUID == f0.SeriesInstanceUID and ds.FrameOfReferenceUID == f0.FrameOfReferenceUID and ds.PatientID == f0.PatientID
                 and int(ds.Rows) == int(f0.Rows) and int(ds.Columns) == int(f0.Columns) for _, ds in filas)
    check(ok_geo, "serie homogénea")
    zs = np.array([z for z, _ in filas])
    d = np.diff(zs)
    dz = float(np.median(d))
    check(dz >= 1e-4 and float(np.max(np.abs(d - dz))) <= max(0.01, dz * 0.01), "espaciado Z uniforme", f"dz={dz:.3f}")
    check(int(f0.Rows) * int(f0.Columns) * len(filas) <= 128 * 1024 * 1024, "<= 128 M vóxeles")
    texto = str(f0).upper()
    check("PET-02" not in texto and "MAGDALENA" not in texto, "sin ID ni ruta del CT original en el exportado")
    ipp = [float(v) for v in f0.ImagePositionPatient]
    ps = [float(v) for v in f0.PixelSpacing]
    return f0, [ipp[0], ipp[1], zs[0]], [ps[1], ps[0], dz], (len(filas), int(f0.Rows), int(f0.Columns))


def verificar_fusion(ct, nm, ct_o, ct_e, ct_n, nm_o, nm_e, nm_n):
    print("\n== Fusión")
    check(ct.PatientID == nm.PatientID, "PatientID compartido", ct.PatientID)
    check(str(getattr(ct, "IssuerOfPatientID", "")) == str(getattr(nm, "IssuerOfPatientID", "")), "IssuerOfPatientID compartido")
    check(ct.FrameOfReferenceUID == nm.FrameOfReferenceUID, "FrameOfReferenceUID compartido")
    check(ct.StudyInstanceUID == nm.StudyInstanceUID, "mismo estudio")
    check(str(ct.PatientSex) == str(nm.PatientSex) and str(getattr(ct, "PatientAge", "")) == str(getattr(nm, "PatientAge", "")), "sexo y edad iguales (MicroDicom)")
    n_ct = [ct_n[2], ct_n[1], ct_n[0]]
    n_nm = [nm_n[2], nm_n[1], nm_n[0]]
    for i, eje in enumerate("xyz"):
        a0, a1 = ct_o[i], ct_o[i] + (n_ct[i] - 1) * ct_e[i]
        b0, b1 = nm_o[i], nm_o[i] + (n_nm[i] - 1) * nm_e[i]
        print(f"      {eje}: CT {a0:8.1f}..{a1:8.1f}   NM {b0:8.1f}..{b1:8.1f}   solape {max(a0, b0):8.1f}..{min(a1, b1):8.1f} mm")
    check(all(max(ct_o[i], nm_o[i]) <= min(ct_o[i] + (n_ct[i] - 1) * ct_e[i], nm_o[i] + (n_nm[i] - 1) * nm_e[i]) for i in range(3)), "los volúmenes se solapan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="inferior-derecho-12")
    a = ap.parse_args()
    carpeta = os.path.join(RAIZ, "salida", "dicom", a.caso)
    ct, ct_o, ct_e, ct_n = verificar_ct(os.path.join(carpeta, "CT"))
    for fase in ("precoz", "tardia"):
        ruta = os.path.join(carpeta, f"NM_{fase}.dcm")
        if os.path.exists(ruta):
            nm, nm_o, nm_e, nm_n = verificar_nm(ruta)
            verificar_fusion(ct, nm, ct_o, ct_e, ct_n, nm_o, nm_e, nm_n)
    print(f"\n{'SIN FALLAS' if not fallas else str(len(fallas)) + ' FALLAS: ' + '; '.join(fallas)}")


if __name__ == "__main__":
    main()
