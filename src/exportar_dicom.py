"""Paso 7: DICOM de salida con una sola identidad de estudio, para el visor_dicom y para MicroDicom.

Por caso y fase se escriben en salida/dicom/<caso>/:
  CT/           serie CT recortada y remuestreada (el fantoma), un archivo por corte
  NM_<fase>.dcm volumen reconstruido (ImageType RECON TOMO, multiframe con Slice Vector), en la misma
                grilla LPS que el CT, con FrameOfReferenceUID y PatientID compartidos (fusión alineada)
  PROY_<fase>.dcm proyecciones crudas (ImageType TOMO) con la geometría de la órbita

Reglas que exige core.js del visor (transcritas en sim_mdp/paso04_verificar.py): 16 bits sin signo,
MONOCHROME2, un detector con ImagePositionPatient/ImageOrientationPatient en (0054,0022), una ventana
de energía, FrameIncrementPointer = SliceVector (0054,0080) consecutivo desde 1, SpacingBetweenSlices,
PixelSpacing en la raíz, Units = CNTS, espaciado Z uniforme en el CT, y PatientID + IssuerOfPatientID
iguales entre CT y NM. Paciente sintético: SIM-PARA-<n>, sin datos del original de TCIA salvo la
referencia a la colección en (0012,0063).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid
from scipy import ndimage

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAIZ_UID = "1.2.826.0.1.3680043.10.1245."      # raíz de UID propia (prefijo 1.2.826.0.1.3680043.10 es de uso libre)


def _base(ds: Dataset, meta: dict, caso: str, modalidad: str, uids: dict, serie_num: int, serie_desc: str):
    hoy = dt.datetime(2026, 10, 7, 9, 0, 0)
    ds.PatientName = f"SIM^PARATIROIDES {caso.upper()}"
    ds.PatientID = uids["paciente"]
    ds.IssuerOfPatientID = "SIM"
    ds.PatientBirthDate = ""
    ds.PatientSex = meta.get("sexo", "") or "O"
    ds.PatientAge = meta.get("edad", "") or ""
    ds.PatientIdentityRemoved = "YES"
    ds.DeidentificationMethod = "Fantoma de CT publico TCIA (PS3.15 AnnexE); actividad simulada"
    ds.add_new((0x0012, 0x0063), "LO", "Simulacion Monte Carlo; CT base TCIA (caso 2 PET/CT docente)")
    ds.StudyInstanceUID = uids["estudio"]
    ds.FrameOfReferenceUID = uids["frame"]
    ds.StudyID = "SIMPARA"
    ds.StudyDate = hoy.strftime("%Y%m%d")
    ds.StudyTime = hoy.strftime("%H%M%S")
    ds.StudyDescription = "SPECT/CT PARATIROIDES SESTAMIBI (SIMULADO)"
    ds.AccessionNumber = ""
    ds.ReferringPhysicianName = ""
    ds.Modality = modalidad
    ds.Manufacturer = "simulador-paratiroides-mc"
    ds.ManufacturerModelName = "Monte Carlo LEHR"
    ds.SeriesInstanceUID = uids[f"serie_{serie_desc}"]
    ds.SeriesNumber = serie_num
    ds.SeriesDescription = serie_desc
    ds.SeriesDate = ds.StudyDate
    ds.SeriesTime = ds.StudyTime
    ds.SpecificCharacterSet = "ISO_IR 100"


def _archivo(sop_class: str, sop_uid: str) -> FileDataset:
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = sop_class
    fm.MediaStorageSOPInstanceUID = sop_uid
    fm.TransferSyntaxUID = ExplicitVRLittleEndian
    fm.ImplementationClassUID = RAIZ_UID + "1"
    ds = FileDataset("", {}, file_meta=fm, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    ds.SOPClassUID = sop_class
    ds.SOPInstanceUID = sop_uid
    return ds


def exportar_ct(hu: np.ndarray, iso: float, origen: list, meta: dict, caso: str, uids: dict, carpeta: str):
    """hu (z,y,x) a iso mm, origen LPS del vóxel (0,0,0). Un archivo por corte, z creciente."""
    os.makedirs(carpeta, exist_ok=True)
    nz = hu.shape[0]
    for k in range(nz):
        ds = _archivo("1.2.840.10008.5.1.4.1.1.2", RAIZ_UID + f"2.{uids['semilla']}.{k + 1}")
        _base(ds, meta, caso, "CT", uids, 2, "CT")
        ds.ImageType = ["DERIVED", "SECONDARY", "AXIAL"]
        ds.InstanceNumber = k + 1
        ds.ImagePositionPatient = [f"{origen[0]:.4f}", f"{origen[1]:.4f}", f"{origen[2] + k * iso:.4f}"]
        ds.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
        ds.SliceLocation = f"{origen[2] + k * iso:.4f}"
        ds.PixelSpacing = [f"{iso:.4f}", f"{iso:.4f}"]
        ds.SliceThickness = f"{iso:.4f}"
        ds.SpacingBetweenSlices = f"{iso:.4f}"
        ds.KVP = "140"
        ds.Rows, ds.Columns = hu.shape[1], hu.shape[2]
        ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
        ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 1
        ds.RescaleIntercept, ds.RescaleSlope, ds.RescaleType = "0", "1", "HU"
        ds.WindowCenter, ds.WindowWidth = "40", "400"
        ds.PixelData = np.ascontiguousarray(hu[k].astype(np.int16)).tobytes()
        ds.save_as(os.path.join(carpeta, f"CT_{k + 1:03d}.dcm"), enforce_file_format=True)
    return nz


def _nm_comun(ds: Dataset, filas: int, cols: int, n: int, pix_mm: float):
    ds.Rows, ds.Columns, ds.NumberOfFrames = filas, cols, n
    ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
    ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 0
    ds.PixelSpacing = [f"{pix_mm:.4f}", f"{pix_mm:.4f}"]
    ds.Units = "CNTS"
    ds.CountsSource = "EMISSION"
    ds.NumberOfEnergyWindows = 1
    ew = Dataset()
    ew.EnergyWindowName = "Tc99m 20%"
    rango = Dataset()
    rango.EnergyWindowLowerLimit, rango.EnergyWindowUpperLimit = "126.45", "154.55"
    ew.EnergyWindowRangeSequence = Sequence([rango])
    ds.EnergyWindowInformationSequence = Sequence([ew])
    rf = Dataset()
    rn = Dataset()
    rn.CodeValue, rn.CodingSchemeDesignator, rn.CodeMeaning = "C-163A8", "SRT", "Technetium^99m"
    rf.RadionuclideCodeSequence = Sequence([rn])
    rf.Radiopharmaceutical = "Tc-99m sestamibi"
    rf.RadionuclideTotalDose = "740000000"
    ds.RadiopharmaceuticalInformationSequence = Sequence([rf])
    ds.NumberOfDetectors = 1


def exportar_nm_recon(vol: np.ndarray, pix_mm: float, origen: list, meta: dict, caso: str, fase: str, uids: dict, ruta: str):
    """vol (z,y,x) en cuentas relativas, grilla del detector (pix_mm isotrópico) con origen LPS dado."""
    n, filas, cols = vol.shape
    escala = 65000.0 / max(float(vol.max()), 1e-6)
    px = np.clip(np.round(vol * escala), 0, 65535).astype(np.uint16)
    ds = _archivo("1.2.840.10008.5.1.4.1.1.20", RAIZ_UID + f"3.{uids['semilla']}.{1 if fase == 'precoz' else 2}")
    _base(ds, meta, caso, "NM", uids, 3 if fase == "precoz" else 4, f"SPECT {fase.upper()} RECON AC")
    ds.ImageType = ["ORIGINAL", "PRIMARY", "RECON TOMO", "EMISSION"]
    ds.InstanceNumber = 1
    _nm_comun(ds, filas, cols, n, pix_mm)
    ds.SliceThickness = f"{pix_mm:.4f}"
    ds.SpacingBetweenSlices = f"{pix_mm:.4f}"
    ds.RescaleIntercept, ds.RescaleSlope = "0", f"{1.0 / escala:.8g}"
    ds.FrameIncrementPointer = (0x0054, 0x0080)
    ds.SliceVector = list(range(1, n + 1))
    ds.NumberOfSlices = n
    det = Dataset()
    det.ImagePositionPatient = [f"{origen[0]:.4f}", f"{origen[1]:.4f}", f"{origen[2]:.4f}"]
    det.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
    det.CollimatorGridName, det.CollimatorType = "LEHR", "PARA"
    ds.DetectorInformationSequence = Sequence([det])
    ds.PixelData = np.ascontiguousarray(px).tobytes()
    ds.save_as(ruta, enforce_file_format=True)
    return escala


def exportar_proyecciones(proy: np.ndarray, pix_mm: float, angulos_rad: np.ndarray, radio_cm: float, tiempo_s: float, meta: dict, caso: str, fase: str, uids: dict, ruta: str):
    n, filas, cols = proy.shape
    ds = _archivo("1.2.840.10008.5.1.4.1.1.20", RAIZ_UID + f"4.{uids['semilla']}.{1 if fase == 'precoz' else 2}")
    _base(ds, meta, caso, "NM", uids, 5 if fase == "precoz" else 6, f"SPECT {fase.upper()} PROYECCIONES")
    ds.ImageType = ["ORIGINAL", "PRIMARY", "TOMO", "EMISSION"]
    ds.InstanceNumber = 1
    _nm_comun(ds, filas, cols, n, pix_mm)
    ds.FrameIncrementPointer = (0x0054, 0x0090)
    ds.AngularViewVector = list(range(1, n + 1))
    ds.TypeOfDetectorMotion = "STEP AND SHOOT"
    rot = Dataset()
    rot.StartAngle = f"{np.degrees(angulos_rad[0]):.2f}"
    rot.AngularStep = f"{np.degrees(angulos_rad[1] - angulos_rad[0]) if n > 1 else 0:.2f}"
    rot.ScanArc = f"{np.degrees(angulos_rad[-1] - angulos_rad[0]) + np.degrees(angulos_rad[1] - angulos_rad[0]) if n > 1 else 0:.1f}"
    rot.RotationDirection = "CW"
    rot.RadialPosition = [f"{radio_cm * 10:.1f}"] * n
    rot.NumberOfFramesInRotation = n
    rot.ActualFrameDuration = str(int(tiempo_s * 1000))
    ds.RotationInformationSequence = Sequence([rot])
    ds.NumberOfRotations = 1
    det = Dataset()
    det.CollimatorGridName, det.CollimatorType = "LEHR", "PARA"
    det.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
    det.ImagePositionPatient = ["0", "0", "0"]
    ds.DetectorInformationSequence = Sequence([det])
    ds.RescaleIntercept, ds.RescaleSlope = "0", "1"
    ds.PixelData = np.ascontiguousarray(np.clip(proy, 0, 65535).astype(np.uint16)).tobytes()
    ds.save_as(ruta, enforce_file_format=True)


def volumen_a_lps(vol: np.ndarray, pix_mm: float, fantoma_meta: dict, fantoma_forma, iso: float):
    """El volumen reconstruido (matriz³ a pix_mm, centrado en el centro del fantoma) -> origen LPS de su vóxel (0,0,0).
    Convención del reconstructor: eje y del volumen = eje u del detector; con el detector en phi=0 mirando desde +x."""
    o = fantoma_meta["origen_mm"]
    centro = [o[0] + fantoma_forma[2] * iso / 2.0, o[1] + fantoma_forma[1] * iso / 2.0, o[2] + fantoma_forma[0] * iso / 2.0]
    n = vol.shape
    return [centro[0] - n[2] * pix_mm / 2.0, centro[1] - n[1] * pix_mm / 2.0, centro[2] - n[0] * pix_mm / 2.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", default="inferior-derecho-12")
    ap.add_argument("--numero", type=int, default=1, help="número de paciente sintético SIM-PARA-<n>")
    ap.add_argument("--fases", default="precoz,tardia")
    a = ap.parse_args()
    f = np.load(os.path.join(RAIZ, "salida", "fantoma.npz"))
    meta = json.load(open(os.path.join(RAIZ, "salida", "fantoma.json"), encoding="utf-8"))
    hu, iso = f["hu"], float(f["iso"])
    carpeta_caso = os.path.join(RAIZ, "salida", "casos", a.caso)
    salida = os.path.join(RAIZ, "salida", "dicom", a.caso)
    os.makedirs(salida, exist_ok=True)
    uids = {"paciente": f"SIM-PARA-{a.numero:02d}", "estudio": RAIZ_UID + f"1.{a.numero}.1", "frame": RAIZ_UID + f"1.{a.numero}.2",
            "serie_CT": RAIZ_UID + f"1.{a.numero}.3", "semilla": a.numero}
    for fase in ("precoz", "tardia"):
        uids[f"serie_SPECT {fase.upper()} RECON AC"] = RAIZ_UID + f"1.{a.numero}.{4 if fase == 'precoz' else 5}"
        uids[f"serie_SPECT {fase.upper()} PROYECCIONES"] = RAIZ_UID + f"1.{a.numero}.{6 if fase == 'precoz' else 7}"
    n_ct = exportar_ct(hu, iso, meta["origen_mm"], meta, a.caso, uids, os.path.join(salida, "CT"))
    resumen = {"paciente": uids["paciente"], "ct_cortes": n_ct, "fases": {}}
    for fase in a.fases.split(","):
        d = np.load(os.path.join(carpeta_caso, f"proyecciones_{fase}.npz"))
        pix_mm = float(d["pixel_mm"])
        exportar_proyecciones(d["ruido"], pix_mm, d["angulos_rad"], float(d["radio_cm"]), 25.0, meta, a.caso, fase, uids, os.path.join(salida, f"PROY_{fase}.dcm"))
        ruta_rec = os.path.join(carpeta_caso, f"recon_{fase}_ac.npy")
        if os.path.exists(ruta_rec):
            vol = np.load(ruta_rec)
            origen = volumen_a_lps(vol, pix_mm, meta, hu.shape, iso)
            esc = exportar_nm_recon(vol, pix_mm, origen, meta, a.caso, fase, uids, os.path.join(salida, f"NM_{fase}.dcm"))
            resumen["fases"][fase] = {"recon": True, "escala_cuentas": esc, "origen_lps": origen}
        else:
            resumen["fases"][fase] = {"recon": False}
    json.dump(resumen, open(os.path.join(salida, "resumen.json"), "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print(f"{a.caso}: CT {n_ct} cortes; fases {resumen['fases']} -> {salida}")


if __name__ == "__main__":
    main()
