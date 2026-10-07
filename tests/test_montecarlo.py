"""Validación del motor contra lo que se sabe de una cámara LEHR y de la física del Tc-99m."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import montecarlo as mc  # noqa: E402


def _fuente_puntual(mu_valor=0.0, lado=21, vox=1.0):
    """Fuente puntual en el centro de un cubo (lado vóxeles de `vox` cm) con mu uniforme."""
    act = np.zeros((lado, lado, lado), np.float32)
    act[lado // 2, lado // 2, lado // 2] = 1000.0
    mu = np.full((lado, lado, lado), mu_valor, np.float32)
    hu = np.zeros((lado, lado, lado), np.int16)
    return act, mu, hu


def test_klein_nishina_normalizada():
    c = np.linspace(-1, 1, 20001)
    for E in (140.5, 100.0):
        p = np.array([mc._kn_pdf(E, ci) for ci in c])
        integral = 2 * math.pi * np.trapezoid(p, c)
        assert abs(integral - 1.0) < 1e-3, integral


def test_sensibilidad_fuente_puntual_en_aire():
    act, mu, hu = _fuente_puntual(0.0)
    cam = mc.Camara(n_angulos=4, radio_cm=10.0, matriz=64, pix_mm=3.3, tiempo_s=1.0)
    esp, _, est = mc.simular(act, mu, hu, 1.0, cam, n_hist=200_000, semilla=3, poisson=False)
    por_foton = esp.sum(axis=(1, 2)) / est["emitidos_por_proyeccion"]
    teorica = mc.G_GEOM * mc.EFIC_CRISTAL * mc._p_ventana(mc.E0)
    assert np.allclose(por_foton, teorica, rtol=0.02), (por_foton, teorica)


def test_resolucion_a_10_cm():
    act, mu, hu = _fuente_puntual(0.0, lado=21, vox=0.1)          # fuente de 1 mm, no de 1 cm
    cam = mc.Camara(n_angulos=1, radio_cm=10.0, matriz=128, pix_mm=1.0, tiempo_s=1.0)
    esp, _, _ = mc.simular(act, mu, hu, 0.1, cam, n_hist=50_000, semilla=1, poisson=False)
    perfil = esp[0].sum(axis=0)
    perfil /= perfil.max()
    arriba = np.nonzero(perfil >= 0.5)[0]
    fwhm = (arriba.max() - arriba.min() + 1) * cam.pix_mm
    fwhm_col = mc.D_AGUJERO * (mc.L_EFF + 100.0) / mc.L_EFF
    esperada = math.hypot(fwhm_col, mc.FWHM_INTRINSECA)        # ~11.3 mm
    assert abs(fwhm - esperada) <= 1.5, (fwhm, esperada)


def test_atenuacion_y_dispersion_en_agua():
    """Fuente en el centro de un cubo de agua de 20 cm: los primarios se atenúan exp(-mu·10) y la dispersión agrega 20-50 %."""
    act, mu, hu = _fuente_puntual(0.1537, lado=21, vox=1.0)
    cam = mc.Camara(n_angulos=4, radio_cm=15.0, matriz=64, pix_mm=3.3, tiempo_s=1.0)
    esp, _, est = mc.simular(act, mu, hu, 1.0, cam, n_hist=300_000, semilla=5, poisson=False)
    por_foton = esp.sum(axis=(1, 2)).mean() / est["emitidos_por_proyeccion"]
    primario = mc.G_GEOM * mc.EFIC_CRISTAL * mc._p_ventana(mc.E0) * math.exp(-0.1537 * 10.5)
    assert 1.15 * primario < por_foton < 1.6 * primario, (por_foton / primario)
    assert est["compton"] > est["fotoelectricos"] * 5       # a 140 keV en agua domina Compton


def test_reproducible():
    act, mu, hu = _fuente_puntual(0.1537, lado=11)
    cam = mc.Camara(n_angulos=2, radio_cm=10.0, matriz=32, pix_mm=3.3, tiempo_s=1.0)
    a, _, _ = mc.simular(act, mu, hu, 1.0, cam, n_hist=20_000, semilla=7, poisson=False, n_threads=2)
    b, _, _ = mc.simular(act, mu, hu, 1.0, cam, n_hist=20_000, semilla=7, poisson=False, n_threads=2)
    assert np.array_equal(a, b)
