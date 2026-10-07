"""Paso 4: Monte Carlo de fotones para SPECT con colimador paralelo (CPU, Numba, en paralelo).

Cada historia es un fotón de 140 keV emitido desde el mapa de actividad. Viaja por el fantoma con el
método de Woodcock (tracking delta: paso con el mu máximo, interacción real con probabilidad
mu(x,E)/mu_max). Interacciones: fotoeléctrico (absorción) y Compton (Klein-Nishina, método de Kahn);
Rayleigh se omite (menos del 3 % a estas energías y no cambia la energía). Los fotones bajo 90 keV
se descartan: no pueden entrar en la ventana de 126-154 keV ni con la resolución del detector.

Detección forzada con múltiples ángulos: en la emisión y en cada vértice Compton, para cada ángulo
de proyección se suma al detector la probabilidad de que el fotón saliera hacia el colimador
(isotrópico: g; tras Compton: p_KN(theta)·4pi·g), multiplicada por la transmisión del rayo hasta el
colimador, la eficiencia del cristal y la probabilidad de caer en la ventana con la resolución
energética; el peso se reparte en el detector con la PSF del colimador a esa distancia más la
intrínseca. Después el fotón sigue su camino real. Es la "detección forzada por convolución" de
SIMIND, con todos los ángulos reutilizando cada historia.

Las cuentas esperadas por píxel se escalan a los fotones realmente emitidos (A·0.889·t) y se
muestrean con Poisson. Semilla explícita y generador propio (xorshift) por historia: reproducible.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

# ---- datos físicos ----
E0 = 140.5                        # keV, Tc-99m
RENDIMIENTO_GAMMA = 0.889         # fotones de 140 keV por desintegración
ME = 511.0                        # keV
# NIST XCOM (cm²/g): energía, incoherente (Compton), fotoeléctrico; agua y hueso cortical ICRU-44
_E = np.array([40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0])
_AGUA_C = np.array([0.1830, 0.1803, 0.1757, 0.1658, 0.1565, 0.1392, 0.1264])
_AGUA_F = np.array([0.0646, 0.0310, 0.0183, 0.0074, 0.0037, 0.0010, 0.0004])
_HUESO_C = np.array([0.1680, 0.1660, 0.1610, 0.1520, 0.1440, 0.1280, 0.1160])
_HUESO_F = np.array([0.4830, 0.2340, 0.1310, 0.0531, 0.0260, 0.0073, 0.0029])
RHO_AGUA, RHO_HUESO = 1.0, 1.92
HU_HUESO = 200                    # vóxeles con HU > 200 usan la tabla de hueso

# ---- cámara (Siemens Symbia LEHR, cristal 9.5 mm NaI) ----
D_AGUJERO, L_AGUJERO, SEPTO = 1.11, 24.05, 0.16          # mm
MU_PB_140 = 2.70                                        # 1/mm
L_EFF = L_AGUJERO - 2.0 / MU_PB_140
G_GEOM = 0.26 ** 2 * (D_AGUJERO / L_EFF) ** 2 * (D_AGUJERO / (D_AGUJERO + SEPTO)) ** 2   # ~1.17e-4
FWHM_INTRINSECA = 3.8                                   # mm
RES_ENERGIA_140 = 0.095                                 # FWHM relativa a 140 keV
EFIC_CRISTAL = 0.88
VENTANA = (126.45, 154.55)                              # keV (ventana del 20 %)


def tablas(n=1024, e_min=30.0, e_max=160.0):
    """Tablas finas en energía: razón mu_total(E)/mu_total(140) y fracción Compton para agua y hueso; mu_max por energía."""
    E = np.linspace(e_min, e_max, n)
    def interp(y):
        return np.exp(np.interp(np.log(E), np.log(_E), np.log(y)))
    ac, af, hc, hf = interp(_AGUA_C), interp(_AGUA_F), interp(_HUESO_C), interp(_HUESO_F)
    a_tot, h_tot = (ac + af) * RHO_AGUA, (hc + hf) * RHO_HUESO
    a140 = np.interp(E0, E, a_tot)
    h140 = np.interp(E0, E, h_tot)
    return E, (a_tot / a140).astype(np.float64), (ac / (ac + af)).astype(np.float64), (h_tot / h140).astype(np.float64), (hc / (hc + hf)).astype(np.float64)


@njit(cache=True, inline="always")
def _xorshift(s):
    s ^= (s << np.uint64(13)) & np.uint64(0xFFFFFFFFFFFFFFFF)
    s ^= s >> np.uint64(7)
    s ^= (s << np.uint64(17)) & np.uint64(0xFFFFFFFFFFFFFFFF)
    return s


@njit(cache=True, inline="always")
def _u01(s):
    s = _xorshift(s)
    return s, (s >> np.uint64(11)) * (1.0 / 9007199254740992.0)


@njit(cache=True)
def _indice_e(E, e_min, de, n):
    i = int((E - e_min) / de)
    if i < 0:
        i = 0
    elif i >= n:
        i = n - 1
    return i


@njit(cache=True)
def _kahn(s, E):
    """Muestreo de Klein-Nishina (Kahn). Devuelve (estado, cos(theta), E')."""
    a = E / ME
    while True:
        s, r1 = _u01(s)
        s, r2 = _u01(s)
        s, r3 = _u01(s)
        if r1 <= (2.0 * a + 1.0) / (2.0 * a + 9.0):
            x = 1.0 + 2.0 * a * r2
            if r3 <= 4.0 * (1.0 / x - 1.0 / (x * x)):
                c = 1.0 - 2.0 * r2
                return s, c, E / x
        else:
            x = (2.0 * a + 1.0) / (1.0 + 2.0 * a * r2)
            c = 1.0 - (x - 1.0) / a
            if r3 <= 0.5 * (c * c + 1.0 / x):
                return s, c, E / x


@njit(cache=True, inline="always")
def _kn_pdf(E, c):
    """Densidad de probabilidad de Klein-Nishina por estereorradián, normalizada a 1 sobre 4pi."""
    a = E / ME
    x = 1.0 / (1.0 + a * (1.0 - c))
    f = x * x * (x + 1.0 / x - (1.0 - c * c))
    # normalización: integral de f sobre dOmega = sigma_KN / (r_e^2/2); aproximación numérica en serie cerrada
    t = 1.0 + 2.0 * a
    sig = 2.0 * math.pi * ((1.0 + a) / (a * a) * (2.0 * (1.0 + a) / t - math.log(t) / a) + math.log(t) / (2.0 * a) - (1.0 + 3.0 * a) / (t * t))
    return 0.5 * f / sig


@njit(cache=True)
def _transmision(mu, nz, ny, nx, vox, x0, y0, z0, dx, dy, dz, razon_e, paso):
    """Integral de mu a lo largo del rayo desde (x0,y0,z0) en dirección (dx,dy,dz) hasta salir de la grilla; razon_e escala mu a la energía."""
    tau = 0.0
    x, y, z = x0, y0, z0
    for _ in range(4000):
        i, j, k = int(z / vox), int(y / vox), int(x / vox)
        if i < 0 or i >= nz or j < 0 or j >= ny or k < 0 or k >= nx:
            break
        tau += mu[i, j, k] * paso
        x += dx * paso
        y += dy * paso
        z += dz * paso
    return math.exp(-tau * razon_e)


@njit(cache=True, inline="always")
def _p_ventana(E):
    sigma = RES_ENERGIA_140 * E0 * math.sqrt(E / E0) / 2.3548
    lo = (VENTANA[0] - E) / (sigma * 1.4142135623730951)
    hi = (VENTANA[1] - E) / (sigma * 1.4142135623730951)
    return 0.5 * (math.erf(hi) - math.erf(lo))


@njit(cache=True)
def _depositar(acum, t, ia, u, v, w, dist_mm, pix_mm, npx, npy):
    """Reparte el peso w en el detector con la PSF (colimador a distancia dist_mm + intrínseca)."""
    if dist_mm < 0.0:
        dist_mm = 0.0
    fwhm_c = D_AGUJERO * (L_EFF + dist_mm) / L_EFF
    fwhm = math.sqrt(fwhm_c * fwhm_c + FWHM_INTRINSECA * FWHM_INTRINSECA)
    sig = fwhm / 2.3548 / pix_mm
    cu, cv = u / pix_mm + npx * 0.5, v / pix_mm + npy * 0.5
    r = int(3.0 * sig) + 1
    iu0, iv0 = int(cu), int(cv)
    norm = 0.0
    for a in range(iu0 - r, iu0 + r + 1):
        for b in range(iv0 - r, iv0 + r + 1):
            g = math.exp(-0.5 * (((a + 0.5 - cu) / sig) ** 2 + ((b + 0.5 - cv) / sig) ** 2))
            norm += g
    if norm <= 0.0:
        return
    for a in range(iu0 - r, iu0 + r + 1):
        if a < 0 or a >= npx:
            continue
        for b in range(iv0 - r, iv0 + r + 1):
            if b < 0 or b >= npy:
                continue
            g = math.exp(-0.5 * (((a + 0.5 - cu) / sig) ** 2 + ((b + 0.5 - cv) / sig) ** 2))
            acum[t, ia, b, a] += w * g / norm


@njit(cache=True)
def _forzar(acum, t, mu, nz, ny, nx, vox, x, y, z, dirx, diry, dirz, tiene_dir, E, w, angulos, radio_cm, cz,
            pix_mm, npx, npy, E_tab, e_min, de, ne, r_agua, r_hueso, es_hueso_max, paso):
    """Contribución forzada a todos los ángulos desde el vértice (x,y,z) en cm (coordenadas de grilla)."""
    cxg, cyg = nx * vox * 0.5, ny * vox * 0.5
    for ia in range(angulos.shape[0]):
        phi = angulos[ia]
        # normal del detector apuntando del centro hacia el detector
        n_x, n_y = math.cos(phi), math.sin(phi)
        if tiene_dir:
            c = dirx * n_x + diry * n_y
            p = _kn_pdf(E, c) * 4.0 * math.pi * G_GEOM
            E_sal = E / (1.0 + (E / ME) * (1.0 - c))     # energía con la que viaja hacia ESE detector
        else:
            p = G_GEOM
            E_sal = E
        pw = _p_ventana(E_sal) * EFIC_CRISTAL
        if pw <= 1e-9:
            continue
        razon = r_agua[_indice_e(E_sal, e_min, de, ne)]
        # distancia del vértice al plano del colimador (a radio_cm del centro de giro)
        d_plano = radio_cm - ((x - cxg) * n_x + (y - cyg) * n_y)
        if d_plano < 0.0:
            continue
        T = _transmision(mu, nz, ny, nx, vox, x, y, z, n_x, n_y, 0.0, razon, paso)
        # coordenadas en el detector: u tangencial (eje -sin, cos), v axial
        u = (-(x - cxg) * n_y + (y - cyg) * n_x) * 10.0
        v = (z - cz) * 10.0
        _depositar(acum, t, ia, u, v, w * p * T * pw, d_plano * 10.0, pix_mm, npx, npy)


@njit(cache=True, parallel=True)
def _simular(n_hist, semilla, cdf, mu, es_hueso, vox, angulos, radio_cm, cz, pix_mm, npx, npy, acum, E_tab, e_min, de, ne,
             r_agua, r_hueso, f_agua, f_hueso, mu_max_tab, paso, n_threads, estad):
    nz, ny, nx = mu.shape
    nvox = nz * ny * nx
    por_hilo = (n_hist + n_threads - 1) // n_threads
    for t in prange(n_threads):
        s = np.uint64(0x9E3779B97F4A7C15) ^ (np.uint64(semilla) * np.uint64(0xBF58476D1CE4E5B9) + np.uint64(t + 1) * np.uint64(0x94D049BB133111EB))
        for _ in range(8):
            s = _xorshift(s)
        n_mio = min(por_hilo, n_hist - t * por_hilo)
        for h in range(n_mio):
            # emisión: vóxel según la actividad, posición uniforme dentro
            s, r = _u01(s)
            lo, hi = 0, nvox - 1
            while lo < hi:
                m = (lo + hi) // 2
                if cdf[m] < r:
                    lo = m + 1
                else:
                    hi = m
            i = lo // (ny * nx)
            j = (lo // nx) % ny
            k = lo % nx
            s, r1 = _u01(s)
            s, r2 = _u01(s)
            s, r3 = _u01(s)
            z, y, x = (i + r1) * vox, (j + r2) * vox, (k + r3) * vox
            E = E0
            w = 1.0
            # dirección isotrópica
            s, r1 = _u01(s)
            s, r2 = _u01(s)
            cth = 1.0 - 2.0 * r1
            sth = math.sqrt(max(0.0, 1.0 - cth * cth))
            ph = 2.0 * math.pi * r2
            dx, dy, dz = sth * math.cos(ph), sth * math.sin(ph), cth
            _forzar(acum, t, mu, nz, ny, nx, vox, x, y, z, 0.0, 0.0, 0.0, False, E, w, angulos, radio_cm, cz, pix_mm, npx, npy,
                    E_tab, e_min, de, ne, r_agua, r_hueso, False, paso)
            estad[t, 0] += 1
            vivo = True
            n_comp = 0
            while vivo:
                ie = _indice_e(E, e_min, de, ne)
                mu_max = mu_max_tab[ie]
                s, r = _u01(s)
                sdist = -math.log(r + 1e-300) / mu_max
                x += dx * sdist
                y += dy * sdist
                z += dz * sdist
                i, j, k = int(z / vox), int(y / vox), int(x / vox)
                if i < 0 or i >= nz or j < 0 or j >= ny or k < 0 or k >= nx or x < 0.0 or y < 0.0 or z < 0.0:
                    vivo = False      # escapó
                    break
                hueso = es_hueso[i, j, k]
                razon = r_hueso[ie] if hueso else r_agua[ie]
                mu_real = mu[i, j, k] * razon
                s, r = _u01(s)
                if r > mu_real / mu_max:
                    continue          # interacción ficticia
                fc = f_hueso[ie] if hueso else f_agua[ie]
                s, r = _u01(s)
                if r > fc:
                    estad[t, 1] += 1
                    vivo = False      # fotoeléctrico
                    break
                # Compton
                s, c, E2 = _kahn(s, E)
                n_comp += 1
                estad[t, 2] += 1
                # nueva dirección: rotar (dx,dy,dz) por theta=acos(c) y azimut uniforme
                s, r = _u01(s)
                az = 2.0 * math.pi * r
                sn = math.sqrt(max(0.0, 1.0 - c * c))
                if abs(dz) > 0.999999:
                    ndx, ndy, ndz = sn * math.cos(az), sn * math.sin(az), c * (1.0 if dz > 0 else -1.0)
                else:
                    den = math.sqrt(1.0 - dz * dz)
                    ndx = dx * c + sn * (dx * dz * math.cos(az) - dy * math.sin(az)) / den
                    ndy = dy * c + sn * (dy * dz * math.cos(az) + dx * math.sin(az)) / den
                    ndz = dz * c - sn * den * math.cos(az)
                # detección forzada con la dirección ANTES de dispersar (el ángulo hacia cada detector se evalúa adentro)
                _forzar(acum, t, mu, nz, ny, nx, vox, x, y, z, dx, dy, dz, True, E, w, angulos, radio_cm, cz, pix_mm, npx, npy,
                        E_tab, e_min, de, ne, r_agua, r_hueso, False, paso)
                dx, dy, dz, E = ndx, ndy, ndz, E2
                if E < 90.0 or n_comp > 10:
                    vivo = False
                    break


class Camara:
    def __init__(self, n_angulos=60, arco_grados=360.0, inicio_grados=0.0, radio_cm=20.0, matriz=128, pix_mm=3.3, tiempo_s=25.0):
        self.angulos = np.deg2rad(inicio_grados + np.arange(n_angulos) * arco_grados / n_angulos).astype(np.float64)
        self.radio_cm, self.matriz, self.pix_mm, self.tiempo_s = radio_cm, matriz, pix_mm, tiempo_s


def simular(actividad_kbq_ml: np.ndarray, mu: np.ndarray, hu: np.ndarray, vox_cm: float, camara: Camara, n_hist: int = 2_000_000,
            semilla: int = 1, paso_cm: float = 0.1, n_threads: int | None = None, poisson: bool = True):
    """Devuelve (proyecciones esperadas float64 [ang, v, u], proyecciones con ruido int32, estadísticas)."""
    import numba
    if n_threads is None:
        n_threads = numba.get_num_threads()
    E_tab, r_agua, f_agua, r_hueso, f_hueso = tablas()
    e_min, de, ne = float(E_tab[0]), float(E_tab[1] - E_tab[0]), len(E_tab)
    es_hueso = (hu > HU_HUESO)
    mu_soft_max = float(mu[~es_hueso].max()) if (~es_hueso).any() else 0.0
    mu_bone_max = float(mu[es_hueso].max()) if es_hueso.any() else 0.0
    mu_max_tab = np.maximum(mu_soft_max * r_agua, mu_bone_max * r_hueso) + 1e-6
    act = actividad_kbq_ml.astype(np.float64).ravel()
    cdf = np.cumsum(act)
    total_kbq = cdf[-1] * (vox_cm ** 3)        # kBq (kBq/mL · mL)
    cdf /= cdf[-1]
    nz, ny, nx = mu.shape
    acum = np.zeros((n_threads, len(camara.angulos), camara.matriz, camara.matriz), np.float64)
    estad = np.zeros((n_threads, 3), np.int64)
    cz = nz * vox_cm * 0.5
    _simular(n_hist, semilla, cdf, mu.astype(np.float64), es_hueso, float(vox_cm), camara.angulos, float(camara.radio_cm), cz,
             float(camara.pix_mm), camara.matriz, camara.matriz, acum, E_tab, e_min, de, ne, r_agua, r_hueso, f_agua, f_hueso,
             mu_max_tab, paso_cm, n_threads, estad)
    prob = acum.sum(axis=0)                    # probabilidad de detección por píxel, por fotón emitido
    emitidos = total_kbq * 1000.0 * RENDIMIENTO_GAMMA * camara.tiempo_s
    esperado = prob * (emitidos / n_hist)
    rng = np.random.default_rng(semilla)
    ruido = rng.poisson(esperado).astype(np.int32) if poisson else np.round(esperado).astype(np.int32)
    e = estad.sum(axis=0)
    return esperado, ruido, {"historias": int(e[0]), "fotoelectricos": int(e[1]), "compton": int(e[2]), "emitidos_por_proyeccion": float(emitidos),
                             "actividad_MBq": total_kbq / 1000.0, "cuentas_medias_por_proyeccion": float(esperado.sum(axis=(1, 2)).mean())}
