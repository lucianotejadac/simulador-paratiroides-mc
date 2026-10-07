# Simulador Monte Carlo de cintigrafía de paratiroides

Genera una adquisición SPECT de paratiroides con Tc-99m sestamibi, de cuello y mediastino hasta el corazón, con panorámicas anteriores, con la física real (atenuación,
dispersión Compton, colimador, resolución y ruido de Poisson), a partir de los cortes axiales de un CT.
El CT da el mapa de atenuación; la biodistribución del radiofármaco se construye sobre su segmentación
y el adenoma se coloca a voluntad, así que cada estudio simulado tiene la verdad conocida.

Página con los cinco casos (cortes fusionados, cine de proyecciones, verdad oculta y DICOM descargables):
<https://lucianotejadac.github.io/simulador-paratiroides-mc/>

Pensado para docencia: proyecciones y volúmenes reconstruidos en DICOM listos para el
[visor_dicom](https://github.com/lucianotejadac/visor_dicom) y para MicroDicom, fusionables con el CT.

## Flujo

```
CT (TCIA, público) ──► fantoma.py ──► mu(140 keV) a 2 mm ──► segmentar.py ──► regiones
                                                                      │
                                   actividad.py ◄─────────────────────┘  (casos: normal + 4 adenomas)
                                        │
                               simular_caso.py  (Monte Carlo, 60 proyecciones × 2 fases)
                                        │
                               reconstruir.py (OSEM con corrección de atenuación)
                                        │
                               exportar_dicom.py (CT + NM RECON TOMO + proyecciones)
```

```bash
python -m venv .venv && .venv/Scripts/pip install numpy pydicom numba==0.61.2 pillow scipy pytest
python src/fantoma.py            # lee el CT, recorta el cuello, remuestrea, HU -> mu
python src/segmentar.py          # regiones por umbral (hueso, vía aérea, salivales)
python src/totalseg.py           # TotalSegmentator en CPU: cuello y tórax
TotalSegmentator -i salida/totalseg/cuello_ct.nii.gz -o salida/totalseg/total.nii.gz --ml -d cpu --roi_subset thyroid_gland esophagus trachea ...
python src/comparar_totalseg.py  # compara y guarda las etiquetas en la grilla del fantoma
python src/regiones_totalseg.py  # tiroides, tráquea, esófago y vasos de TotalSegmentator
python src/actividad.py          # mapas de actividad por caso y fase
python src/simular_caso.py --caso inferior-derecho-12 --historias 2000000
python src/reconstruir.py --caso inferior-derecho-12 --fase precoz
python src/ct_alta.py              # CT original (0.98 mm) para mostrar y exportar
python src/exportar_dicom.py --caso inferior-derecho-12 --numero 1
python -m pytest -q tests        # validación del motor contra lo conocido de una LEHR
```

## Contornos

Tiroides, tráquea, esófago, vasos del cuello, vértebras y clavículas vienen de
[TotalSegmentator](https://github.com/wasserth/TotalSegmentator) (Wasserthal et al., 2023), corrido en CPU.
Los adenomas se colocan donde la esfera entera cabe en tejido blando libre, lo más cerca del sitio
anatómico pedido. Requisitos extra: `torch==2.8.0` y `torchvision==0.23.0` de CPU, `totalsegmentator`,
`pandas==2.2.3`, `connected-components-3d==3.18.0` (versiones que pasan el Control de aplicaciones de Windows).

## El motor (`src/montecarlo.py`)

- Cada historia es un fotón de 140 keV emitido desde el mapa de actividad. Transporte por el método de
  Woodcock en la grilla de μ; interacciones fotoeléctrica y Compton (Klein-Nishina por el método de
  Kahn) con tablas de NIST para agua y hueso cortical; Rayleigh omitido.
- Detección forzada hacia todos los ángulos desde la emisión y desde cada vértice Compton, con la
  energía que el fotón tendría al salir hacia cada detector, la transmisión del rayo, la eficiencia
  del cristal, la ventana de energía con resolución gaussiana y la PSF del colimador a esa distancia.
- Colimador LEHR (Siemens): eficiencia geométrica 1.17 × 10⁻⁴, resolución 7 mm a 10 cm.
- Numba en paralelo, generador de azar propio por hilo con semilla: resultados reproducibles bit a bit.
- 2 millones de historias con 60 ángulos: ~1 minuto en una laptop de 16 hilos.

## Resultados de esta ronda

Cinco casos sobre el mismo cuello (CT caso 2 de la entrega PET/CT de TCIA): normal, adenoma inferior
derecho de 12 mm, inferior izquierdo de 8 mm, retroesofágico de 10 mm y mediastínico de 15 mm, cada
uno con fase precoz (15 min) y tardía (2 h). Detalles, trampas y validación en [BITACORA.md](BITACORA.md).

## Datos y licencias

- CT base: The Cancer Imaging Archive, colección de la entrega docente PET/CT (desidentificado según
  DICOM PS 3.15 Anexo E). Los estudios simulados conservan la anatomía de ese paciente anónimo.
- Código: MIT. Los DICOM generados llevan paciente sintético `SIM-PARA-nn`.

## Estructura

```
src/fantoma.py        CT -> fantoma de atenuación (npz + json + montaje png)
src/segmentar.py      regiones por umbral (npz) y montajes de control
src/comparar_totalseg.py  TotalSegmentator vs umbral; posición de los adenomas
src/regiones_totalseg.py  regiones finales con las estructuras de TotalSegmentator
src/actividad.py      casos y mapas de actividad (npy + verdad.json)
src/montecarlo.py     motor Monte Carlo (Numba)
src/simular_caso.py   adquisición de un caso (proyecciones npz + png + adquisicion.json)
src/reconstruir.py    OSEM con atenuación y filtro posterior
src/exportar_dicom.py DICOM CT + NM (reglas del visor_dicom)
src/verificar.py      comprueba los DICOM contra las reglas del visor y de la fusión
src/exportar_web.py   datos compactos y ZIP para la página (docs/)
src/ct_alta.py        CT a la resolución del tomógrafo para la página y los DICOM
src/contornos_web.py  contornos de las regiones para la página
docs/                 página de GitHub Pages con los cinco casos
tests/                validación del motor
salida/               resultados (los pesados no se versionan)
```
