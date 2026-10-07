# Bitácora de decisiones

Una entrada por ronda de trabajo, al estilo ADR: qué se decidió, por qué, qué se vio y qué queda
abierto. Las trampas (resultados que parecían bien y no lo estaban) se anotan antes de corregirlas.

---

## 0001 · 2026-10-07 · De un CT público a una cintigrafía de paratiroides simulada

**Pregunta.** ¿Se puede producir una adquisición SPECT de paratiroides con sestamibi, con la física
real (atenuación, dispersión, colimador, ruido), a partir de los cortes axiales de un CT, de modo que
el estudio sirva en docencia con la verdad conocida?

**Decisiones de partida** (acordadas antes de escribir código):
1. CT base: caso 2 de los PET/CT de TCIA ya descargados (hombre de 62 años, GE Discovery STE, 140 kVp,
   sin contraste, paso 3.3 mm). Elegido entre cinco porque es el único con cuello recto, sin artefactos
   dentales y con 15 cm de cabeza por encima; el caso 4 también cubre el cuello pero con cifosis,
   metal y la cabeza inclinada.
2. Fantoma: recorte desde la mandíbula hasta 12 cm bajo el mínimo del cuello, remuestreo a 2 mm
   isotrópicos (96 × 132 × 240), conversión bilineal HU → μ a 140 keV (agua 0.1537 /cm, hueso cortical
   0.286 /cm en 1000 HU).
3. Segmentación por densidad y referencias anatómicas (`src/segmentar.py`): tiroides por su
   hiperdensidad en CT sin contraste (70–200 HU) a ambos lados de la tráquea; tráquea por el aire
   central; hueso y vía aérea por umbral; glándulas salivales como elipsoides referidos a la tráquea
   (solo importan como captación de fondo).
4. Biodistribución (kBq/mL, precoz/tardía): tejido blando 4/2.5, hueso 3/2, tiroides 40/16, salivales
   30/18, adenoma = tiroides precoz × relación (2 a 3), con retención 0.8 en la tardía. 21.5 MBq en
   el campo en la fase precoz.
5. Casos: normal, inferior derecho 12 mm (relación 2.5), inferior izquierdo 8 mm (2.0), retroesofágico
   10 mm (2.5) y mediastínico 15 mm (3.0), con la verdad en `verdad.json`.
6. Cámara: LEHR de Siemens (agujero 1.11 mm, largo 24.05, septo 0.16), cristal de 9.5 mm, resolución
   energética 9.5 % a 140 keV, ventana 126–154 keV, 60 proyecciones en 360°, radio 20 cm, matriz 128
   de 3.3 mm, 25 s por proyección.
7. Motor: Monte Carlo en CPU con Numba (`src/montecarlo.py`), paralelo por hilos, generador propio por
   hilo con semilla explícita (reproducible bit a bit, verificado en test). Woodcock para el transporte,
   fotoeléctrico y Compton (Kahn), Rayleigh omitido, detección forzada por convolución hacia todos los
   ángulos desde cada vértice. La versión WebGPU queda para después de validar esta.

**Trampas encontradas y corregidas.**
- **Relleno de huecos en 3D.** La máscara corporal con `binary_fill_holes` en 3D dejaba fuera la
  tráquea y los vértices pulmonares, porque en 3D están conectados con el exterior por la boca y por
  el borde inferior del recorte. El relleno corte a corte (2D) los recupera: la tráquea pasó de 2 a
  53 mL.
- **El cartílago tiroides también mide 70–200 HU.** La primera segmentación de la tiroides incluía
  las láminas calcificadas del cartílago, 3 cm por encima de la glándula. Acotada a 3–5.4 cm bajo el
  mínimo del cuello: 11.8 mL, dentro del rango normal, en los cortes correctos (revisado a ojo en
  `salida/tiroides_zoom.png`).
- **Energía equivocada en la detección forzada.** La primera versión usaba la energía del fotón
  *antes* de dispersar para calcular la ventana y la transmisión hacia cada detector; así toda la
  dispersión entraba en ventana y la fracción de dispersión salía del 120 % (literatura: 25–35 %
  para Tc-99m a 10 cm de profundidad con ventana del 20 %). Corregido: para cada ángulo se calcula la
  energía que tendría el fotón al salir hacia *ese* detector. Test `test_atenuacion_y_dispersion_en_agua`.
- **Klein-Nishina sin normalizar.** Faltaba el factor 1/2 de la sección diferencial; el test de
  integral igual a 1 sobre la esfera lo fija.
- **Fuente puntual de un centímetro.** El test de resolución medía 10 mm en vez de 7 porque la
  "fuente puntual" ocupaba un vóxel de 1 cm; con vóxel de 1 mm da 7.0 mm a 10 cm, igual a la fórmula
  del colimador más la intrínseca.
- **El reconstructor giraba al revés.** `ndimage.rotate` con `axes=(2,1)` gira en sentido contrario al
  que yo había supuesto: las proyecciones del Monte Carlo y las del proyector de OSEM coincidían en
  0° y 180° y estaban espejadas en 90° y 270°. La reconstrucción mostraba el contorno del cuello
  (simétrico) y ningún lóbulo (asimétricos, desdibujados en un anillo). Detectado con una fuente
  puntual descentrada proyectada por los dos caminos; corregido el signo.

**Validación del motor** (`tests/test_montecarlo.py`, 5 pruebas):
- Sensibilidad de fuente puntual en aire: igual a g·ε·P_ventana dentro del 2 % (1.02 × 10⁻⁴ cuentas
  por fotón emitido; LEHR nominal ~1.1 × 10⁻⁴).
- Resolución a 10 cm: 7.0 mm (colimador 5.9 ⊕ intrínseca 3.8).
- Fuente a 10 cm de profundidad en agua: primarios atenuados exp(−μ·10) y total entre 1.15 y 1.6 veces
  los primarios; Compton domina sobre fotoeléctrico por más de 5 a 1.
- Reproducibilidad exacta con la misma semilla.

**Rendimiento.** 2 millones de historias con 60 ángulos de detección forzada: ~60 s en los 16 hilos
del Ryzen 9 8945HS; cada caso con dos fases, unos 2 minutos. Proyecciones con 22 000 cuentas de
media (precoz), del orden de lo clínico para la actividad en el campo.

**Qué queda abierto.**
- Validar contra SIMIND (sensibilidad, fracción de dispersión, perfiles) cuando esté instalado.
- La transmisión del rayo hacia el detector escala μ con la tabla del agua también en hueso (error
  pequeño a 140 keV; a 100 keV el hueso atenúa 10 % más de lo calculado).
- Penetración septal y Rayleigh, irrelevantes para Tc-99m con LEHR, necesarios para I-131.
- Glándulas salivales como elipsoides: reemplazables por segmentación manual.
- Órbita de contorno (la circular de 20 cm deja el cuello lejos del colimador; la clínica acerca).
- Versión WebGPU para el navegador.

---

## 0002 · 2026-10-07 · Contorneo con TotalSegmentator: los adenomas «inferiores» eran tiroideos

**Pregunta del usuario.** ¿Cómo se validó la ubicación de los órganos? Respuesta honesta: la
tiroides se revisó a ojo, el resto eran umbrales y elipsoides sin validación anatómica. Medida la
posición de cada adenoma contra el CT, el retroesofágico tenía el 19 % de su volumen dentro de C7.

**Decisión.** Contornear con una herramienta independiente: TotalSegmentator 2 (Wasserthal et al.,
Radiology: AI 2023), red nnU-Net entrenada sobre CT, corrida localmente en CPU sobre el fantoma
exportado a NIfTI con la geometría exacta (la salida cae en la misma grilla). Estructuras: tiroides,
tráquea, esófago, C3–T2, clavículas, lóbulos superiores, carótidas, subclavias, tronco y venas
braquiocefálicas. 85 s de inferencia.

**Lo que mostró.**
- Mi tiroides por umbral estaba corrida 16 mm hacia craneal: cubría los cortes 15–32 y la de
  TotalSegmentator los 7–30; los cortes superiores de la mía eran cartílago laríngeo. Volúmenes
  parecidos (11.8 y 13.2 mL) pero Dice 0.54.
- Por eso los adenomas «inferiores», colocados bajo mi polo inferior, caían dentro de la tiroides
  real (98 de 107 vóxeles en el derecho): simulaban nódulos tiroideos, no paratiroides.
- **Lección:** un volumen plausible no valida una segmentación; hay que medir solape contra una
  referencia independiente. Y una segmentación por umbral de densidad confunde tejidos con la
  misma densidad (cartílago calcificado y tiroides con yodo).

**Cambios.**
- `src/regiones_totalseg.py`: tiroides, tráquea, esófago, vasos (pool sanguíneo, 6/3 kBq/mL),
  vértebras y clavículas de TotalSegmentator; se conserva del umbral solo lo que la red no cubre
  (laringe y vía aérea alta, salivales, resto del hueso). La versión anterior queda en
  `salida/regiones_umbral.npz`.
- `src/actividad.py`: el adenoma se coloca en el vóxel más cercano al objetivo anatómico donde la
  esfera entera cabe en tejido blando libre (−200 a 150 HU, sin tiroides, esófago, tráquea, vasos,
  hueso ni pulmón), calculado como erosión del tejido libre por la esfera. La verdad registra el
  objetivo y el desvío.
- Resultado: inferior derecho en el surco traqueoesofágico bajo el polo (a 3 mm del objetivo),
  inferior izquierdo lateral a la tráquea junto al polo (4 mm), retroesofágico posterolateral
  izquierdo al esófago, prevertebral (14 mm: detrás del esófago no hay 10 mm libres en este
  cuello), mediastínico superior anterior a nivel de la escotadura esternal (13 mm). Ninguno toca
  estructuras vecinas (bordes a 0.3–1.6 mm).

**Entorno (trampas del Control de aplicaciones de Windows, otra vez).** Bloqueados: pandas 3.0.6
(pasa 2.2.3), connected-components-3d 4.1.0 (pasa 3.18.0), la extensión compilada de torchvision
0.29 (pasan torch 2.8.0 + torchvision 0.23.0 de CPU). torch 2.9 también cargaba, pero nnU-Net lo
excluye explícitamente; se usó 2.8.

**Queda abierto.** Glándulas salivales siguen como elipsoides (la tarea de cabeza de
TotalSegmentator pide licencia académica); revisión de los contornos por un profesional en 3D
Slicer; un adenoma en el polo superior (sitio más frecuente de las paratiroides superiores) que hoy
no está entre los casos.

---

## 0003 · 2026-10-07 · El CT a la resolución del tomógrafo

**Pedido.** El CT se veía borroso: en la página se mostraba remuestreado a la grilla del SPECT
(3.3 mm) y en los DICOM a la del fantoma (2 mm).

**Decisión.** El fantoma de 2 mm sigue siendo el mapa de atenuación del Monte Carlo (la física no
cambia); para mostrar y exportar se usa el CT original, con el mismo recorte del cuello, la misma
máscara corporal y las mismas coordenadas LPS (`src/ct_alta.py`).
- DICOM: 59 cortes de 0.98 mm en el plano y 3.27 mm entre cortes (los del tomógrafo, espaciado
  uniforme verificado), misma identidad y `FrameOfReferenceUID` que el SPECT: la fusión sigue
  alineada. Los cinco casos pasan `verificar.py`.
- Página: el CT se muestrea a 4 sub-píxeles por píxel SPECT (0.825 mm, 512 × 512 por corte) en los
  mismos cortes de la grilla, se guarda una sola vez (`docs/datos/ct_hd.bin`, uint16 HU+1024,
  recortado al cuerpo, 17.5 MB) y el SPECT se superpone como capa suavizada. Se agregaron ventanas
  de partes blandas y hueso, «solo CT» y zoom (1.6× por defecto); el umbral por defecto del SPECT
  sube para que el fondo corporal no tiña todo el corte.

**Agregado (mismo día): contornos en la página.** `src/contornos_web.py` lleva las regiones finales
(TotalSegmentator + umbral) a los cortes de la página a 1.65 mm (vecino más cercano desde la grilla
de 2 mm), 2.2 MB. El visor dibuja el borde de cada estructura sobre el CT con un interruptor por
estructura: tiroides, esófago, tráquea y vasos encendidos; salivales (marcadas «aprox.», son
elipsoides), hueso y pulmón apagados. Revisados a ojo en un corte medio tiroideo: lóbulos alrededor de
la tráquea, esófago posterior y carótidas laterales.
