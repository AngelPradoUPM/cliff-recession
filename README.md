# cliff-recession

Modelo de erosión y retroceso de acantilados de arcillas (tills) en lagos sin
marea (lago Erie), de R. Castedo (UPM). La tesis doctoral relacionada está en
`tesis/`.

## Contenido

| Fichero | Descripción |
|---|---|
| `Great_lake_v.1.f90` | Código original en Fortran (Intel Visual Fortran). |
| `great_lake.py` | Traducción a Python del código Fortran (motor de cálculo). |
| `great_lake_gui.py` | Herramienta gráfica de escritorio sobre `great_lake.py`. |
| `ejemplos/caso_base/data.dat` | Caso base: 1260 meses, sin playa, terreno saturado, un estrato. |
| `ejemplos/caso_simple/data.dat` | Caso simple y realista: acantilado de 10 m a 45° (estable al inicio), 20 años, playa pequeña, nivel freático, dos estratos y un punto de estudio. Parámetros dentro de los rangos de la tesis (cap. 10). |
| `ejemplos/caso_extremo/data.dat` | Caso de comprobación con parámetros desfavorables (fuera de los rangos de la tesis). Retrocede muy rápido, pero ejercita todas las partes del código. |
| `tools/restore_fortran.py` | Adapta el `.f90` para compilarlo con gfortran (cambios mínimos de sintaxis). |
| `tools/compare_outputs.py` | Compara valor a valor las salidas de ambas versiones. |
| `tools/verify.sh` | Compila el Fortran, ejecuta las dos versiones y compara los resultados. |

## Uso de la versión en Python

Requiere Python 3 y NumPy.

```bash
python great_lake.py ejemplos/caso_simple
```

Lee `data.dat` del directorio indicado y escribe allí `timerecesion.out`,
`Matriz_dibujo.out`, `Matriz_roturas.out`, `Matriz_recesion.out` y
`Matriz_seepage.out`, con el mismo formato que el programa Fortran.

### Nivel del lago variable (opcional)

Por defecto el nivel del lago es constante, como en el Fortran original. El
original tiene comentada una variación del nivel, y se puede activar:

```bash
python great_lake.py ejemplos/caso_simple --nivel-estacional --nivel-interanual
```

- `--nivel-estacional`: aplica los niveles de la tabla 10.iii de la tesis,
  +0,40 m en primavera y +0,30 m en verano respecto a invierno y otoño.
- `--nivel-interanual`: añade un ciclo de 15 años de ±5·Δz (±0,5 m con
  Δz = 0,1 m), con cambios cada 9 meses desde el mes 69.

Con estas opciones activadas, los resultados ya no coinciden con el Fortran.

El ciclo interanual comentado en el original perdía los cambios que caían en
invierno o durante un salto por derrubios, y el nivel derivaba hacia arriba
(+1,4 m en 100 años). Aquí se aplican todos los cambios pendientes y el nivel
oscila correctamente.

## Herramienta gráfica

```bash
pip install -r requirements.txt
python great_lake_gui.py ejemplos/caso_simple/data.dat
```

![Herramienta gráfica con el caso base](docs/captura_gui.png)

- **Panel «Datos»** (izquierda):
  - edición del perfil inicial (cota, distancia y tipo de punto);
  - edición de los parámetros y de los estratos;
  - lista de los perfiles de salida, con casillas «Ver».
- **Vista central**: el perfil con el nivel del lago, el nivel freático y los
  puntos singulares. Se actualiza al editar los datos. Tras calcular, muestra
  los perfiles en cada instante de salida (escala de color por año). Tiene dos
  vistas, «Acantilado» y «Perfil completo», y zoom con la barra de matplotlib.
- **Pestañas de la derecha**:
  - *Recesión*: gráficas de retroceso y de volumen, que se pueden exportar a
    PNG o CSV.
  - *Roturas*: el perfil antes y después de cada rotura.
  - *Oleaje*: tabla editable del oleaje por estación. Por defecto tiene los
    valores del código original; si se modifica, los resultados dejan de
    coincidir con el Fortran.
  - *Nivel / mareas*: activa la variación del nivel del lago en el cálculo
    (estacional, con los desniveles editables, y ciclo interanual), y muestra
    el nivel usado en cada mes. También carga y representa una serie de
    niveles externa (por ejemplo, `marea.dat`), que solo se visualiza.
- **Línea temporal** (debajo de la vista central): una barra que se arrastra
  para ver el perfil mes a mes. También se maneja con las flechas (un mes),
  RePág/AvPág (un año) y Inicio/Fin. Tiene botones de reproducir y avanzar o
  retroceder un mes, y un selector de velocidad. En la barra se marcan los
  meses con roturas (rojo), los meses de invierno, que no se calculan, y los
  meses de protección por derrubios. En los meses con rotura se dibuja, con
  línea discontinua, el perfil antes de cada rotura. Si el nivel es variable,
  también se dibuja el nivel del lago de ese mes.
- **Calcular**: se ejecuta en segundo plano, con barra de progreso y botón
  para detener. Los ficheros de salida se guardan en la carpeta `resultados/`,
  junto al `data.dat`. Los resultados son idénticos a los de `great_lake.py`
  por línea de órdenes.

## Verificación

```bash
tools/verify.sh ejemplos/caso_base
tools/verify.sh ejemplos/caso_simple
tools/verify.sh ejemplos/caso_extremo
```

En los tres casos las cinco salidas son **idénticas** a las del Fortran
compilado con gfortran (diferencia máxima 0), con el nivel del lago constante
del original. Entre los tres casos se ejercitan todos los bloques del modelo:

- erosión por esfuerzo cortante, por impacto del oleaje y por run-up;
- formación de playa y berma;
- roturas circulares, incluidas las descartadas por "crear material", y
  roturas por vuelco;
- salto temporal por protección de derrubios;
- nivel freático de los tipos 1 y 2.

La traducción es literal a propósito para reproducir exactamente el
comportamiento numérico del original:

- Los índices empiezan en 1 y se respetan los bucles DO y los GO TO.
- Se respeta la división entera.
- Las constantes escritas sin `d0` en el original son de precisión simple,
  incluido `pi = 3.141592653589793`, que vale 3.1415927410…
- Las funciones elementales son las de la libm de C.

La versión en Python es más lenta que la de Fortran: unos 60 s frente a ~1 s
en el caso base.

## Notas sobre el código original

- `lakelevel` es INTEGER, pero `data.dat` trae `0.0`. gfortran no acepta esa
  lectura, así que la referencia lee un real y lo trunca, igual que la versión
  Python. Intel Fortran sí la acepta.
- Varias constantes no llevan `d0` (`pi`, `1.56`, `9.81` en `freshwater`, las
  frecuencias del oleaje…), por lo que se evalúan en precisión simple. El
  efecto es pequeño, pero conviene saberlo.
- La variación del nivel del lago está desactivada en el original:
  `changelakelevel` vale 0 y las líneas `TL = controllevel…` están comentadas.
  En la versión Python se puede activar como opción (ver más arriba).
- La tesis da un coeficiente de amortiguación de la fuerza de las olas en la
  zona de rompientes χ = 0,01, mientras que el código usa 0,05
  (`dexp(-SZ*0.05)`).
- Tasas de retroceso de referencia (tesis, cap. 10): en el lago Erie, entre
  0,1 m y unos pocos metros al año. En los casos con playa pequeña y nivel
  variable, 0,25–0,45 m/año. El caso base da ~1,2 m/año en el nivel del lago y
  `caso_simple` ~0,6 m/año en la cresta.
- Con gfortran, `open(6, file='Matriz_roturas.out')` hace que los mensajes de
  pantalla se escriban también en ese fichero, porque la unidad 6 es la salida
  estándar. `tools/verify.sh` los filtra. Con Intel Fortran no ocurre.
- `marea.dat` no se usa: el `open` correspondiente está comentado.
