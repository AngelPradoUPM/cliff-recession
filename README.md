# cliff-fall

Modelo de erosión y retroceso de acantilados de arcillas (tills) en lagos sin
marea (lago Erie), de R. Castedo (UPM). La tesis doctoral relacionada está en
`tesis/`.

## Contenido

| Fichero | Descripción |
|---|---|
| `Great_lake_v.1.f90` | Código original en Fortran (tal como se entregó). |
| `great_lake.py` | Traducción a Python del código Fortran. |
| `ejemplos/caso_base/data.dat` | Fichero de datos de entrada del caso base. |
| `tools/restore_fortran.py` | Convierte el `.f90` original en Fortran compilable con gfortran (referencia). |
| `tools/compare_outputs.py` | Compara, valor a valor, las salidas de ambas versiones. |
| `tools/verify.sh` | Compila el Fortran, ejecuta las dos versiones y compara los resultados. |

## Uso de la versión en Python

Requiere Python 3 y NumPy.

```bash
python great_lake.py ejemplos/caso_base
```

Lee `data.dat` del directorio indicado y escribe allí `timerecesion.out`,
`Matriz_dibujo.out`, `Matriz_roturas.out`, `Matriz_recesion.out` y
`Matriz_seepage.out`, con el mismo formato que el programa Fortran.

## Verificación

```bash
tools/verify.sh                 # caso base
tools/verify.sh otro/directorio # cualquier otro data.dat
```

Con el caso base (1260 meses) las cinco salidas son **idénticas** a las del
Fortran (más de 300 000 valores, diferencia máxima 0). El perfil completo
también coincide bit a bit tras cada paso de tiempo. También se comprobaron dos
variantes con playa (`Qb > 0`) y nivel freático (`groundwater = 1` y `0`), con
resultados idénticos.

La traducción es literal a propósito para reproducir exactamente el
comportamiento numérico del original:

- Los índices empiezan en 1 y se respetan los bucles DO y los GO TO.
- Las divisiones enteras de Fortran se mantienen.
- Las constantes sin sufijo `d0` (`9.81`, `1.56`, `3.1416`, `0.05`…) se tratan
  como de precisión simple, incluido `pi`.
- Las funciones elementales son las de la libm de C.

## Aspectos del código original a revisar

Se han conservado tal cual en la traducción para no cambiar los resultados,
pero conviene revisarlos:

1. **Esfuerzo cortante nulo**: en `SS = (1/2)*Fw*seawater*(Uo**2)`, `(1/2)` es
   una división entera y vale 0. La erosión por cortante es
   `No*Kss*(-SFmin)**p`. Con `p = 1` sale negativa y se anula, así que nunca
   actúa. Con un `p` no entero el resultado sería NaN. Seguramente debería ser
   `0.5`.
2. **`pi` en precisión simple**: `real*8, parameter :: pi = 3.141592653589793`
   sin `d0` vale 3.1415927410…
3. **`WBD = (int(WH(j)*1.28*10))/10`**: la profundidad de rotura se trunca a
   metros enteros (división entera), no a decímetros.
4. **`Volacumulafor`**: el fichero entregado usa esta variable, que no está
   declarada y no compilaría con `IMPLICIT NONE`. Es un residuo del reemplazo
   automático `do` → `for`; se ha interpretado como `Volacumulado`, que
   controla el tiempo de protección por derrubios.
5. **Variación interanual del nivel del lago desactivada**: `changelakelevel`
   está desactivada y las líneas `TL = controllevel…` están comentadas, así que
   el nivel del lago no cambia.
6. **`lakelevel` es INTEGER** pero `data.dat` trae `0.0`. Con gfortran esa
   lectura falla. La versión Python lee un real y lo trunca.
7. **Con gfortran**, `open(6, file='Matriz_roturas.out')` hace que los mensajes
   de pantalla se escriban también en ese fichero, porque la unidad 6 es la
   salida estándar. Con Intel Fortran, el compilador original, no ocurre.

`marea.dat` no se usa: el `open` correspondiente está comentado en el código.
