#!/usr/bin/env bash
# Verifica que great_lake.py da los mismos resultados que el Fortran original.
# Uso: tools/verify.sh [directorio_con_data.dat]   (por defecto ejemplos/caso_base)
# Requiere gfortran, python3 y numpy.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CASE="$(cd "${1:-$ROOT/ejemplos/caso_base}" && pwd)"
WORK="$(mktemp -d)"
mkdir -p "$WORK/fortran" "$WORK/python"
python3 "$ROOT/tools/restore_fortran.py" "$ROOT/Great_lake_v.1.f90" "$WORK/gl.f90"
gfortran -std=legacy -O0 -ffree-line-length-none -fcheck=all -o "$WORK/gl" "$WORK/gl.f90"
cp "$CASE/data.dat" "$WORK/fortran/"; cp "$CASE/data.dat" "$WORK/python/"
(cd "$WORK/fortran" && "$WORK/gl" > /dev/null 2>&1)
# gfortran conecta la unidad 6 a la salida estandar: al abrir la unidad 6 como
# Matriz_roturas.out, los mensajes de pantalla acaban en ese fichero. Se filtran.
grep -vE "Great_lake_v.1|Time steep|tiempo CPU" "$WORK/fortran/Matriz_roturas.out" > "$WORK/tmp" || true
mv "$WORK/tmp" "$WORK/fortran/Matriz_roturas.out"
python3 "$ROOT/great_lake.py" "$WORK/python" > /dev/null
python3 "$ROOT/tools/compare_outputs.py" "$WORK/fortran" "$WORK/python"
echo "Resultados en: $WORK"
