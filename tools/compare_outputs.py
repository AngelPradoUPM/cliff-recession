#!/usr/bin/env python3
"""Compara las salidas del programa Fortran y de la traduccion a Python.

Uso: python compare_outputs.py DIR_FORTRAN DIR_PYTHON

Para cada fichero de salida compara linea a linea y valor a valor. Ignora el
tiempo de CPU (fila 1, columna 2 de Matriz_dibujo.out).
"""
import math
import sys
from pathlib import Path

FILES = ["timerecesion.out", "Matriz_dibujo.out", "Matriz_roturas.out",
         "Matriz_recesion.out", "Matriz_seepage.out"]


def parse(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        toks = line.split()
        rows.append([tok for tok in toks])
    return rows


def to_float(tok):
    try:
        return float(tok)
    except ValueError:
        if tok.lower() in ("nan",):
            return math.nan
        if tok.lower().lstrip("+") in ("infinity", "inf"):
            return math.inf if not tok.startswith("-") else -math.inf
        return None  # p.ej. '*********' (desbordamiento del formato)


def main(dir_f, dir_p):
    ok = True
    for name in FILES:
        rf, rp = parse(Path(dir_f) / name), parse(Path(dir_p) / name)
        if len(rf) != len(rp):
            print(f"{name}: numero de lineas distinto ({len(rf)} vs {len(rp)})")
            ok = False
            continue
        n = ndiff_text = 0
        maxdiff = 0.0
        where = None
        for i, (a, b) in enumerate(zip(rf, rp)):
            if len(a) != len(b):
                print(f"{name}: linea {i+1} con distinto numero de valores")
                ok = False
                break
            for k, (x, y) in enumerate(zip(a, b)):
                if name == "Matriz_dibujo.out" and i == 0 and k == 1:
                    continue  # tiempo de CPU
                n += 1
                if x != y:
                    ndiff_text += 1
                fx, fy = to_float(x), to_float(y)
                if fx is None or fy is None:
                    if x != y:
                        maxdiff = math.inf
                    continue
                if math.isnan(fx) and math.isnan(fy):
                    continue
                d = abs(fx - fy)
                if d > maxdiff:
                    maxdiff, where = d, (i + 1, k + 1, x, y)
        status = "IDENTICO" if ndiff_text == 0 else "DIFERENTE"
        if ndiff_text:
            ok = False
        print(f"{name:22s} {n:8d} valores  {status:9s}  textos distintos: {ndiff_text:6d}"
              f"  dif. max: {maxdiff:.3e}" + (f"  en {where}" if where and maxdiff else ""))
    return ok


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
