#!/usr/bin/env python3
"""Genera una version del original Great_lake_v.1.f90 compilable con gfortran.

El original esta escrito para Intel/Compaq Visual Fortran y usa algunas
extensiones que gfortran no acepta. Este script hace solo los cambios minimos,
sin tocar la logica ni las constantes, para usarlo como referencia al
verificar la traduccion a Python:

  * real*8 x[allocatable](:)  ->  real*8, allocatable :: x(:)
  * parameter z=500,...        ->  integer, parameter :: z=500,...
  * se comenta "call seed(iseed)" (extension DEC; no se usa)
  * lakelevel (INTEGER) se lee a traves de un REAL, porque data.dat trae "0.0"
    y gfortran no acepta ese valor para un entero
  * write(2,11) Tglob -> dble(Tglob) (el formato 11 usa F para un entero)

Uso: python restore_fortran.py Great_lake_v.1.f90 salida.f90
"""
import re
import sys


def split_decl(decl):
    items, depth, cur = [], 0, ""
    for ch in decl:
        if ch == "," and depth == 0:
            items.append(cur.strip())
            cur = ""
            continue
        depth += ch == "("
        depth -= ch == ")"
        cur += ch
    items.append(cur.strip())
    return items


def main(src, dst):
    lines = open(src, encoding="latin-1").read().replace("\r", "").splitlines()
    out = []
    for line in lines:
        code, sep, com = line.partition("!")
        low = code.strip().lower()
        if low.startswith("real*8") and "[allocatable]" in code:
            for it in split_decl(code.strip()[len("real*8"):].strip()):
                if "[allocatable]" in it:
                    out.append("real*8, allocatable :: " + it.replace("[allocatable]", ""))
                else:
                    out.append("real*8 " + it)
            continue
        if re.match(r"\s*parameter\s+z\s*=", code, re.I):
            code = re.sub(r"parameter", "integer, parameter ::", code, count=1, flags=re.I)
        if re.match(r"\s*call\s+seed", code, re.I):
            code = "! " + code
        if re.match(r"\s*read\(1,\*\)\s*lakelevel\s*$", code, re.I):
            code = "read(1,*) rlakelevel; lakelevel = int(rlakelevel)"
        if re.match(r"\s*real\*8\s+Zc\s*$", code, re.I):
            code = code.rstrip() + ", rlakelevel"
        code = re.sub(r"write\(2,11\)\s*Tglob,", "write(2,11) dble(Tglob),", code, flags=re.I)
        out.append(code + (sep + com if sep else ""))
    open(dst, "w", encoding="latin-1").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
