#!/usr/bin/env python3
"""Genera una version compilable con gfortran de Great_lake_v.1.f90.

El fichero original tiene sintaxis mezclada (comentarios con %, "for" en lugar
de "do", "&&", round/fix, IF sin THEN, declaraciones [allocatable] de Compaq/
Intel Visual Fortran...). Este script lo devuelve a Fortran valido sin cambiar
la logica, para usarlo como referencia al verificar la traduccion a Python.

Cambios aplicados:
  * % -> !, for -> do, end for -> end do, && -> .and., anade THEN
  * round() -> anint(), fix() -> int()
  * real*8 x[allocatable](:) -> real*8, allocatable :: x(:)
  * Volacumulafor -> Volacumulado (variable sin declarar con IMPLICIT NONE;
    resto del reemplazo automatico do -> for)
  * se comenta "call seed(iseed)" (extension DEC, no se usa)
  * lakelevel (INTEGER) se lee a traves de un REAL, porque data.dat trae "0.0"
  * write(2,11) Tglob -> dble(Tglob) (formato F con un entero)

Uso: python restore_fortran.py Great_lake_v.1.f90 salida.f90
"""
import re, sys
src = open(sys.argv[1], encoding='latin-1').read().splitlines()
out = []
def close_paren(s, start):
    d = 0
    for i in range(start, len(s)):
        if s[i] == '(':
            d += 1
        elif s[i] == ')':
            d -= 1
            if d == 0:
                return i
    return -1
for line in src:
    line = line.replace('%', '!')
    code, sep, com = line.partition('!')
    code = code.replace('Volacumulafor', 'Volacumulado')
    code = re.sub(r'\bend\s+for\b', 'end do', code)
    code = re.sub(r'\bfor\b', 'do', code)
    code = code.replace('&&', '.and.')
    code = re.sub(r'\bround\(', 'anint(', code)
    code = re.sub(r'\bfix\(', 'int(', code)
    code = re.sub(r'^\s*parameter z=', 'integer, parameter :: z=', code)
    code = code.replace('real*8 Zc', 'real*8 Zc, rlakelevel')
    if re.match(r'\s*call seed', code):
        code = '! ' + code
    code = re.sub(r'^read\(1,\*\) lakelevel\s*$', 'read(1,*) rlakelevel; lakelevel = int(rlakelevel)', code)
    code = code.replace('write(2,11) Tglob,', 'write(2,11) dble(Tglob),')
    m = re.match(r'(\s*(?:\d+\s+)?)(if|elseif|else if|IF)\s*\(', code, re.I)
    if m:
        st = code.index('(', m.end(2))
        e = close_paren(code, st)
        if code[e+1:].strip() == '':
            code = code[:e+1] + ' then'
    # allocatable declarations (CVF syntax)
    if re.match(r'\s*real\*8\s', code) and '[allocatable]' in code:
        decl = code.strip()[len('real*8'):].strip()
        items, depth, cur = [], 0, ''
        for ch in decl:
            if ch == ',' and depth == 0:
                items.append(cur.strip()); cur = ''
                continue
            depth += ch == '('
            depth -= ch == ')'
            cur += ch
        items.append(cur.strip())
        lines = []
        for it in items:
            if '[allocatable]' in it:
                lines.append('real*8, allocatable :: ' + it.replace('[allocatable]', ''))
            else:
                lines.append('real*8 ' + it)
        out.extend(lines)
        continue
    out.append(code + (('!' + com) if sep else ''))
open(sys.argv[2], 'w').write('\n'.join(out) + '\n')
