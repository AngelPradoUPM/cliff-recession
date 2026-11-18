#!/usr/bin/env python3
"""Pruebas del modo mar y de los vuelcos del motor de la herramienta.

No hay un Fortran marino con el que comparar, asi que se comprueban propiedades:

1. Fracciones de tiempo de la marea sinusoidal: suman 1, son simetricas y el
   agua pasa mas tiempo cerca de pleamar y bajamar que en el nivel medio.
2. Reduccion al lago: en modo mar con amplitud 0, sin subida del nivel y con
   el oleaje del lago (invierno sin olas), el resultado es IDENTICO al del modo
   lago (y por tanto al del Fortran). Esto comprueba la reestructuracion del
   bucle mensual, incluidos los saltos por proteccion de derrubios.
3. Con marea, la erosion de la plataforma se reparte por el rango de marea.
4. Los vuelcos registrados tienen volumen positivo y el perfil "despues" no
   queda por delante del "antes".

Uso: python tools/test_mar.py
"""
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)
import great_lake_motor as gm  # noqa: E402


def ejecutar(caso, **kw):
    t = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(RAIZ, "ejemplos", caso, "data.dat"), t)
        return gm.run(t, verbose=False, **kw)
    finally:
        shutil.rmtree(t)


class TestMarea(unittest.TestCase):
    def test_fracciones(self):
        for amp, dz, n in ((3.0, 0.1, 15), (1.0, 0.1, 21), (0.5, 0.05, 81)):
            b = gm.bandas_marea(amp, dz, n)
            fr = [f for _, f in b]
            self.assertAlmostEqual(sum(fr), 1.0, places=12)
            self.assertTrue(all(abs(fr[i] - fr[-1 - i]) < 1e-12 for i in range(len(fr))))
            self.assertGreater(fr[0], fr[len(fr) // 2])
            self.assertLessEqual(len(b), n)

    def test_sin_marea(self):
        self.assertEqual(gm.bandas_marea(0.0, 0.1, 15), [(0, 1.0)])


class TestReduccionAlLago(unittest.TestCase):
    def test_mar_sin_marea_igual_que_lago(self):
        for caso in ("caso_simple", "caso_extremo"):
            lago = ejecutar(caso)
            mar = ejecutar(caso, mar={"amplitud": 0.0, "subida_mm_ano": 0.0},
                           seasons=gm.SEASONS)
            a, b = lago["DibMatriz"].copy(), mar["DibMatriz"].copy()
            a[1, 2] = b[1, 2] = 0.0   # tiempo de CPU
            self.assertTrue(np.array_equal(a, b), caso)
            self.assertEqual(len(lago["roturas"]), len(mar["roturas"]))
            self.assertEqual(len(lago["vuelcos"]), len(mar["vuelcos"]))


class TestMarConMarea(unittest.TestCase):
    def test_erosion_repartida_en_el_rango_de_marea(self):
        d = gm.read_data(os.path.join(RAIZ, "ejemplos", "caso_mar", "data.dat"))
        d["tsim"] = 3   # anos
        t = tempfile.mkdtemp()
        try:
            res = {}
            for amp in (0.0, 3.0):
                sub = os.path.join(t, f"a{amp:g}")
                os.makedirs(sub)
                gm.write_data(os.path.join(sub, "data.dat"), d)
                r = gm.run(sub, verbose=False, mar={"amplitud": amp, "subida_mm_ano": 0.0})
                P0 = r["profile_inicial"]
                N = r["meta"]["N"]
                res[amp] = (P0[1, 1:N + 1], P0[2, 1:N + 1] - r["DibMatriz"][2:N + 2, -1])
            z, ret0 = res[0.0]
            _, ret3 = res[3.0]
            # sin marea la plataforma solo se erosiona por debajo del nivel medio;
            # con marea tambien entre el nivel medio y la pleamar (0 < z < 2.5)
            banda = (z > 0.3) & (z < 2.3)
            self.assertAlmostEqual(float(np.max(ret0[banda])), 0.0, places=6)
            self.assertGreater(float(np.max(ret3[banda])), 0.05)
            self.assertTrue(np.all(np.isfinite(ret3)))
        finally:
            shutil.rmtree(t)

    def test_vuelcos_coherentes(self):
        r = ejecutar("caso_simple")
        self.assertGreater(len(r["vuelcos"]), 0)
        for mes, vol, antes, despues in r["vuelcos"]:
            self.assertGreater(vol, 0)
            self.assertTrue(np.all(despues <= antes + 1e-9))


if __name__ == "__main__":
    unittest.main(verbosity=2)
