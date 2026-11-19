#!/usr/bin/env python3
"""Pruebas de la erosion marina de Castedo (tesis, cap. 5, 8 y 9).

No existe el codigo original de la variante marina, asi que se comprueba la
implementacion frente a lo que publica la tesis:

1. Funcion de forma pw (tabla 5.iv): nula fuera de 0 <= profundidad/Db <= 1 y
   con un maximo de ~12,5 hacia profundidad/Db = 0,3 (figura 5.3).
2. Altura y profundidad de rotura (Komar y Gaughan): en el ejemplo del
   apartado 5.3.2 (H0 = 1,5 m, T = 1 s, marea 1,46 m) la profundidad Am + Db
   vale 2,55 m (linea discontinua de la figura 5.4).
3. La regla de Simpson con 100 pasos por marea (dt = 448,56 s, cap. 8)
   coincide con una integracion muy fina (a menos de un 5 % del maximo).
4. La erosion de un ciclo de marea solo actua entre la pleamar y Am + Db por
   debajo del nivel medio.
5. Configuracion de Holderness de la tesis (ola media 0,98 m y 7,66 s, marea
   3 m, K = 102, till de 182 kPa): retroceso del orden del medido, 1,55 m/ano.
6. Mas K (mas resistencia) => menos retroceso.

Uso: python tools/test_castedo.py
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

OLA_MEDIA = gm.SEASONS_MAR


def retroceso(K=None, anos=20):
    d = gm.read_data(os.path.join(RAIZ, "ejemplos", "caso_mar", "data.dat"))
    d["tsim"] = anos
    t = tempfile.mkdtemp()
    try:
        gm.write_data(os.path.join(t, "data.dat"), d)
        ero = {"metodo": "castedo"}
        if K is not None:
            ero["K"] = K
        r = gm.run(t, verbose=False, mar=dict(gm.MAR_TESIS), erosion=ero, seasons=OLA_MEDIA)
        P0 = r["profile_inicial"]
        N = r["meta"]["N"]
        z = P0[1, 1:N + 1]
        ret = P0[2, 1:N + 1] - r["DibMatriz"][2:N + 2, -1]
        k8 = int(np.argmin(np.abs(z - 8.0)))
        return ret[k8] / anos, r
    finally:
        shutil.rmtree(t)


class TestFormulas(unittest.TestCase):
    def test_funcion_de_forma(self):
        r = np.linspace(-0.5, 1.5, 4001)
        v = gm.pw_castedo(r)
        self.assertTrue(np.all(v[(r < 0) | (r > 1)] == 0))
        self.assertTrue(np.all(v >= 0))
        self.assertAlmostEqual(float(v.max()), 12.53, delta=0.05)
        self.assertTrue(0.25 <= float(r[np.argmax(v)]) <= 0.35)

    def test_rotura_ejemplo_tesis(self):
        hb = gm.altura_rotura(1.5, 1.0)
        db = hb / 0.78
        self.assertAlmostEqual(1.46 + db, 2.55, delta=0.01)

    def test_simpson_100_pasos(self):
        zrel = np.linspace(-4.5, 3.5, 81)
        a = gm.integral_marea(zrel, 3.0, 1.75, pasos=100)
        t = np.linspace(0, gm.PERIODO_MAREA_H, 200001)
        w = 3.0 * np.sin(2 * np.pi * t / gm.PERIODO_MAREA_H)
        fino = np.array([np.trapezoid(gm.pw_castedo((w - z) / 1.75), t) for z in zrel])
        # pw cae bruscamente a 0 en profundidad = Db, y ahi las reglas de
        # integracion pierden precision: con el paso elegido en la tesis
        # (dt = 448,56 s) la diferencia es de ~3 % del maximo.
        self.assertLess(np.max(np.abs(a - fino)) / np.max(fino), 0.05)

    def test_zona_erosionada(self):
        am, db = 1.46, gm.altura_rotura(1.5, 1.0) / 0.78
        zrel = np.linspace(-4, 3, 701)
        i = gm.integral_marea(zrel, am, db)
        self.assertTrue(np.all(i[zrel > am + 1e-9] == 0))
        self.assertTrue(np.all(i[zrel < -(am + db) - 0.02] == 0))
        self.assertTrue(np.all(i[(zrel > -am) & (zrel < am - 0.1)] > 0))

    def test_sin_marea(self):
        zrel = np.array([-1.0, -0.4, 0.2])
        np.testing.assert_allclose(gm.integral_marea(zrel, 0.0, 1.0),
                                   gm.pw_castedo(-zrel) * gm.PERIODO_MAREA_H)


class TestHolderness(unittest.TestCase):
    def test_retroceso_del_orden_del_medido(self):
        tasa, r = retroceso()
        self.assertEqual(r["erosion"]["K"], 102.0)
        self.assertTrue(1.0 <= tasa <= 2.2, f"{tasa:.2f} m/ano")

    def test_mas_K_menos_retroceso(self):
        t1, _ = retroceso(K=100.0, anos=5)
        t2, _ = retroceso(K=400.0, anos=5)
        self.assertGreater(t1, t2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
