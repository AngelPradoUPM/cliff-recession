#!/usr/bin/env python3
"""Pruebas de los metodos de estabilidad del motor de la herramienta.

Comprueba propiedades teoricas que deben cumplir los metodos de rebanadas:

1. Suelo puramente cohesivo (phi = 0): Fellenius, Bishop y Spencer dan
   exactamente el mismo FS, porque las fuerzas entre rebanadas no influyen.
2. Bishop simplificado es el caso particular de la ecuacion de momentos de
   Spencer con theta = 0.
3. La solucion de Spencer cumple a la vez el equilibrio de fuerzas y el de
   momentos (residuos nulos) con el mismo FS.
4. Con rozamiento, Bishop >= Fellenius, y Spencer queda muy cerca de Bishop
   (en superficies circulares suelen diferir menos de un 2-3 %).
5. Con las opciones por defecto, el motor de la herramienta reproduce
   exactamente a great_lake.py (comprobado tambien en tools/verify.sh).

Uso: python tools/test_estabilidad.py
"""
import math
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, RAIZ)
import great_lake_motor as gm  # noqa: E402


def rebanadas_talud(altura=10.0, angulo=45.0, xc=-2.0, zc=14.0, n=30, c=15.0, phi=25.0,
                    gamma=20.0, ru=0.0):
    """Rebanadas de un circulo (centro xc, zc) que corta un talud simple.

    Talud: terreno a cota `altura` para y < -altura/tan(angulo), cara inclinada
    hasta el pie en (0, 0) y fondo horizontal a cota 0 para y > 0. El circulo
    pasa por el pie. `ru` es el coeficiente de presion intersticial u = ru*gamma*h.
    """
    yc_cresta = -altura / math.tan(math.radians(angulo))
    r = math.hypot(xc - 0.0, zc - 0.0)

    def superficie(y):
        if y <= yc_cresta:
            return altura
        if y >= 0:
            return 0.0
        return altura * (-y) / (-yc_cresta)

    # extremo superior: interseccion del circulo con la coronacion
    y_sup = xc - math.sqrt(r ** 2 - (altura - zc) ** 2)
    reb = []
    ys = np.linspace(y_sup, 0.0, n + 1)
    for y0, y1 in zip(ys[:-1], ys[1:]):
        b = y1 - y0
        ym = 0.5 * (y0 + y1)
        zb = zc - math.sqrt(r ** 2 - (ym - xc) ** 2)
        h = superficie(ym) - zb
        if h <= 0:
            continue
        alfa = math.atan((xc - ym) / (zc - zb))   # positiva: empuja hacia el pie (+y)
        W = gamma * h * b
        u = ru * gamma * h
        reb.append((b, alfa, W, u, c, math.radians(phi)))
    return reb


class TestMetodos(unittest.TestCase):
    def test_cohesivo_puro_todos_iguales(self):
        reb = rebanadas_talud(c=30.0, phi=0.0)
        f = gm.fs_fellenius(reb)
        b = gm.fs_bishop(reb, f)
        s = gm.fs_spencer(reb, f)
        self.assertAlmostEqual(b, f, places=9)
        self.assertAlmostEqual(s, f, places=6)

    def test_bishop_es_spencer_con_theta_cero(self):
        for ru in (0.0, 0.3):
            reb = rebanadas_talud(ru=ru)
            b = gm.fs_bishop(reb, 1.0)
            fm0 = gm.spencer_fm(reb, 0.0)
            self.assertAlmostEqual(fm0, b, places=6)

    def test_spencer_cumple_ambos_equilibrios(self):
        for ru in (0.0, 0.25):
            reb = rebanadas_talud(ru=ru)
            f, th = gm.fs_spencer(reb, 1.0, devolver_theta=True)
            self.assertTrue(math.isfinite(f))
            q = gm._q_spencer(f, th, reb)
            fuerzas = sum(qi for qi, _ in q)
            momentos = sum(qi * math.cos(a - th) for qi, a in q)
            escala = sum(abs(qi) for qi, _ in q)
            self.assertLess(abs(fuerzas) / escala, 1e-5)
            self.assertLess(abs(momentos) / escala, 1e-5)

    def test_orden_tipico_de_resultados(self):
        for ru in (0.0, 0.2, 0.4):
            reb = rebanadas_talud(ru=ru)
            f = gm.fs_fellenius(reb)
            b = gm.fs_bishop(reb, f)
            s = gm.fs_spencer(reb, f)
            self.assertGreaterEqual(b, f - 1e-9)
            self.assertLess(abs(s - b) / b, 0.03)

    def test_fellenius_igual_a_formula_del_modelo(self):
        # La formula del modelo (Num/Deno) es la de fs_fellenius con u*l.
        reb = rebanadas_talud(ru=0.2, n=5)
        num = den = 0.0
        for b, a, W, u, c, phi in reb:
            num += c * (b / math.cos(a)) + (W * math.cos(a) - u * b / math.cos(a)) * math.tan(phi)
            den += W * math.sin(a)
        self.assertAlmostEqual(gm.fs_fellenius(reb), num / den, places=12)


class TestMotorPorDefecto(unittest.TestCase):
    def test_opciones_explicitas_por_defecto_igual_que_sin_opciones(self):
        caso = os.path.join(RAIZ, "ejemplos", "caso_simple", "data.dat")
        t = tempfile.mkdtemp()
        try:
            for sub in ("a", "b"):
                os.makedirs(os.path.join(t, sub))
                shutil.copy(caso, os.path.join(t, sub, "data.dat"))
            ra = gm.run(os.path.join(t, "a"), verbose=False)
            rb = gm.run(os.path.join(t, "b"), verbose=False,
                        estabilidad={"metodo": "fellenius", "rebanadas": 5, "criterio": "todas"})
            da, db = ra["DibMatriz"].copy(), rb["DibMatriz"].copy()
            da[1, 2] = db[1, 2] = 0.0   # tiempo de CPU
            self.assertTrue(np.array_equal(da, db))
        finally:
            shutil.rmtree(t)

    def test_metodos_alternativos_se_ejecutan(self):
        caso = os.path.join(RAIZ, "ejemplos", "caso_simple", "data.dat")
        t = tempfile.mkdtemp()
        try:
            for metodo in ("bishop", "spencer"):
                for criterio in ("todas", "profunda", "minimo"):
                    d = os.path.join(t, f"{metodo}_{criterio}")
                    os.makedirs(d)
                    shutil.copy(caso, os.path.join(d, "data.dat"))
                    r = gm.run(d, verbose=False,
                               estabilidad={"metodo": metodo, "rebanadas": 10,
                                            "criterio": criterio})
                    self.assertTrue(np.all(np.isfinite(r["DibMatriz"][2:, 3:])))
                    for fs in r["roturas_fs"]:
                        self.assertLess(fs, 1.0)
        finally:
            shutil.rmtree(t)


if __name__ == "__main__":
    unittest.main(verbosity=2)
