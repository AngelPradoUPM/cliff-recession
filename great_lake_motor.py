#!/usr/bin/env python3
"""
Great Lake - motor de calculo de la herramienta grafica.

Parte de ``great_lake.py`` (traduccion fiel del Fortran Great_lake_v.1) y
anade opciones que NO estan en el modelo original. Con las opciones por
defecto da exactamente los mismos resultados que ``great_lake.py`` (y que el
Fortran); ``tools/verify.sh`` lo comprueba.

Opciones adicionales de ``run``:

* ``nivel_lago``: variacion estacional e interanual del nivel del lago.
* ``estabilidad``: metodo de equilibrio limite para las roturas circulares
  ("fellenius" -original-, "bishop" simplificado o "spencer"), numero de
  rebanadas (original: 5) y criterio de rotura ("todas" -original-,
  "profunda" o "minimo").
* ``seasons``: tabla de oleaje por estacion.
* ``progress``: funcion de progreso / cancelacion.

Uso por linea de ordenes:
    python great_lake_motor.py DIRECTORIO [--nivel-estacional] [--nivel-interanual]
        [--metodo {fellenius,bishop,spencer}] [--rebanadas N]
        [--criterio {todas,profunda,minimo}]
"""

import math
import os
import sys
import time

import numpy as np

np.seterr(all="ignore")

F64 = np.float64
# En el original "real*8, parameter :: pi= 3.141592653589793" es un literal de
# precision simple (sin sufijo d0): pi vale 3.1415927410125732.
pi = F64(np.float32(3.141592653589793))


# ---------------------------------------------------------------------------
# Utilidades para reproducir la aritmetica de Fortran
# ---------------------------------------------------------------------------
def f32(x):
    """Constante real por defecto de Fortran (REAL*4) promocionada a REAL*8."""
    return F64(np.float32(x))


def idiv(a, b):
    """Division entera de Fortran (trunca hacia cero)."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b > 0) else -q


def fix(x):
    """INT() de Fortran: trunca hacia cero."""
    return int(x)


def anint(x):
    """ANINT() de Fortran: redondeo al entero mas proximo, mitades hacia fuera."""
    x = float(x)
    if math.isnan(x) or math.isinf(x):
        return F64(x)
    ax = abs(x)
    r = math.floor(ax)
    if ax - r >= 0.5:
        r += 1
    return F64(math.copysign(r, x))


# Funciones elementales de la libm de C (las mismas que usa gfortran). Las de
# NumPy usan implementaciones vectoriales propias que pueden diferir en 1 ulp.
# Los envoltorios devuelven Inf/NaN en lugar de lanzar excepciones (IEEE).
def exp(x):
    try:
        return F64(math.exp(x))
    except OverflowError:
        return F64(math.inf)


def sqrt(x):
    x = float(x)
    if x < 0 or math.isnan(x):
        return F64(math.nan)
    return F64(math.sqrt(x))


def _wrap(fn):
    def f(x):
        try:
            return F64(fn(x))
        except ValueError:
            return F64(math.nan)
    return f


tan = _wrap(math.tan)
atan = _wrap(math.atan)
cos = _wrap(math.cos)
sin = _wrap(math.sin)

# Constantes de precision simple usadas en el codigo original
G = F64(9.81)                    # 9.81d0
G981_S = f32(9.81)               # freshwater = 9.81 (REAL*4)
SQRT_G981 = F64(np.sqrt(np.float32(9.81)))  # 9.81**0.5 se evalua en REAL*4
C128 = F64(1.28)                 # 1.28d0
C005 = F64(0.05)                 # 0.05d0
C3_1416 = F64(3.1416)            # 3.1416d0
TWO_314 = F64(2.0) * f32(3.14)   # 2.d0*3.14 (3.14 es REAL*4)
TWO_31416 = 2 * C3_1416          # 2*3.1416d0
FOUR_31416 = 4 * C3_1416         # 2.d0*2.d0*3.1416d0


def ffmt(v, w, d):
    """Descriptor Fw.d de Fortran."""
    v = float(v)
    if math.isnan(v):
        return "NaN".rjust(w)
    if math.isinf(v):
        s = "Infinity" if v > 0 else "-Infinity"
        return s.rjust(w) if len(s) <= w else ("Inf" if v > 0 else "-Inf").rjust(w)
    s = f"{v:.{d}f}"
    if len(s) > w:
        if s.startswith("0."):
            s = s[1:]
        elif s.startswith("-0."):
            s = "-" + s[2:]
    if len(s) > w:
        return "*" * w
    return s.rjust(w)


# ---------------------------------------------------------------------------
# Condiciones de oleaje por estacion
# ---------------------------------------------------------------------------
SEASONS = {
    "spring": ([0.5, 1, 2, 3], [2.0, 3.5, 4.5, 5.5], [33.03, 36.33, 2.23, 0.1]),
    "summer": ([0.5, 1, 2, 3], [2.0, 3.5, 4.5, 6.5], [31.27, 41.73, 1.9, 0.33]),
    "fall": ([0.5, 1, 2, 3, 4, 5, 6], [2.0, 3.8, 4.5, 5.5, 6.5, 7.0, 7.5],
             [8.87, 36.27, 16.57, 9.40, 1.5, 0.2, 0.1]),
    "winter": ([0], [0], [0]),
}


def set_waves(season, seasons=None):
    """Devuelve Waves, H, T, L, WF, Wn (vectores con indice desde 1).

    ``seasons`` permite sustituir la tabla SEASONS (oleaje por estacion:
    alturas H [m], periodos T [s] y frecuencias WF [%]). Por defecto se usan
    los valores del codigo original.
    """
    Hs, Ts, WFs = (seasons or SEASONS)[season]
    waves = len(Hs)
    H = np.array([0.0] + [f32(v) for v in Hs], dtype=F64)
    T = np.array([0.0] + [f32(v) for v in Ts], dtype=F64)
    WF = np.array([0.0] + [f32(v) for v in WFs], dtype=F64)
    L = np.zeros(waves + 1, dtype=F64)
    Wn = np.zeros(waves + 1, dtype=F64)
    for i in range(1, waves + 1):
        L[i] = f32(1.56) * (T[i] ** 2)   # CERC (2002)
        Wn[i] = F64(60 * 60) / T[i]      # olas por hora
    return waves, H, T, L, WF, Wn


# ---------------------------------------------------------------------------
# Lectura de datos
# ---------------------------------------------------------------------------
class _Reader:
    """Lectura secuencial de registros, al estilo de READ(1,*)."""

    def __init__(self, path):
        with open(path, encoding="latin-1") as fh:
            self.lines = fh.read().splitlines()
        self.pos = 0

    def skip(self):
        self.pos += 1

    def values(self, n):
        vals = []
        while len(vals) < n:
            if self.pos >= len(self.lines):
                raise EOFError
            toks = self.lines[self.pos].replace(",", " ").split()
            self.pos += 1
            vals.extend(toks)
        return vals[:n]


def read_data(path):
    r = _Reader(path)
    d = {}
    r.skip()
    r.skip(); d["rn"] = int(float(r.values(1)[0]))
    r.skip(); d["Utsim"] = int(float(r.values(1)[0]))
    r.skip(); d["tsim"] = int(float(r.values(1)[0]))
    r.skip(); d["UTdib"] = int(float(r.values(1)[0]))
    r.skip(); d["Tdib"] = int(float(r.values(1)[0]))
    r.skip(); d["Dz"] = F64(r.values(1)[0])
    r.skip(); d["lakelevel"] = int(float(r.values(1)[0]))   # INTEGER en Fortran
    r.skip(); d["Qb"] = F64(r.values(1)[0])
    r.skip(); d["Gs"] = F64(r.values(1)[0])
    r.skip()
    v = r.values(2)
    d["groundwater"] = int(float(v[0]))
    d["profagua"] = F64(v[1])
    r.skip(); d["ndatosresis"] = int(float(r.values(1)[0]))
    r.skip(); d["SFmin"], d["Kss"], d["p"] = (F64(x) for x in r.values(3))
    r.skip(); d["WFmin"], d["MMM"] = (F64(x) for x in r.values(2))
    r.skip(); d["Rsc"], d["Cb"], d["KRs"], d["Bmat"] = (F64(x) for x in r.values(4))
    n = d["ndatosresis"]
    names = ["sigmac", "Shearstrength", "Ecohesion", "Ephi", "pesounit", "potestrato"]
    for nm in names:
        d[nm] = np.zeros(n + 1, dtype=F64)
    r.skip()
    for i in range(1, n + 1):
        vals = r.values(6)
        for nm, x in zip(names, vals):
            d[nm][i] = F64(x)
    r.skip(); d["Mf"] = F64(r.values(1)[0])
    r.skip()
    pts = []
    while r.pos < len(r.lines):
        toks = r.lines[r.pos].replace(",", " ").split()
        r.pos += 1
        if not toks:
            continue
        if len(toks) < 3:
            break
        pts.append([F64(t) for t in toks[:3]])
    d["Pzy"] = pts
    return d


def _num(x):
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return f"{x:.1f}"
    return repr(x)


def write_data(path, d):
    """Escribe un data.dat con el formato que espera el programa."""
    L = ["Data entry for Great_lake_v.", "Enter run number", str(int(d["rn"])),
         "Units of simulation time (1:months; 2:years)_integer", str(int(d["Utsim"])),
         "Simulation time", str(int(d["tsim"])),
         "Units per result output (1:months; 2:years)_integer", str(int(d["UTdib"])),
         "Result output time", str(int(d["Tdib"])),
         "Vertical discretization size (m)", _num(d["Dz"]),
         "Tidal level (m above/below[in negative] zero)", str(int(d["lakelevel"])),
         "Sand amount (m3/m)", _num(d["Qb"]),
         "Sand grain diameter (m)", _num(d["Gs"]),
         "Groundwater (0:no; 1:yes; 2: fully saturated); Groundwater depth landward (m)",
         f"{int(d['groundwater'])},{_num(d['profagua'])}",
         "Number of geological layers (from top to bottom)", str(int(d["ndatosresis"])),
         "Clay erosion treshold value [0.5-20 Pa]; Coef. for shear stress [Kss]; Calibration coef. [p]",
         f"{_num(d['SFmin'])} {_num(d['Kss'])} {_num(d['p'])}",
         "Threshold surf force [50-500 Kg/m^2]; Coef. for bluff erosion [MMM]",
         f"{_num(d['WFmin'])} {_num(d['MMM'])}",
         "Treshold run-up distance; Coef. for run-up [Cb]; Coef. for Rs to shear [KRs]; Beach material factor [Bmat]",
         f"{_num(d['Rsc'])} {_num(d['Cb'])} {_num(d['KRs'])} {_num(d['Bmat'])}",
         "UCS (KPa); Shear strength (KPa); Effective cohesion (KPa); Effective friction angle (deg); "
         "Unit weight (KN/m**3); Strata depth (m)"]
    names = ["sigmac", "Shearstrength", "Ecohesion", "Ephi", "pesounit", "potestrato"]
    for i in range(1, int(d["ndatosresis"]) + 1):
        L.append(" ".join(_num(d[nm][i]) for nm in names))
    L += ["Debris Mobility Coeff. (m**2/day)", _num(d["Mf"]),
          "Bluff initial profile (m) [1:Study point; 2:Talus wedge; 3:Seepage point; "
          "4:Talus base point => roturanotch]"]
    for a, b, c in d["Pzy"]:
        L.append(f"{_num(a)} {_num(b)} {_num(c)}")
    with open(path, "w", encoding="latin-1") as fh:
        fh.write("\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def prepare_profile(Pzy_points, Dz, lakelevel):
    """Discretiza el perfil inicial y localiza los puntos singulares.

    Devuelve un diccionario con el perfil discretizado ``P`` (P[1,k] cota,
    P[2,k] posicion, k = 1..N) y los indices de los puntos: ``TL`` (nivel del
    lago), ``zerolevel``, ``crestatalud``, ``seepagelevel``, ``roturanotch`` y
    ``label`` (puntos de estudio).
    """
    # Perfil inicial: Pzy[k] = (cota, distancia, etiqueta), k = 1..NG
    NG = len(Pzy_points)
    Pzy = np.zeros((4, NG + 1), dtype=F64)
    nlabel = 0
    for k, (a, b, c) in enumerate(Pzy_points, start=1):
        Pzy[1, k], Pzy[2, k], Pzy[3, k] = a, b, c
        if c == 1.0:
            nlabel += 1

    # Factor de orientacion
    Factor = 0
    if Pzy[2, 1] < Pzy[2, 2]:
        Factor = 1
    elif Pzy[2, 1] > Pzy[2, 2]:
        Factor = -1

    # Numero de puntos de discretizacion
    NPrefino = fix(abs(Pzy[1, 1] - Pzy[1, NG]) / Dz)
    N = NPrefino

    # P[1, k] = cota (z), P[2, k] = posicion horizontal (y)
    P = np.zeros((3, N + 1), dtype=F64)
    P[1, 1] = Pzy[1, 1]; P[2, 1] = Pzy[2, 1]
    P[1, N] = Pzy[1, NG]; P[2, N] = Pzy[2, NG]
    for k in range(2, N):
        P[1, k] = P[1, k - 1] - Dz
        for ib in range(2, NG + 1):
            if Pzy[1, ib] < P[1, k]:
                if Pzy[2, ib] == Pzy[2, ib - 1]:
                    P[2, k] = Pzy[2, ib]
                else:
                    P[2, k] = ((P[1, k] - Pzy[1, ib - 1]) / (Pzy[1, ib] - Pzy[1, ib - 1])
                               * (Pzy[2, ib] - Pzy[2, ib - 1]) + Pzy[2, ib - 1])
                break

    # 1: punto de estudio; 2: cresta del talud; 3: punto de rezume;
    # 4: base del talud => roturanotch
    label = [0] * 11
    ilabel = 0
    crestatalud = 0
    seepagelevel = 0
    roturanotch = 0
    for ib in range(2, NG + 1):
        tag = Pzy[3, ib]
        if tag in (1.0, 2.0, 3.0, 4.0):
            for k in range(1, N + 1):
                if P[1, k] < Pzy[1, ib]:
                    if tag == 1.0:
                        ilabel += 1
                        label[ilabel] = k
                    elif tag == 2.0:
                        crestatalud = k
                    elif tag == 3.0:
                        seepagelevel = k
                    else:
                        roturanotch = k - 1
                    break

    # Nivel cero
    for j in range(1, N + 1):
        if P[1, j] <= 0:
            break
    else:
        j = N + 1
    zerolevel = j - 1

    # Nivel del lago
    for j in range(1, N + 1):
        if P[1, j] <= (P[1, zerolevel] + lakelevel):
            break
    else:
        j = N + 1
    TL = j
    return dict(Pzy=Pzy, NG=NG, nlabel=nlabel, Factor=Factor, N=N, P=P,
                label=label, crestatalud=crestatalud, seepagelevel=seepagelevel,
                roturanotch=roturanotch, zerolevel=zerolevel, TL=TL)


class Cancelled(Exception):
    """Calculo interrumpido desde la funcion de progreso."""


# ---------------------------------------------------------------------------
# Factor de seguridad de una superficie circular (metodos de rebanadas)
# ---------------------------------------------------------------------------
# Cada rebanada: (b, alfa, W, u, c, phi) con b ancho [m], alfa inclinacion de la
# base [rad] (positiva si empuja hacia el lago), W peso [kN/m], u presion
# intersticial en la base [kPa], c cohesion efectiva [kPa], phi angulo de
# rozamiento efectivo [rad].
METODOS_ESTABILIDAD = ("fellenius", "bishop", "spencer")
CRITERIOS_ROTURA = ("todas", "profunda", "minimo")


def fs_fellenius(reb):
    """Metodo ordinario de rebanadas (Fellenius, 1936)."""
    num = den = 0.0
    for b, a, W, u, c, phi in reb:
        l = b / math.cos(a)
        num += c * l + (W * math.cos(a) - u * l) * math.tan(phi)
        den += W * math.sin(a)
    return num / den if den else math.inf


def fs_bishop(reb, f0=1.0, tol=1e-9, itmax=500):
    """Bishop simplificado (1955): equilibrio de momentos y de fuerzas verticales."""
    den = sum(W * math.sin(a) for _, a, W, _, _, _ in reb)
    if den <= 0:
        return math.inf
    f = f0 if (f0 and math.isfinite(f0) and f0 > 0) else 1.0
    for _ in range(itmax):
        num = 0.0
        for b, a, W, u, c, phi in reb:
            t = math.tan(phi)
            m = math.cos(a) + math.sin(a) * t / f
            if m <= 1e-6:
                return math.nan
            num += (c * b + (W - u * b) * t) / m
        nf = num / den
        if abs(nf - f) < tol * max(1.0, abs(nf)):
            return nf
        f = nf
    return f


def _q_spencer(f, th, reb):
    """Resultantes Q de las fuerzas entre rebanadas (Spencer, 1967)."""
    qs = []
    for b, a, W, u, c, phi in reb:
        t = math.tan(phi)
        ca = math.cos(a)
        num = c * b / (f * ca) + (W * ca - u * b / ca) * t / f - W * math.sin(a)
        den = math.cos(a - th) * (1.0 + math.tan(a - th) * t / f)
        if abs(den) < 1e-9:
            return None
        qs.append((num / den, a))
    return qs


def _raiz(g, a, b, n=80, cerca=1.0):
    """Raiz de g en [a, b].

    Busca cambios de signo en una malla logaritmica, los refina por biseccion
    y descarta los que son polos (donde g crece en vez de anularse). Si hay
    varias raices devuelve la mas proxima a `cerca`.
    """
    raices = []
    prev = None
    for x in np.geomspace(a, b, n):
        x = float(x)
        v = g(x)
        if v is None or not math.isfinite(v):
            prev = None
            continue
        if v == 0:
            raices.append(x)
        elif prev is not None and (v > 0) != (prev[1] > 0):
            lo, hi, glo = prev[0], x, prev[1]
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                gm_ = g(mid)
                if gm_ is None or not math.isfinite(gm_):
                    break
                if (gm_ > 0) == (glo > 0):
                    lo, glo = mid, gm_
                else:
                    hi = mid
            mid = 0.5 * (lo + hi)
            gm_ = g(mid)
            if gm_ is not None and math.isfinite(gm_) and \
                    abs(gm_) <= 1e-6 * max(abs(prev[1]), abs(v)):
                raices.append(mid)
        prev = (x, v)
    if not raices:
        return None
    if cerca is None or not math.isfinite(cerca):
        return raices[0]
    return min(raices, key=lambda r: abs(r - cerca))


def spencer_fm(reb, th, cerca=1.0):
    """Factor de seguridad de momentos de Spencer para una inclinacion theta."""
    def g(f):
        q = _q_spencer(f, th, reb)
        return None if q is None else sum(qi * math.cos(a - th) for qi, a in q)
    return _raiz(g, 0.02, 50.0, cerca=cerca)


def spencer_ff(reb, th, cerca=1.0):
    """Factor de seguridad de fuerzas de Spencer para una inclinacion theta."""
    def g(f):
        q = _q_spencer(f, th, reb)
        return None if q is None else sum(qi for qi, _ in q)
    return _raiz(g, 0.02, 50.0, cerca=cerca)


def _secante(g, x0, x1, tol=1e-10, itmax=60):
    """Raiz de g por el metodo de la secante. Devuelve None si no converge."""
    g0, g1 = g(x0), g(x1)
    for _ in range(itmax):
        if g0 is None or g1 is None or not (math.isfinite(g0) and math.isfinite(g1)):
            return None
        if g1 == g0:
            return x1 if g1 == 0 else None
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        if not math.isfinite(x2):
            return None
        if abs(x2 - x1) <= tol * max(1.0, abs(x2)):
            return x2
        x0, g0, x1, g1 = x1, g1, x2, g(x2)
    return None


def _spencer_rapido(reb, ref):
    """Spencer por secantes anidadas partiendo de Bishop (theta = 0)."""
    def gm_(f, th):
        q = _q_spencer(f, th, reb)
        return None if q is None else sum(qi * math.cos(a - th) for qi, a in q)

    def gf_(f, th):
        q = _q_spencer(f, th, reb)
        return None if q is None else sum(qi for qi, _ in q)

    def d(th):
        fm = _secante(lambda f: gm_(f, th), ref, ref * 1.01)
        ff = _secante(lambda f: gf_(f, th), ref, ref * 1.01)
        if fm is None or ff is None or fm <= 0 or ff <= 0:
            return None
        return fm - ff

    th = _secante(d, 0.0, 0.1, tol=1e-10, itmax=40)
    if th is None or abs(th) > math.radians(80):
        return None
    f = _secante(lambda f: gm_(f, th), ref, ref * 1.01)
    if f is None or f <= 0:
        return None
    # comprobacion: ambos equilibrios con el mismo FS
    q = _q_spencer(f, th, reb)
    if q is None:
        return None
    escala = sum(abs(qi) for qi, _ in q) or 1.0
    if (abs(sum(qi for qi, _ in q)) > 1e-6 * escala
            or abs(sum(qi * math.cos(a - th) for qi, a in q)) > 1e-6 * escala):
        return None
    return f, th


def fs_spencer(reb, f0=1.0, devolver_theta=False):
    """Spencer (1967): equilibrio completo de fuerzas y momentos.

    Busca la inclinacion theta (constante) de las fuerzas entre rebanadas para
    la que el factor de seguridad de momentos (Fm) y el de fuerzas (Ff)
    coinciden, partiendo de Bishop simplificado (Fm con theta = 0). Primero se
    resuelve por secantes (rapido); si no converge y el circulo es
    potencialmente inestable (Bishop < 1.3), con una busqueda en malla
    (robusta, lenta); si es claramente estable, devuelve Bishop. Si no hay
    solucion devuelve NaN.
    """
    ref = fs_bishop(reb, f0)
    if not math.isfinite(ref) or ref <= 0:
        ref = f0 if (f0 and math.isfinite(f0) and f0 > 0) else 1.0
    def aceptable(f):
        # En superficies circulares Spencer y Bishop difieren unas pocas
        # unidades por ciento; una "solucion" muy alejada de Bishop suele
        # corresponder a una inclinacion theta sin sentido fisico.
        return f is not None and math.isfinite(f) and f > 0 and abs(f / ref - 1.0) <= 0.15

    rapido = _spencer_rapido(reb, ref)
    if rapido is not None and aceptable(rapido[0]):
        f, th = rapido
        return (f, th) if devolver_theta else f
    if ref >= 1.3:
        # Circulo claramente estable: Spencer difiere de Bishop en pocas
        # unidades por ciento, asi que no cambia ninguna decision; se evita
        # la busqueda lenta y se devuelve Bishop.
        return (ref, math.nan) if devolver_theta else ref
    res = _fs_spencer_malla(reb, ref, True)
    if not aceptable(res[0]):
        res = (math.nan, math.nan)
    return res if devolver_theta else res[0]


def _fs_spencer_malla(reb, ref, devolver_theta=False):
    """Spencer con busqueda de theta en malla (respaldo robusto)."""
    def d(th):
        a = spencer_fm(reb, th, cerca=ref)
        b = spencer_ff(reb, th, cerca=ref)
        return None if (a is None or b is None) else a - b

    nan = (math.nan, math.nan) if devolver_theta else math.nan
    # malla de theta ordenada desde 0 hacia fuera: se toma la solucion mas cercana a 0
    grados = sorted(np.arange(-40.0, 60.0 + 1e-9, 2.5), key=abs)
    vals = {}
    for gdeg in grados:
        th = math.radians(float(gdeg))
        vals[gdeg] = d(th)
    ordenados = sorted(vals)
    soluciones = []
    for g0, g1 in zip(ordenados[:-1], ordenados[1:]):
        v0, v1 = vals[g0], vals[g1]
        if v0 is None or v1 is None:
            continue
        if v0 == 0:
            soluciones.append(math.radians(g0))
            continue
        if (v0 > 0) == (v1 > 0):
            continue
        lo, hi, dlo = math.radians(g0), math.radians(g1), v0
        ok = True
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            dm = d(mid)
            if dm is None:
                ok = False
                break
            if (dm > 0) == (dlo > 0):
                lo, dlo = mid, dm
            else:
                hi = mid
        if not ok:
            continue
        mid = 0.5 * (lo + hi)
        dm = d(mid)
        # se descartan los cambios de signo debidos a polos
        if dm is not None and abs(dm) <= 1e-6 * max(abs(v0), abs(v1)):
            soluciones.append(mid)
    if not soluciones:
        return nan
    th = min(soluciones, key=abs)
    f = spencer_fm(reb, th, cerca=ref)
    if f is None:
        return nan
    return (f, th) if devolver_theta else f


# Modo mar: valores por defecto de la tesis para la costa de Holderness
# (cap. 9.5.1): ola media en aguas profundas H0 = 0.98 m y periodo 7.66 s (boya
# de Hornsea), amplitud de marea 3 m, periodo 12.46 h y subida relativa del
# nivel del mar de 6 mm/ano. Es el oleaje que usa la erosion de Castedo, que no
# tiene umbrales.
SEASONS_MAR = {k: ([0.98], [7.66], [100.0]) for k in ("winter", "spring", "summer", "fall")}
# La erosion de Trenhaile solo erosiona por encima de unos umbrales, asi que con
# una unica ola media apenas actua: para ella se ofrece una distribucion
# ILUSTRATIVA con temporales cuya altura media anual coincide con esos 0.98 m.
SEASONS_MAR_TEMPORALES = {
    "winter": ([0.5, 1, 2, 3, 4], [5.0, 6.5, 7.66, 8.5, 9.0], [30, 30, 20, 8, 2]),
    "spring": ([0.5, 1, 2, 3, 4], [5.0, 6.5, 7.66, 8.5, 9.0], [40, 30, 15, 4, 1]),
    "summer": ([0.5, 1, 2, 3, 4], [5.0, 6.5, 7.66, 8.5, 9.0], [50, 30, 8, 2, 0.5]),
    "fall": ([0.5, 1, 2, 3, 4], [5.0, 6.5, 7.66, 8.5, 9.0], [35, 30, 18, 6, 1.5]),
}
MAR_TESIS = {"amplitud": 3.0, "subida_mm_ano": 6.0, "niveles": 15}
ESTACION_DEL_MES = ("winter", "winter", "winter", "spring", "spring", "spring",
                    "summer", "summer", "summer", "fall", "fall", "fall")


def bandas_marea(amplitud, dz, niveles_max=21):
    """Niveles de agua de una marea sinusoidal y fraccion de tiempo en cada uno.

    Divide el rango de marea [-A, +A] en bandas de anchura multiplo de dz
    (como mucho `niveles_max` bandas) y devuelve [(desplazamiento en indices
    del perfil, fraccion de tiempo)], con desplazamiento positivo = nivel mas
    alto que el medio. Para w(t) = A sin(2 pi t / T), la fraccion de tiempo con
    el agua entre z1 y z2 es (asin(z2/A) - asin(z1/A)) / pi.
    """
    if amplitud <= 0:
        return [(0, 1.0)]
    pasos = max(1, int(round(2 * amplitud / dz)))
    k = max(1, math.ceil(pasos / max(1, niveles_max - 1)))
    h = k * dz
    n = int(round(amplitud / h))
    out = []
    for b in range(-n, n + 1):
        z1 = max(-1.0, min(1.0, (b * h - h / 2) / amplitud))
        z2 = max(-1.0, min(1.0, (b * h + h / 2) / amplitud))
        fr = (math.asin(z2) - math.asin(z1)) / math.pi
        if fr > 0:
            out.append((b * k, fr))
    total = sum(fr for _, fr in out)
    return [(d, fr / total) for d, fr in out]


# ---------------------------------------------------------------------------
# Erosion marina de Castedo (tesis, cap. 5, 8 y 9)
# ---------------------------------------------------------------------------
# Para cada ciclo de marea T (ec. 8.1):
#   y(z, t+T) = y(z, t) - Hb^(13/4) Tb^(3/2) / (K sigma_c(z)) * tan(beta)
#                           * integral_0^T pw((w(t) - z) / Db) dt
# con Hb de Komar y Gaughan (1972), Db = Hb / 0.78, tan(beta) = (dy/dz)^-1 la
# pendiente local (limitada a `pendiente_max`), w(t) = A sen(2 pi t / T) la
# marea y pw la funcion de forma ajustada a los ensayos de Skafel y Bishop
# (1994): polinomio racional de Chebyshev de orden 5/6 (tabla 5.iv), con un
# maximo de 12,5 en Depth/Db = 0,3. La integral se evalua con la regla de
# Simpson y dt = 448,56 s (100 pasos por marea), como en el cap. 8.
# La tesis no fija del todo las unidades: aqui sigma_c va en MPa (columna UCS
# de los estratos, en kPa / 1000), pw en mm/h por unidad de pendiente y la
# integral en horas (mm), que se pasa a metros. K se calibra con datos reales.
PW_NUM = (-1.069838130676276, -0.01221236978255141, 0.4610744748787188,
          0.07844716807628178, 0.652841133940385, -0.08926720346045524)
PW_DEN = (1.0, -1.441825011828062, 1.303412347447682, -0.9727210591773651,
          0.5803220956908709, -0.2877147080827893, 0.1702945908127132)
PERIODO_MAREA_H = 12.46
EROSION_METODOS = ("trenhaile", "castedo")
# K = 102 (tesis, cap. 9, Holderness). Con esta convencion de unidades, el caso
# de Holderness de la tesis (ola media, marea 3 m, till de 182 kPa) da ~1,3-1,5
# m/ano, frente a 1,55 m/ano medidos. Pendiente maxima 60 grados (cap. 9) y
# derrubios: 70 % del volumen de un deslizamiento y 0 % de un vuelco (cap. 9).
CASTEDO_DEFECTO = {"K": 102.0, "pendiente_max": 60.0,
                   "derrubios_deslizamiento": 0.7, "derrubios_vuelco": 0.0}


def _chebyshev(coef, x):
    t0 = np.ones_like(x)
    t1 = x
    s_ = coef[0] * t0 + coef[1] * t1
    for c in coef[2:]:
        t0, t1 = t1, 2 * x * t1 - t0
        s_ = s_ + c * t1
    return s_


def pw_castedo(r):
    """Funcion de forma pw (>= 0) para r = profundidad bajo el agua / Db."""
    r = np.asarray(r, dtype=float)
    x = (r - 0.5111604219506374) / 0.5110522089614126
    with np.errstate(all="ignore"):
        v = -_chebyshev(PW_NUM, x) / _chebyshev(PW_DEN, x)
    return np.where((r >= 0.0) & (r <= 1.0), np.maximum(v, 0.0), 0.0)


def altura_rotura(h0, t0):
    """Altura de ola en rotura (Komar y Gaughan, 1972)."""
    return 0.39 * 9.81 ** 0.2 * (t0 * h0 ** 2) ** 0.4


def integral_marea(zrel, amplitud, db, pasos=100):
    """Integral de pw sobre un ciclo de marea (en horas), regla de Simpson.

    zrel: cotas respecto al nivel medio; amplitud de la marea (0 = nivel fijo).
    """
    zrel = np.asarray(zrel, dtype=float)
    if amplitud <= 0:
        return pw_castedo(-zrel / db) * PERIODO_MAREA_H
    t = np.linspace(0.0, PERIODO_MAREA_H, pasos + 1)
    w = amplitud * np.sin(2 * np.pi * t / PERIODO_MAREA_H)
    v = pw_castedo((w[:, None] - zrel[None, :]) / db)
    h = PERIODO_MAREA_H / pasos
    pesos = np.ones(pasos + 1)
    pesos[1:-1:2] = 4.0
    pesos[2:-1:2] = 2.0
    return (h / 3.0) * (pesos[:, None] * v).sum(axis=0)


NIVEL_ESTACIONAL_TESIS = {"spring": 0.4, "summer": 0.3, "fall": 0.0, "winter": 0.0}


def run(workdir=".", verbose=True, progress=None, seasons=None, nivel_lago=None,
        estabilidad=None, mar=None, erosion=None):
    """Ejecuta el modelo con el data.dat de ``workdir``.

    progress:   funcion opcional progress(Tglob, tsim); si devuelve False se
                interrumpe el calculo (excepcion Cancelled).
    seasons:    tabla de oleaje por estacion (ver SEASONS); None = original.
    nivel_lago: None (por defecto) = nivel del lago constante, como en el
                Fortran original. Si es un diccionario, activa la variacion del
                nivel que el original tiene comentada (Tabla 10.iii de la tesis):
                  "estacional": {"spring": m, "summer": m, "fall": m, "winter": m}
                                subida del nivel en cada estacion respecto a
                                lakelevel (tesis: 0.4, 0.3, 0, 0), o None;
                  "interanual": True para el ciclo de 15 anos (+-0.5 m con
                                Dz = 0.1, cambios de un Dz cada 9 meses a partir
                                del mes 69).
    estabilidad: None (por defecto) = como el original. Diccionario con:
                  "metodo":    "fellenius" (original), "bishop" o "spencer";
                  "rebanadas": numero de rebanadas (original: 5);
                  "criterio":  "todas" (original: rompe cada circulo con
                               FS < 1, de dentro hacia fuera), "profunda" (en
                               cada paso rompe el circulo inestable mas
                               profundo, como describe la tesis) o "minimo"
                               (el de menor FS).
    mar:        None (por defecto) = lago, como el original. Si es un
                diccionario, modo MAR (adaptacion, no esta en el original):
                  "amplitud":      amplitud de la marea sinusoidal (m);
                  "subida_mm_ano": subida del nivel medio del mar (mm/ano);
                  "niveles":       numero maximo de niveles de marea.
                En modo mar no hay meses de hielo (se calculan los 12 meses con
                el oleaje de su estacion), la erosion de cada ola se reparte
                entre los niveles de la marea segun el tiempo que el agua pasa
                en cada uno, y el nivel medio sube con el tiempo. Por defecto el
                oleaje es SEASONS_MAR (Holderness, tesis cap. 9).
    erosion:    None (por defecto) = erosion de Trenhaile del original. Si es
                un diccionario con "metodo": "castedo", usa la erosion marina de
                Castedo (tesis cap. 5, 8, 9) con "K" (calibracion),
                "pendiente_max" (grados) y las fracciones del volumen de cada
                rotura que quedan como derrubios al pie
                ("derrubios_deslizamiento", "derrubios_vuelco"; original: 1).
    """
    np.seterr(all="ignore")   # la configuracion de NumPy es por hilo
    tic = time.process_time()
    d = read_data(os.path.join(workdir, "data.dat"))

    rn = d["rn"]; Utsim = d["Utsim"]; tsim = d["tsim"]; UTdib = d["UTdib"]
    Tdib = d["Tdib"]; Dz = d["Dz"]; lakelevel = d["lakelevel"]; Qb = d["Qb"]
    Gs = d["Gs"]; groundwater = d["groundwater"]; Profagua = d["profagua"]
    SFmin = d["SFmin"]; Kss = d["Kss"]; p = d["p"]; WFmin = d["WFmin"]
    MMM = d["MMM"]; Rsc = d["Rsc"]; Cb = d["Cb"]; KRs = d["KRs"]; Bmat = d["Bmat"]
    Shearstrength = d["Shearstrength"]; Ecohesion = d["Ecohesion"]
    Ephi = d["Ephi"]; pesounit = d["pesounit"]; potestrato = d["potestrato"]
    ndatosresis = d["ndatosresis"]; Mf = d["Mf"]

    if verbose:
        print("Great_lake_v.1 - Bluff clay lakes profiles recesion")

    if Utsim == 2:
        tsim = tsim * 12
    if UTdib == 2:
        Tdib = Tdib * 12

    g = prepare_profile(d["Pzy"], Dz, lakelevel)
    Pzy = g["Pzy"]; NG = g["NG"]; nlabel = g["nlabel"]; Factor = g["Factor"]
    N = NPrefino = g["N"]; P = g["P"]; label = g["label"]
    crestatalud = g["crestatalud"]; seepagelevel = g["seepagelevel"]
    roturanotch = g["roturanotch"]; zerolevel = g["zerolevel"]; TL = g["TL"]
    P_inicial = P.copy()
    controllevel = TL
    if erosion and (erosion.get("metodo") == "castedo"):
        # resistencia a compresion (MPa) de cada punto del perfil segun su estrato
        sig_pt = np.zeros(N + 1)
        for k_ in range(1, N + 1):
            tp_ = 1
            acum = float(d["potestrato"][1])
            while P[1, 1] - acum > P[1, k_] and tp_ < int(d["ndatosresis"]):
                tp_ += 1
                acum += float(d["potestrato"][tp_])
            sig_pt[k_] = float(d["sigmac"][tp_]) / 1000.0
        if np.any(sig_pt[1:] <= 0):
            raise ValueError("La erosion de Castedo necesita la resistencia a compresion "
                             "(UCS) de todos los estratos")
    if mar is not None:
        _bandas = bandas_marea(float(mar.get("amplitud", MAR_TESIS["amplitud"])), float(Dz),
                               int(mar.get("niveles", MAR_TESIS["niveles"])))
        _subida = float(mar.get("subida_mm_ano", MAR_TESIS["subida_mm_ano"]))
        _tab = seasons if seasons is not None else SEASONS_MAR
        _tmax = max([float(t) for v in _tab.values() for t in v[1]] + [0.0])
        _db = fix(F64(f32(1.56) * _tmax ** 2 / 2) / Dz)
        _arriba = max(d for d, _ in _bandas) + int(_subida / 1000 * (tsim / 12) / float(Dz)) + 1
        _abajo = -min(d for d, _ in _bandas)
        if TL - _arriba < 2 or TL + _abajo + _db > N:
            raise ValueError(
                "El perfil no abarca el rango de la marea: debe llegar al menos "
                f"{(_abajo + _db) * float(Dz):.1f} m por debajo del nivel medio y "
                f"{_arriba * float(Dz):.1f} m por encima.")

    jdibujo = 1
    DimDiby = idiv(tsim, Tdib)

    # Matrices de resultados (indices desde 1)
    DibMatriz = np.zeros((N + 2, DimDiby + 3), dtype=F64)
    DibMatriz[1, 1] = rn
    DibMatriz[1, 2] = 0.0
    DibMatriz[2:N + 2, 1] = P[1, 1:N + 1]
    DibMatriz[2:N + 2, 2] = P[2, 1:N + 1]

    # Matrizrotura / Matrizrecesion / MatrizHwdraw: lista de columnas
    # (columna 0 ficticia); cada columna tiene indice de fila desde 1.
    r1 = np.zeros(N + 2, dtype=F64); r2 = np.zeros(N + 2, dtype=F64)
    r1[1] = rn; r2[1] = 0.0
    r1[2:N + 2] = P[1, 1:N + 1]   # columna 1: cotas
    r2[2:N + 2] = P[2, 1:N + 1]   # columna 2: posiciones
    Matrizrotura = [None, r1, r2]

    c1 = np.zeros(N + 3, dtype=F64); c2 = np.zeros(N + 3, dtype=F64)
    c1[1] = 1.0; c2[1] = 0.0
    c1[2] = 0.0; c2[2] = 0.0
    c1[3:N + 3] = P[1, 1:N + 1]; c2[3:N + 3] = P[2, 1:N + 1]
    Matrizrecesion = [None, c1, c2]

    freshwater = G981_S   # peso especifico del agua dulce (kN/m3)
    seawater = F64(1030)  # densidad (kg/m3)

    # Inicializacion de variables
    volacumulado = F64(0); volacumulado1 = F64(0); volacumulado2 = F64(0)
    changelakelevel = 1 if (nivel_lago and nivel_lago.get("interanual")) else 0
    estacional = (nivel_lago or {}).get("estacional") or {}
    desplaz = {k: int(round(float(estacional.get(k, 0.0)) / float(Dz)))
               for k in ("spring", "summer", "fall", "winter")}
    est_op = estabilidad or {}
    metodo = est_op.get("metodo", "fellenius")
    nreb = int(est_op.get("rebanadas", 5))
    criterio = est_op.get("criterio", "todas")
    if metodo not in METODOS_ESTABILIDAD:
        raise ValueError(f"metodo de estabilidad desconocido: {metodo}")
    if criterio not in CRITERIOS_ROTURA:
        raise ValueError(f"criterio de rotura desconocido: {criterio}")
    if nreb < 1:
        raise ValueError("el numero de rebanadas debe ser >= 1")
    MAX_ROTURAS_MES = 50   # salvaguarda para los criterios "profunda" y "minimo"
    ero_op = dict(erosion or {})
    metodo_erosion = ero_op.get("metodo", "trenhaile")
    if metodo_erosion not in EROSION_METODOS:
        raise ValueError(f"metodo de erosion desconocido: {metodo_erosion}")
    castedo = metodo_erosion == "castedo"
    if castedo:
        for k_, v_ in CASTEDO_DEFECTO.items():
            ero_op.setdefault(k_, v_)
    frac_desl = float(ero_op.get("derrubios_deslizamiento", 1.0)) if erosion else None
    frac_vuel = float(ero_op.get("derrubios_vuelco", 1.0)) if erosion else None
    modo_mar = mar is not None
    if modo_mar:
        if seasons is None:
            seasons = SEASONS_MAR
        amp_marea = float(mar.get("amplitud", MAR_TESIS["amplitud"]))
        subida_mar = float(mar.get("subida_mm_ano", MAR_TESIS["subida_mm_ano"]))
        bandas = bandas_marea(amp_marea, float(Dz), int(mar.get("niveles", MAR_TESIS["niveles"])))
    ciclo = 0
    Dprot = F64(0); Pmonth = 0; Drem = F64(0)
    parametrorotura = 0
    retroceso = F64(0.5)
    salto = 1
    counter = 1
    Draw13 = 0; Draw14 = 0
    invierno = 1
    pendientetgPcirc1 = tan(F64(1.0) * pi / F64(180.0))
    mcirculo = F64(0)
    nslides = F64(nreb) if estabilidad else F64(5)
    dimmatrizrotura = 0; dimmatrizrecesion = 0; dimmatrizHwdraw = 0

    # Variables que en Fortran conservan su valor entre iteraciones
    WBD = F64(0); HB = F64(0); Uo = F64(0); LL = F64(0)
    Playa = 0; DC = 0; BD = 0; DB = 0
    XDistPcirc2 = F64(0); YDistPcirc2 = F64(0); PtocorteCirc2 = 0
    YPtomedio = F64(0); DistPcirc2 = F64(0)
    CXc = F64(0); CYc = F64(0); radio = F64(0)
    LimWBmax = F64(0)

    distancia = np.zeros(1001, dtype=F64)
    Hwdraw = np.zeros(1001, dtype=F64)

    def groundwater_curve():
        Profseepage = abs(P[1, seepagelevel] - P[1, 1])
        Longitud = abs(P[2, seepagelevel] - P[2, 1])
        dist = Longitud / F64(1000.0)
        for i in range(1, 1001):
            distancia[i] = P[2, seepagelevel] - (i - 1) * dist
        for i in range(1, 1001):
            Hwdraw[i] = sqrt((abs(distancia[i] - P[2, seepagelevel]) / Longitud)
                             * (Profseepage - Profagua) ** 2) + P[1, seepagelevel]

    if groundwater == 1:
        groundwater_curve()

    h1 = np.zeros(1002, dtype=F64); h2 = np.zeros(1002, dtype=F64)
    h1[1] = rn; h2[1] = 0
    h1[2:1002] = Hwdraw[1:1001]; h2[2:1002] = distancia[1:1001]
    MatrizHwdraw = [None, h1, h2]

    # Condiciones de oleaje iniciales (invierno)
    Waves, H, T, L, WF, Wn = set_waves("winter", seasons)

    fout2 = open(os.path.join(workdir, "timerecesion.out"), "w")
    trec = []      # filas de timerecesion.out
    roturas = []   # (mes, volumen, perfil antes, perfil despues)
    roturas_fs = []      # factor de seguridad del circulo de cada rotura
    vuelcos = []         # (mes, volumen, perfil antes, perfil despues)
    n_fs_fallidos = 0    # circulos con FS < 1.5 en los que Spencer no converge (se usa Bishop)
    primer_paso = None   # (mes, perfil) al final del primer mes calculado
    meses_calc = []      # meses en los que se calcula erosion
    nivel_mes = []       # cota del nivel del lago en cada mes calculado
    perfiles_mes = []    # perfil al final de cada mes calculado
    meses_invierno = []  # meses sin calculo (invierno)
    meses_derrubios = [] # meses sin erosion por proteccion de derrubios

    def record(jt):
        """Guarda el perfil en DibMatriz si jt es un instante de salida."""
        nonlocal jdibujo
        if jt == jdibujo * Tdib and jdibujo * Tdib <= tsim:
            DibMatriz[1, jdibujo + 2] = jdibujo * Tdib
            DibMatriz[2:N + 2, jdibujo + 2] = P[2, 1:N + 1]
            jdibujo += 1

    # Ciclo interanual del nivel del lago (opcion; comentado en el original).
    # Cambios de un Dz cada 9 meses a partir del mes 60+9: 5 subidas, 10 bajadas
    # y 5 subidas (ciclo de 15 anos). El original solo comprobaba la fecha
    # exacta y perdia los cambios que caian en invierno o dentro de un salto
    # por derrubios, por lo que el nivel derivaba; aqui se aplican todos los
    # cambios pendientes hasta el mes t.
    n_cambios = 0

    def aplicar_cambios_nivel(t):
        nonlocal controllevel, n_cambios
        while 5 * 12 + 9 * (n_cambios + 1) <= t:
            n_cambios += 1
            n = (n_cambios - 1) % 20 + 1
            controllevel += -1 if (n <= 5 or n >= 16) else 1

    def change_level_exact(t):
        aplicar_cambios_nivel(t)

    def change_level_range(_t0, t1):
        aplicar_cambios_nivel(t1)

    Tglob = 1
    end_sim = False

    # ======================================================================
    # Bucle de erosion - recesion (etiqueta 500 en Fortran)
    # ======================================================================
    while Tglob <= tsim:
        if verbose:
            print(f"Time steep: {Tglob:8d}")
        if progress is not None and progress(Tglob, tsim) is False:
            fout2.close()
            raise Cancelled()
        goto = 1001

        if modo_mar:
            # oleaje de la estacion del calendario (no hay meses de hielo)
            Waves, H, T, L, WF, Wn = set_waves(ESTACION_DEL_MES[(Tglob - 1) % 12], seasons)
            mar_sin_oleaje = not np.any(T[1:] > 0)

        if invierno == 1 and not modo_mar:
            t_inv = Tglob
            salto = salto + counter
            if Tglob == 3 * (salto - 2) + 1:
                Tglob += 2
                Draw13 = 1; Draw14 = 0
            elif (3 * (salto - 2) + 1) < Tglob < 3 * (salto - 1):
                Tglob += 1
                Draw13 = 0; Draw14 = 1
            else:
                Draw13 = 0; Draw14 = 0
            if counter == 4:
                counter = 1
            meses_invierno.extend(range(t_inv, min(Tglob, tsim) + 1))
            Waves, H, T, L, WF, Wn = set_waves("spring", seasons)
            if nivel_lago:
                TL = controllevel - desplaz["spring"]

            last = tsim if Tglob >= tsim else Tglob
            if Draw13 == 1:
                for jt in range(Tglob - 2, last + 1):
                    record(jt)
            elif Draw14 == 1:
                for jt in range(Tglob - 1, last + 1):
                    record(jt)
            else:
                record(Tglob)
            if Tglob >= tsim:
                end_sim = True
                break
            goto = 1001
        elif modo_mar and mar_sin_oleaje:
            # estacion sin oleaje: como el invierno del lago, el mes se salta
            meses_invierno.append(Tglob)
            record(Tglob)
            goto = 1001
        else:
            # ------------------------------------------------------------
            # EROSION - Trenhaile (2009)
            # ------------------------------------------------------------
            BEACH = np.zeros(Waves + 1, dtype=F64)
            TER = np.zeros(N + 1, dtype=F64)
            TRunup = np.zeros(N + 1, dtype=F64)
            aux5 = P[2].copy()
            Ey = F64(0)

            if modo_mar:
                # nivel medio del mar (sube con el tiempo) y niveles de la marea
                TL = controllevel - int(subida_mar / 1000.0 * (Tglob / 12.0) / float(Dz))
                niveles_mes = [(TL - d, fr) for d, fr in bandas]
            else:
                niveles_mes = [(TL, 1.0)]
            TL_ref = TL
            ey_niveles = []
            if castedo:
                # erosion marina de Castedo: un paso por ciclo de marea
                niveles_mes = []
                amp_c = amp_marea if modo_mar else 0.0
                zrel = P[1, 1:N + 1] - P[1, TL]
                e_ciclo = np.zeros(N)
                db_max = 0.0
                for i in range(1, Waves + 1):
                    if T[i] > 0 and H[i] > 0 and WF[i] > 0:
                        hb_ = altura_rotura(float(H[i]), float(T[i]))
                        db_ = hb_ / 0.78
                        db_max = max(db_max, db_)
                        e_ciclo += ((float(WF[i]) / 100.0) * hb_ ** 3.25 * float(T[i]) ** 1.5
                                    * integral_marea(zrel, amp_c, db_))
                e_ciclo = e_ciclo / (ero_op["K"] * sig_pt[1:N + 1]) / 1000.0
                n_ciclos = int(round(float(30 - Dprot) * 24.0 / PERIODO_MAREA_H))
                tan_max = math.tan(math.radians(ero_op["pendiente_max"]))
                activo = e_ciclo > 0
                for _ in range(n_ciclos):
                    dydz = np.zeros(N)
                    dydz[:-1] = (P[2, 1:N] - P[2, 2:N + 1]) / float(Dz)
                    with np.errstate(divide="ignore"):
                        tanb = np.where(dydz != 0, 1.0 / np.abs(dydz), tan_max)
                    tanb = np.minimum(tanb, tan_max)
                    P[2, 1:N + 1] = np.where(activo, P[2, 1:N + 1] - e_ciclo * tanb,
                                             P[2, 1:N + 1])
                DB = int((amp_c + db_max) / float(Dz)) + 2
            for TL, f_marea in niveles_mes:
                Ey = F64(0)
                for i in range(1, Waves + 1):
                    # Altura de ola en asomeramiento
                    DA = L[i] / 2
                    DB = fix(DA / Dz)
                    DC = TL + DB
                    WH = np.zeros(N + 1, dtype=F64)

                    for j in range(DC, TL - 1, -1):
                        hL = (j - TL) * Dz / L[i]
                        goto62 = False
                        if hL > 0.25:
                            LL = L[i]
                            WH[j] = H[i]
                            goto62 = True
                        elif hL <= 0.25 and hL > C005:
                            DD = (exp(2 * pi * hL) - exp(-2 * pi * hL)) / (exp(2 * pi * hL) + exp(-2 * pi * hL))
                            LL = L[i] * sqrt(DD)
                        elif hL <= C005:
                            LL = T[i] * sqrt(G * (j - TL) * Dz)
                        if not goto62:
                            OO = 2 * 2 * pi * (j - TL) * Dz / LL
                            NN = (exp(OO) - exp(-OO)) / 2
                            n = F64(0.5) * (1 + (OO / NN))
                            if hL <= 0.25 and hL > C005:
                                xx = TWO_314 * (j - TL) * Dz / LL
                                DD1 = (exp(xx) - exp(-xx)) / (exp(xx) + exp(-xx))
                                WH[j] = H[i] * sqrt(L[i] / (2 * T[i] * n * (G * T[i] * DD1 / TWO_314)))
                            elif hL <= C005:
                                WH[j] = H[i] * sqrt(L[i] / (2 * T[i] * n * sqrt(G * (j - TL) * Dz)))
                        # 62
                        if WH[j] >= (j - TL) * Dz / C128:
                            WBD = int(WH[j] * C128 * F64(10.0)) / F64(10.0)
                            HB = WH[j]
                            break
                    # 63: longitud de la zona de rompientes
                    BD = fix(WBD / Dz) + TL
                    SZ = P[2, BD] - P[2, TL]

                    # Altura de ola en la zona de rompientes
                    for j in range(BD - 1, TL - 1, -1):
                        GR = F64(0)
                        if SZ != 0:
                            u1 = P[2, BD] - P[2, j]
                            u2 = P[1, BD] - P[1, j]
                            GR = abs(u2 / u1)
                        if P[2, j] == P[2, BD]:
                            WH[j] = HB
                        else:
                            WH[j] = HB * exp((F64(-35.2) * (abs(P[2, BD] - P[2, j])) * GR)
                                             / (T[i] * sqrt(G * WBD)))

                    # Pendiente de la playa - Sunamura (1989)
                    VVV = HB / (SQRT_G981 * (Gs ** F64(0.5)) * T[i])
                    BEACH[i] = atan(f32(0.12) / (VVV ** F64(0.5)))

                    # Forma y volumen de la playa
                    LimWB = anint((P[1, TL] + HB) * 10)   # REAL*8 en el original
                    LimWBmax = LimWB / 10
                    for j in range(1, BD + 1):
                        aux2 = anint(P[1, j] * 10)
                        aux2 = aux2 / 10
                        if aux2 <= LimWBmax:
                            break
                    else:
                        j = BD + 1
                    WBmax = j

                    Beachpoint = 0
                    BaseBeach = 0
                    Berm = 0

                    if Qb == 0:
                        Beachpoint = 0
                        BaseBeach = 0
                        Playa = 2
                    elif Qb != 0:
                        for j in range(BD, TL - 1, -1):
                            if P[2, j] == P[2, j - 1]:
                                Plat = pi / 2
                                Beachpoint = 0; BaseBeach = 0; Playa = 2
                                continue
                            else:
                                Plat = atan(Dz / (P[2, j] - P[2, j - 1]))
                            HWB = abs(P[1, WBmax] - P[1, j])
                            DWB = abs(P[2, WBmax] - P[2, j])
                            if BEACH[i] <= Plat:
                                Beachpoint = 0; BaseBeach = 0; Playa = 2
                                continue

                            Inter = 2
                            for k in range(j - 1, WBmax, -1):
                                u1 = P[2, k - 1] - P[2, k]
                                u2 = P[1, k - 1] - P[1, k]
                                mcliff = F64(0) if u1 == 0 else u2 / u1
                                mbeach = -atan(BEACH[i])
                                Yc = (((P[1, k] - mcliff * P[2, k]) - (P[1, j] - mbeach * P[2, j]))
                                      / (mbeach - mcliff))
                                Zc = mbeach * Yc + P[1, j] - mbeach * P[2, j]
                                if mcliff == 0:
                                    lam = (Zc - P[1, k]) / u2
                                else:
                                    lam = (Yc - P[2, k]) / u1
                                if lam >= 0.0 and lam < 1.0:
                                    Inter = k
                                    break
                                elif lam == 1.0:
                                    Inter = k - 1
                                    break

                            if (BEACH[i] * DWB < HWB) or (P[1, Inter] < P[1, WBmax]):
                                BE = P[2, Inter]
                                HWB = abs(P[1, Inter] - P[1, j])
                                DWB = abs(P[2, Inter] - P[2, j])
                                Beachpoint = Inter
                            elif (BEACH[i] * DWB == HWB) or (P[1, Inter] == P[1, WBmax]):
                                Beachpoint = WBmax
                                BE = P[2, WBmax]
                            elif (BEACH[i] * DWB > HWB) or (P[1, Inter] > P[1, WBmax]):
                                Beachpoint = WBmax
                                BE = P[2, j] - (HWB / tan(BEACH[i]))
                                Berm = 1

                            # Volumen
                            a = DWB * HWB
                            b = F64(0.5) * HWB * (P[2, j] - BE)
                            e = F64(0.5) * Dz * abs(P[2, j - 1] - P[2, j])
                            TB = F64(0)
                            CT = 0
                            M = 0
                            for k in range(j - 1, Beachpoint, -1):
                                if k < (M - 1):
                                    continue
                                elif k == (M - 1):
                                    M = 0
                                elif P[2, k] == P[2, k - 1]:
                                    jump29 = False
                                    for M in range(k + 1, Beachpoint - 1, -1):
                                        if P[2, M] != P[2, k]:
                                            jump29 = True
                                            break
                                        else:
                                            CT = M
                                    else:
                                        M = Beachpoint - 1
                                    if not jump29:
                                        continue
                                # 29
                                BX = (F64(0.5) * Dz * abs(P[2, k] - P[2, k - 1])
                                      + (j - k) * Dz * abs(P[2, k] - P[2, k - 1]))
                                TB = TB + BX

                            TOTAL = a - (b + e + TB)
                            BaseBeach = j
                            Playa = 1
                            if Qb >= TOTAL:
                                break
                            elif Qb < TOTAL:
                                Berm = 0
                    # 30
                    if BaseBeach == Beachpoint + 1:
                        Playa = 2
                        Berm = 0

                    # Numero de olas de cada tipo
                    No = Wn[i] * (WF[i] / 100) * 24 * (30 - Dprot)
                    if modo_mar:
                        No = No * f_marea   # tiempo que el agua pasa en este nivel

                    ER = np.zeros(N + 1, dtype=F64)
                    Runup = np.zeros(N + 1, dtype=F64)

                    # Erosion por esfuerzo cortante
                    Fw = f32(0.0034)
                    for j in range(DC, TL, -1):
                        Su = F64(0)
                        Sd = F64(0)
                        hL = (j - TL) * Dz / L[i]
                        if hL > 0.25:
                            ER[j] = 0
                            continue
                        elif hL <= 0.25 and hL > C005:
                            if j > BD:
                                XX = TWO_31416 * (j - TL) * Dz / L[i]
                                YY = (exp(XX) - exp(-XX)) / (exp(XX) + exp(-XX))
                                Lb = L[i] * sqrt(YY)
                                HH = FOUR_31416 * (j - TL) * Dz / Lb
                                II = (exp(HH) - exp(-HH)) / 2
                                Sd = C3_1416 * (WH[j] ** 2) / (4 * Lb * II)
                                GG = TWO_31416 * ((j - TL) * Dz - Sd) / Lb
                                EE = (exp(GG) - exp(-GG)) / 2
                                Uo = C3_1416 * WH[j] / (T[i] * EE)
                        elif hL <= C005:
                            Lb = T[i] * sqrt(G * (j - TL) * Dz)
                            HH = FOUR_31416 * (j - TL) * Dz / Lb
                            II = (exp(HH) - exp(-HH)) / 2
                            Sd = C3_1416 * (WH[j] ** 2) / (4 * Lb * II)
                            Uo = (WH[j] / 2) * sqrt(G / ((j - TL) * Dz - Sd))

                        if j < BD and P[2, j] < P[2, BD]:
                            Su = F64(0.185767) * (WBD / SZ) * (P[2, BD] - P[2, j])
                            Uo = (WH[j] / 2) * sqrt(G / ((j - TL) * Dz + Su))
                        elif j < BD and P[2, j] == P[2, BD]:
                            Su = F64(0)
                            Uo = (WH[j] / 2) * sqrt(G / ((j - TL) * Dz))

                        def shear():
                            SS = (1 / F64(2.0)) * Fw * seawater * (Uo ** 2)
                            return No * Kss * ((SS - SFmin) ** p)

                        if Playa == 1:
                            if TL < Beachpoint:
                                if j <= DC and j >= BaseBeach:
                                    ER[j] = shear()
                                elif j < Beachpoint and j > TL:
                                    ER[j] = shear()
                            elif TL >= Beachpoint:
                                if j <= DC and j >= BaseBeach:
                                    ER[j] = shear()
                        elif Playa == 2:
                            if j > TL:
                                ER[j] = shear()

                    # Erosion por impacto del oleaje + run-up
                    skip32 = False
                    if Playa == 1:
                        aux2 = anint(P[1, Beachpoint] * 10)
                        aux2 = aux2 / 10
                        if aux2 >= LimWBmax:
                            Runup[Beachpoint] = 0
                        elif aux2 < LimWBmax:
                            Rs = (LimWBmax - P[1, Beachpoint]) / tan(BEACH[i])
                            Runup[Beachpoint] = No * Cb * Bmat * (KRs * Rs - Rsc)
                        if TL > Beachpoint and TL < BaseBeach:
                            Ey = Ey + 0
                            skip32 = True
                        elif TL <= Beachpoint:
                            Runup[Beachpoint] = 0
                            if Berm == 1:
                                Ey = Ey + 0
                                skip32 = True
                    elif Playa == 2:
                        Runup[roturanotch] = 0
                    if not skip32:
                        if SZ == 0:
                            Fb = f32(512.61) * (HB / f32(0.78))
                        else:
                            Fb = f32(512.61) * (HB / f32(0.78)) * exp(-(SZ * f32(0.05)))
                        Energy = No * MMM * (Fb - WFmin)
                        if Energy < 0:
                            Energy = F64(0)
                        Ey = Ey + Energy
                    # 32: valores acumulados
                    for j in range(1, N + 1):
                        if ER[j] <= 0:
                            ER[j] = 0
                        if Runup[j] <= 0:
                            Runup[j] = 0
                        TER[j] = TER[j] + ER[j]
                        TRunup[j] = TRunup[j] + Runup[j]

                ey_niveles.append((TL, Ey))
            TL = TL_ref
            if modo_mar or castedo:
                DC = min(N - 1, TL + DB)

            # Nuevas posiciones tras la erosion
            for j in range(1, N + 1):
                P[2, j] = P[2, j] - TER[j] - TRunup[j]
            for tl_, ey_ in ey_niveles:
                P[2, tl_] = P[2, tl_] - ey_

            # --------------------------------------------------------
            # Puntos criticos
            # --------------------------------------------------------
            ptosup = crestatalud
            ptoinf = TL
            goto34 = False
            for j in range(crestatalud, TL):
                if P[2, j + 1] >= P[2, j]:
                    ptosup = j
                else:
                    ptosup = j + 1
                    break
            else:
                if ptosup == TL:
                    goto34 = True
            if not goto34:
                for j in range(TL, crestatalud, -1):
                    if P[2, j - 1] <= P[2, j]:
                        ptoinf = j - 1
                    else:
                        ptoinf = j
                        break
            # 34
            if ptosup == ptoinf:
                Upperpoint = ptoinf
                if ptosup == TL:
                    Upperpoint = roturanotch
                    if P[2, roturanotch] > P[2, ptosup]:
                        Upperpoint = ptoinf
            else:
                seg = P[2, ptosup:ptoinf + 1]
                if seg.size == 0:
                    Upperpoint = ptosup - 1      # MINLOC de un vector vacio = 0
                else:
                    Upperpoint = int(np.argmin(seg)) + ptosup

            # Punto critico inferior bajo el nivel del lago
            Lowerpoint = DC
            for j in range(DC, TL, -1):
                if P[2, j - 1] <= P[2, j]:
                    continue
                else:
                    Lowerpoint = j
                    break

            aux3 = P[2].copy()
            aux4 = P[2].copy()
            aux6 = P[2].copy()

            # Nuevo perfil si es necesario
            if Lowerpoint == DC:
                roturanotch = Upperpoint
            else:
                if Lowerpoint == Upperpoint and Lowerpoint <= roturanotch:
                    roturanotch = Lowerpoint
                elif Lowerpoint == Upperpoint and Lowerpoint > roturanotch:
                    if P[2, roturanotch] < P[2, Lowerpoint]:
                        pass
                    elif P[2, roturanotch] == P[2, Lowerpoint]:
                        roturanotch = Lowerpoint
                    elif P[2, roturanotch] > P[2, Lowerpoint]:
                        roturanotch = Lowerpoint
                elif Lowerpoint != Upperpoint:
                    if P[2, Upperpoint] < P[2, Lowerpoint]:
                        roturanotch = Upperpoint
                        mrecta = (P[1, Upperpoint] - P[1, Lowerpoint]) / (P[2, Upperpoint] - P[2, Lowerpoint])
                        for j in range(Lowerpoint, Upperpoint - 1, -1):
                            aux3[j] = (P[1, j] + mrecta * P[2, Lowerpoint] - P[1, Lowerpoint]) / mrecta
                    elif P[2, Upperpoint] == P[2, Lowerpoint]:
                        roturanotch = Lowerpoint
                        for j in range(Lowerpoint, Upperpoint - 1, -1):
                            aux3[j] = P[2, Lowerpoint]
                    elif P[2, Upperpoint] > P[2, Lowerpoint]:
                        roturanotch = Lowerpoint
                        mrecta = (P[1, Upperpoint] - P[1, Lowerpoint]) / (P[2, Upperpoint] - P[2, Lowerpoint])
                        for j in range(Lowerpoint, Upperpoint - 1, -1):
                            aux3[j] = (P[1, j] + mrecta * P[2, Lowerpoint] - P[1, Lowerpoint]) / mrecta
            # 38: correccion si el nuevo perfil "crea" material en la plataforma
            if Lowerpoint != Upperpoint:
                for j in range(Lowerpoint, Upperpoint - 1, -1):
                    if aux3[j] > P[2, j]:
                        mrecta = (P[1, zerolevel] - P[1, Lowerpoint]) / (P[2, zerolevel] - P[2, Lowerpoint])
                        if Lowerpoint == zerolevel:
                            mrecta = F64(1)
                        for ii in range(Lowerpoint, zerolevel - 1, -1):
                            P[2, ii] = (P[1, ii] + mrecta * P[2, Lowerpoint] - P[1, Lowerpoint]) / mrecta
                            if aux6[ii] < P[2, ii] and ii != zerolevel and ii != Lowerpoint:
                                # El bucle original sale en la primera iteracion
                                P[2, Lowerpoint] = aux6[Lowerpoint]
                                break
                        break
                    else:
                        P[2, j] = aux3[j]
            # 900
            # Si ptosup y ptoinf no coinciden, rotura entre ellos
            if Upperpoint != ptoinf and P[2, ptoinf] > P[2, Upperpoint]:
                mrecta = (P[1, Upperpoint] - P[1, ptoinf]) / (P[2, Upperpoint] - P[2, ptoinf])
                for ii in range(ptoinf, Upperpoint - 1, -1):
                    aux4[ii] = (P[1, ii] + mrecta * P[2, ptoinf] - P[1, ptoinf]) / mrecta
                Dibujo = 1
                if aux5[Upperpoint] == aux4[Upperpoint]:
                    Dibujo = 0
                    roturanotch = ptoinf
                else:
                    for j in range(ptoinf, Upperpoint - 1, -1):
                        if aux4[j] > P[2, j]:
                            Dibujo = 0
                            roturanotch = ptoinf
                            break
                    if Dibujo == 1:
                        for ii in range(ptoinf, Upperpoint - 1, -1):
                            P[2, ii] = aux4[ii]
            # 901

            # --------------------------------------------------------
            # ROTURA CIRCULAR (deslizamiento) - etiqueta 1
            # --------------------------------------------------------
            def evaluar_circulo(retro):
                """FS del circulo cuyo punto superior esta 'retro' m tierra adentro.

                Devuelve None si se supera el limite de busqueda (etiqueta 41).
                Mantiene, como el original, el valor de las variables del
                circulo entre llamadas.
                """
                nonlocal DistPcirc2, XDistPcirc2, YDistPcirc2, PtocorteCirc2
                nonlocal CXc, CYc, radio, n_fs_fallidos
                DistPcirc2 = P[2, roturanotch] - retro
                for i in range(roturanotch - 1, 1, -1):
                    u1 = P[2, i - 1] - P[2, i]
                    u2 = P[1, i - 1] - P[1, i]
                    if u1 == 0.0:
                        indvertcliff = 1
                        mcliff = F64(0.0)
                    else:
                        indvertcliff = 0
                        mcliff = u2 / u1
                    if mcliff == mcirculo:
                        continue
                    elif mcirculo == 0:
                        XDistPcirc2 = DistPcirc2
                        YDistPcirc2 = mcliff * XDistPcirc2 + P[1, i] - mcliff * P[2, i]
                    if indvertcliff == 0:
                        lam = (XDistPcirc2 - P[2, i]) / u1
                    else:
                        lam = (YDistPcirc2 - P[1, i]) / u2
                    if lam >= 0.0 and lam < 1.0:
                        PtocorteCirc2 = i
                        break
                    elif lam == 1.0:
                        PtocorteCirc2 = i - 1
                        break
                # 40
                Htalud = abs(P[1, crestatalud] - P[1, roturanotch])
                if DistPcirc2 < P[2, crestatalud] - Htalud:
                    return None  # -> 41

                munion = (YDistPcirc2 - P[1, roturanotch]) / (XDistPcirc2 - P[2, roturanotch])
                CXc = (((-P[2, roturanotch] / pendientetgPcirc1) - (P[1, roturanotch] / 2.0)
                        + (YDistPcirc2 / 2.0) - ((-1 / munion) * P[2, roturanotch] / 2.0)
                        - ((-1 / munion) * XDistPcirc2 / 2.0))
                       / ((-1 / pendientetgPcirc1) + (1 / munion)))
                CYc = (-1 / pendientetgPcirc1) * (CXc - P[2, roturanotch]) + P[1, roturanotch]
                radio = sqrt((CXc - P[2, roturanotch]) ** 2 + (CYc - P[1, roturanotch]) ** 2)
                anchoslide = abs(XDistPcirc2 - P[2, roturanotch]) / nslides

                Num = F64(0)
                Deno = F64(0)
                reb = []

                def lowy(X):
                    disc = sqrt(4 * CYc ** 2 + 4 * (radio ** 2 - CXc ** 2 - X ** 2 + 2 * CXc * X - CYc ** 2))
                    y1 = (2 * CYc + disc) / 2.0
                    y2 = (2 * CYc - disc) / 2.0
                    return y1 if y1 < y2 else y2

                for i in range(1, (nreb if estabilidad else 5) + 1):
                    Xca = P[2, roturanotch] - (i - 1) * anchoslide * Factor
                    Xcb = P[2, roturanotch] - i * anchoslide * Factor
                    Xmedio = (Xca + Xcb) / 2.0
                    Ymedio = lowy(Xmedio)
                    Mmedio = (CXc - Xmedio) / (Ymedio - CYc)
                    gamma = abs((180 / pi) * atan(Mmedio))

                    for j in range(roturanotch, 1, -1):
                        u1 = P[2, j - 1] - P[2, j]
                        u2 = P[1, j - 1] - P[1, j]
                        if u1 == 0.0:
                            indvertcliff = 1
                            mcliff = F64(0.0)
                        else:
                            indvertcliff = 0
                            mcliff = u2 / u1
                        XPtomedio = Xmedio
                        YPtomedio = mcliff * XPtomedio + P[1, j] - mcliff * P[2, j]
                        if indvertcliff == 0:
                            lam = (XPtomedio - P[2, j]) / u1
                        else:
                            lam = (YPtomedio - P[1, j]) / u2
                        if lam >= 0.0 and lam <= 1.0:
                            break
                    Hslide = abs(YPtomedio - Ymedio)

                    # Nivel freatico -> Hw
                    Hw = F64(0)
                    if groundwater == 1:
                        if Xmedio >= P[2, seepagelevel]:
                            Hw = Hslide
                        elif Xmedio < P[2, seepagelevel]:
                            Profseepage = abs(P[1, seepagelevel] - P[1, 1])
                            Longitud = abs(P[2, seepagelevel] - P[2, 1])
                            distY = abs(P[2, seepagelevel] - Xmedio)
                            curvah = sqrt((distY / Longitud) * (Profseepage - Profagua) ** 2)
                            Hw = abs(curvah + P[1, seepagelevel]) - Ymedio
                            if Hw <= 0:
                                Hw = F64(0)
                    elif groundwater == 0:
                        Hw = F64(0)
                    elif groundwater == 2:
                        Hw = Hslide

                    # Estrato geologico en la base de la rebanada
                    tp = 1
                    potacumulada = potestrato[tp]
                    while True:
                        if ((P[1, 1] - potacumulada) <= Ymedio
                                and Ymedio <= P[1, 1] - (potacumulada - potestrato[tp])):
                            Ecohesionrot = Ecohesion[tp]
                            Ephirot = Ephi[tp]
                            Pesounitrot = pesounit[tp]
                            break
                        tp += 1
                        if tp > ndatosresis:
                            raise RuntimeError("Ymedio fuera de los estratos definidos "
                                               "(en Fortran seria un acceso fuera de limites)")
                        potacumulada = potacumulada + potestrato[tp]

                    cg = cos(gamma * pi / 180)
                    Num = (Num + Ecohesionrot * (anchoslide / cg)
                           + (Hslide * Pesounitrot * anchoslide * cg
                              - Hw * freshwater * anchoslide / cg) * tan(Ephirot * pi / 180))
                    Deno = Deno + Hslide * Pesounitrot * anchoslide * sin(gamma * pi / 180)
                    if metodo != "fellenius":
                        alfa = math.atan(float((CXc - Xmedio) / (CYc - Ymedio))) * Factor
                        reb.append((float(anchoslide), alfa,
                                    float(Hslide * Pesounitrot * anchoslide),
                                    float(Hw * freshwater), float(Ecohesionrot),
                                    math.radians(float(Ephirot))))

                FS = Num / Deno
                if metodo != "fellenius":
                    f_b = fs_bishop(reb, float(FS))
                    f_alt = f_b if metodo == "bishop" else fs_spencer(reb, float(FS))
                    if not (math.isfinite(f_alt) and f_alt > 0):
                        # Spencer sin solucion: se usa Bishop (si tampoco, Fellenius)
                        f_alt = f_b
                        if math.isfinite(f_b) and f_b < 1.5:
                            n_fs_fallidos += 1
                    if math.isfinite(f_alt) and f_alt > 0:
                        FS = F64(f_alt)
                return FS

            pendiente = False      # criterio "profunda"/"minimo": circulo elegido
            r_sel = None
            descartados = set()
            n_sel = 0
            while True:
                if criterio != "todas" and not pendiente:
                    guard = (DistPcirc2, XDistPcirc2, YDistPcirc2, PtocorteCirc2,
                             CXc, CYc, radio)
                    cands = []
                    r = F64(0.5)
                    while True:
                        fs_r = evaluar_circulo(r)
                        if fs_r is None:
                            break
                        if fs_r < 1 and float(r) not in descartados:
                            cands.append((float(r), float(fs_r)))
                        r = r + F64(0.5)
                    (DistPcirc2, XDistPcirc2, YDistPcirc2, PtocorteCirc2,
                     CXc, CYc, radio) = guard
                    if not cands or n_sel >= MAX_ROTURAS_MES:
                        break  # -> 41
                    if criterio == "profunda":
                        r_sel = max(cands)[0]
                    else:
                        r_sel = min(cands, key=lambda cfs: cfs[1])[0]
                    retroceso = F64(r_sel)
                    pendiente = True
                    n_sel += 1

                FS = evaluar_circulo(retroceso)
                if FS is None:
                    break  # -> 41

                if FS >= 1:
                    retroceso = retroceso + F64(0.5)
                    if pendiente:
                        descartados.add(r_sel)
                        pendiente = False
                        retroceso = F64(0.5)
                    continue
                elif FS < 1:
                    retroceso = retroceso + F64(0.5)
                    if parametrorotura == 0:
                        dimmatrizrotura = len(Matrizrotura) - 1
                        c = np.zeros(N + 2, dtype=F64)
                        c[1] = Tglob
                        c[2:N + 2] = P[2, 1:N + 1]
                        Matrizrotura.append(c)

                    # Rotura circular desde roturanotch hasta PtocorteCirc2
                    for ii in range(roturanotch - 1, PtocorteCirc2 - 1, -1):
                        disc = sqrt(4 * CXc ** 2 + 4 * (radio ** 2 - CXc ** 2 - CYc ** 2
                                                        + 2 * CYc * P[1, ii] - (P[1, ii]) ** 2))
                        aux = (2 * CXc - disc) / 2.0
                        if aux <= P[2, ii]:
                            P[2, ii] = aux
                        elif aux > P[2, ii]:
                            aux = (2 * CXc + disc) / 2.0
                            P[2, ii] = aux

                    # Corregir geometria
                    restart = False
                    ref = Matrizrotura[dimmatrizrotura + 1]
                    for ii in range(roturanotch - 1, 1, -1):
                        if P[2, ii - 1] <= P[2, ii]:
                            pass
                        elif P[2, ii - 1] > P[2, ii]:
                            P[2, ii - 1] = P[2, ii]
                        else:
                            continue
                        if P[2, ii] > ref[ii + 1]:
                            parametrorotura = 0
                            P[2, 1:N + 1] = ref[2:N + 2]
                            del Matrizrotura[dimmatrizrotura + 1:]
                            restart = True
                            break
                    if restart:
                        if pendiente:
                            descartados.add(r_sel)
                            pendiente = False
                            retroceso = F64(0.5)
                        continue  # GO TO 1

                    # 48: guardar perfil tras la rotura
                    parametrorotura = 1
                    dimmatrizrotura = len(Matrizrotura) - 1
                    c = np.zeros(N + 2, dtype=F64)
                    c[1] = Tglob
                    c[2:N + 2] = P[2, 1:N + 1]
                    Matrizrotura.append(c)

                    # Matriz de recesion
                    dimmatrizrecesion = len(Matrizrecesion) - 1
                    c = np.zeros(N + 3, dtype=F64)
                    c[1] = Tglob
                    new = Matrizrotura[dimmatrizrotura + 1]
                    old = Matrizrotura[dimmatrizrotura]
                    for ii in range(2, N + 2):
                        c[ii + 1] = new[ii] - old[ii]
                    sumaincrementoy = F64(0)
                    for ii in range(1, N + 1):
                        sumaincrementoy = sumaincrementoy + abs(c[ii + 2])
                    c[2] = sumaincrementoy * Dz
                    Matrizrecesion.append(c)
                    volacumulado1 = volacumulado1 + sumaincrementoy * Dz
                    line = (ffmt(Tglob, 10, 1) + "  " + ffmt(sumaincrementoy * Dz, 9, 5) + "  "
                            + ffmt(volacumulado1, 9, 6) + "  "
                            + "".join(ffmt(abs(c[label[il] + 2]), 9, 6) + "  "
                                      for il in range(1, nlabel + 1)))
                    fout2.write(line.rstrip() + "\n")
                    roturas.append((Tglob, float(sumaincrementoy * Dz),
                                    old[2:N + 2].astype(float), new[2:N + 2].astype(float)))
                    roturas_fs.append(float(FS))
                    trec.append((Tglob, float(sumaincrementoy * Dz), float(volacumulado1),
                                 [float(abs(c[label[il] + 2])) for il in range(1, nlabel + 1)]))

                    # Matriz de niveles freaticos
                    if groundwater == 1:
                        groundwater_curve()
                        dimmatrizHwdraw = len(MatrizHwdraw) - 1
                        h1 = np.zeros(1002, dtype=F64); h2 = np.zeros(1002, dtype=F64)
                        h1[1] = Tglob; h2[1] = Tglob
                        h1[2:1002] = Hwdraw[1:1001]; h2[2:1002] = distancia[1:1001]
                        MatrizHwdraw.append(h1)
                        MatrizHwdraw.append(h2)
                    if pendiente:
                        pendiente = False
                        retroceso = F64(0.5)
                    continue  # GO TO 1
            # 41
            if parametrorotura == 1:
                volacumulado2 = F64(0)
                volacumulado = volacumulado1
                if frac_desl is not None:
                    volacumulado = volacumulado * frac_desl
            elif parametrorotura == 0:
                volacumulado1 = F64(0)
                # ROTURA POR VUELCO
                minestabilidad = F64(0)
                Ptocorte = 0
                for j in range(roturanotch - 1, 1, -1):
                    u1 = P[2, j - 1] - P[2, j]
                    u2 = P[1, j - 1] - P[1, j]
                    if u1 == 0.0:
                        indvertcliff = 1
                        mcliff = F64(0.0)
                    else:
                        indvertcliff = 0
                        mcliff = u2 / u1
                    if mcliff == minestabilidad:
                        continue
                    Yc = P[2, roturanotch]
                    Zc = mcliff * Yc + P[1, j] - mcliff * P[2, j]
                    if indvertcliff == 0:
                        lam = (Yc - P[2, j]) / u1
                    else:
                        lam = (Zc - P[1, j]) / u2
                    if lam >= 0.0 and lam < 1.0:
                        Ptocorte = j
                        break
                    elif lam == 1.0:
                        Ptocorte = j - 1
                        break
                if Ptocorte != 0:
                    vvr = P[2].copy()
                    nrot = 0
                    CoordX = F64(0)
                    for ii in range(roturanotch, Ptocorte - 1, -1):
                        vvr[ii] = P[2, roturanotch]
                        nrot += 1
                        CoordX = CoordX + abs(P[2, ii])
                    CoordX = CoordX / nrot
                    k = 1
                    potacumulada = potestrato[k]
                    while True:
                        if (P[1, 1] - potacumulada) <= P[1, roturanotch]:
                            Sigmatensilerock = Shearstrength[k]
                            Pesounitrot = pesounit[k]
                            break
                        k += 1
                        if k > ndatosresis:
                            raise RuntimeError("roturanotch fuera de los estratos definidos")
                        potacumulada = potacumulada + potestrato[k]
                    sumaacumuladarot = F64(0)
                    for ii in range(roturanotch, Ptocorte - 1, -1):
                        sumaacumuladarot = sumaacumuladarot + abs(vvr[ii] - P[2, ii])
                    volacumulado2 = sumaacumuladarot * Dz
                    sumaacumuladarot = volacumulado2 * Pesounitrot
                    hrot = abs(P[1, roturanotch] - P[1, Ptocorte])
                    Lcdg = abs(P[2, roturanotch]) - CoordX
                    Sigmatensilemax = (F64(6.0) * sumaacumuladarot * Lcdg) / (hrot ** 2)
                    if Sigmatensilemax >= Sigmatensilerock:
                        antes_v = P[2, 1:N + 1].astype(float)
                        for ii in range(Ptocorte, roturanotch + 1):
                            P[2, ii] = vvr[ii]
                        volacumulado = volacumulado2
                        if frac_vuel is not None:
                            volacumulado = volacumulado * frac_vuel
                        vuelcos.append((Tglob, float(volacumulado2), antes_v,
                                        P[2, 1:N + 1].astype(float)))
                    elif Sigmatensilemax < Sigmatensilerock:
                        volacumulado = F64(0)
            # 54
            if primer_paso is None:
                primer_paso = (Tglob, P[2, 1:N + 1].astype(float))
            meses_calc.append(Tglob)
            nivel_mes.append(float(P[1, TL]))
            perfiles_mes.append(P[2, 1:N + 1].astype(np.float32))
            record(Tglob)

            # Tiempo sin erosion por proteccion de los derrubios
            if volacumulado == 0:
                Pmonth = 0
                Dprot = F64(0)
                Drem = F64(0)
            else:
                Drem = volacumulado / Mf
                Pmonth = fix(Drem / 30.0)
                Dprot = abs(Pmonth * 30 - Drem)

            # Estacion del ano
            if modo_mar:
                # modo mar: calendario fijo; solo el salto por derrubios
                if Pmonth != 0:
                    Draw = Tglob
                    Tglob = Tglob + Pmonth
                    meses_derrubios.extend(range(Draw + 1, min(Tglob, tsim) + 1))
                    if Tglob >= tsim:
                        for jt in range(Draw, tsim + 1):
                            record(jt)
                        end_sim = True
                        break
                    for jt in range(Draw, Tglob + 1):
                        record(jt)
                goto = 1001
            elif Pmonth == 0:
                if changelakelevel == 1:
                    change_level_exact(Tglob)
                if Tglob == 3 * salto:
                    Waves, H, T, L, WF, Wn = set_waves("summer", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["summer"]
                    goto = 1001
                elif Tglob == 3 * (salto + 1):
                    Waves, H, T, L, WF, Wn = set_waves("fall", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["fall"]
                    goto = 1001
                elif Tglob == 3 * (salto + 2):
                    invierno = 1
                    counter = counter + 3
                    Waves, H, T, L, WF, Wn = set_waves("winter", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["winter"]
                    goto = 1002
            elif Pmonth != 0:
                Draw = Tglob
                Tglob = Tglob + Pmonth
                meses_derrubios.extend(range(Draw + 1, min(Tglob, tsim) + 1))
                if Tglob > Draw + 12:
                    yearjump = idiv(Tglob - Draw, 12)
                    salto = salto + 4 * yearjump
                if changelakelevel == 1:
                    change_level_range(Draw, Tglob)
                if Tglob >= tsim:
                    for jt in range(Draw, tsim + 1):
                        record(jt)
                    end_sim = True
                    break
                else:
                    for jt in range(Draw, Tglob + 1):
                        record(jt)
                if 3 * salto <= Tglob < 3 * (salto + 1):
                    Waves, H, T, L, WF, Wn = set_waves("summer", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["summer"]
                    goto = 1001
                elif 3 * (salto + 1) <= Tglob < 3 * (salto + 2):
                    Waves, H, T, L, WF, Wn = set_waves("fall", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["fall"]
                    goto = 1001
                elif 3 * (salto + 2) <= Tglob < 3 * (salto + 3):
                    invierno = 1
                    counter = counter + 3
                    Waves, H, T, L, WF, Wn = set_waves("winter", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["winter"]
                    goto = 1002
                elif Tglob == 3 * (salto + 3):
                    counter = 1
                    salto = salto + 4
                    Waves, H, T, L, WF, Wn = set_waves("spring", seasons)
                    if nivel_lago:
                        TL = controllevel - desplaz["spring"]
                    goto = 1001

        # 1001 / 1002
        if goto == 1001:
            invierno = 0
        Pmonth = 0
        Drem = F64(0)
        volacumulado1 = F64(0)
        volacumulado2 = F64(0)
        volacumulado = F64(0)
        parametrorotura = 0
        DistPcirc2 = F64(0)
        retroceso = F64(0.5)
        Tglob = Tglob + 1
    # 1000
    fout2.close()

    tac = time.process_time()
    if verbose:
        print("tiempo CPU:  ", tac - tic, "  seg.")
    DibMatriz[1, 2] = tac - tic

    def write_matrix(fname, cols, nrows):
        with open(os.path.join(workdir, fname), "w") as fh:
            for r in range(1, nrows + 1):
                vals = [ffmt(c[r], 15, 7) + " " for c in cols]
                # formato '(200(F15.7,x))': salto de registro cada 200 valores
                for s in range(0, len(vals), 200):
                    fh.write("".join(vals[s:s + 200]).rstrip() + "\n")

    write_matrix("Matriz_dibujo.out",
                 [DibMatriz[:, k] for k in range(1, DimDiby + 3)], N + 1)
    write_matrix("Matriz_roturas.out", Matrizrotura[1:], N + 1)
    write_matrix("Matriz_recesion.out", Matrizrecesion[1:], N + 2)
    write_matrix("Matriz_seepage.out", MatrizHwdraw[1:], 1001)

    return {
        "DibMatriz": DibMatriz,
        "Matrizrotura": Matrizrotura,
        "Matrizrecesion": Matrizrecesion,
        "MatrizHwdraw": MatrizHwdraw,
        "profile": P,
        "profile_inicial": P_inicial,
        "timerecesion": trec,
        "roturas": roturas,
        "roturas_fs": roturas_fs,
        "vuelcos": vuelcos,
        "modo": "mar" if modo_mar else "lago",
        "erosion": dict(ero_op, metodo=metodo_erosion),
        "mar": (dict(amplitud=amp_marea, subida_mm_ano=subida_mar, bandas=bandas)
                if modo_mar else None),
        "estabilidad": dict(metodo=metodo, rebanadas=(nreb if estabilidad else 5),
                            criterio=criterio, no_convergidos=n_fs_fallidos),
        "primer_paso": primer_paso,
        "mensual": dict(meses=np.array(meses_calc, dtype=int),
                        perfiles=(np.array(perfiles_mes) if perfiles_mes
                                  else np.zeros((0, N), dtype=np.float32)),
                        nivel=np.array(nivel_mes),
                        invierno=sorted(set(meses_invierno)),
                        derrubios=sorted(set(meses_derrubios))),
        "meta": dict(N=N, tsim=tsim, Tdib=Tdib, Dz=float(Dz), TL=TL,
                     zerolevel=zerolevel, crestatalud=crestatalud,
                     seepagelevel=seepagelevel, label=label[1:nlabel + 1],
                     groundwater=groundwater, profagua=float(Profagua)),
    }


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Modelo Great_lake_v.1 (recesion de acantilados)")
    ap.add_argument("directorio", nargs="?", default=".", help="carpeta con data.dat")
    ap.add_argument("--nivel-estacional", action="store_true",
                    help="nivel del lago variable por estaciones (Tabla 10.iii de la tesis)")
    ap.add_argument("--nivel-interanual", action="store_true",
                    help="ciclo interanual de 15 anos del nivel del lago")
    ap.add_argument("--metodo", choices=METODOS_ESTABILIDAD, default="fellenius",
                    help="metodo de estabilidad de las roturas circulares (original: fellenius)")
    ap.add_argument("--rebanadas", type=int, default=5,
                    help="numero de rebanadas (original: 5)")
    ap.add_argument("--criterio", choices=CRITERIOS_ROTURA, default="todas",
                    help="criterio de rotura (original: todas)")
    ap.add_argument("--erosion", choices=EROSION_METODOS, default="trenhaile",
                    help="erosion: trenhaile (original del lago) o castedo (marina, tesis)")
    ap.add_argument("--K", type=float, default=CASTEDO_DEFECTO["K"],
                    help="constante de calibracion de la erosion de Castedo")
    ap.add_argument("--mar", action="store_true",
                    help="modo mar: marea, sin hielo y subida del nivel del mar")
    ap.add_argument("--amplitud-marea", type=float, default=MAR_TESIS["amplitud"],
                    help="amplitud de la marea en m (modo mar; tesis Holderness: 3)")
    ap.add_argument("--subida-mar", type=float, default=MAR_TESIS["subida_mm_ano"],
                    help="subida del nivel del mar en mm/ano (modo mar; tesis: 6)")
    a = ap.parse_args()
    nl = None
    if a.nivel_estacional or a.nivel_interanual:
        nl = {"estacional": NIVEL_ESTACIONAL_TESIS if a.nivel_estacional else None,
              "interanual": a.nivel_interanual}
    est = None
    if (a.metodo, a.rebanadas, a.criterio) != ("fellenius", 5, "todas"):
        est = {"metodo": a.metodo, "rebanadas": a.rebanadas, "criterio": a.criterio}
    mr = ({"amplitud": a.amplitud_marea, "subida_mm_ano": a.subida_mar} if a.mar else None)
    er = {"metodo": "castedo", "K": a.K} if a.erosion == "castedo" else None
    run(a.directorio, nivel_lago=nl, estabilidad=est, mar=mr, erosion=er)
