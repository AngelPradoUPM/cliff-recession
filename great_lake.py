#!/usr/bin/env python3
"""
Great_lake_v.1 - Bluff clay lakes profiles recession (LAKE ERIE)

Traduccion a Python del programa Fortran ``Great_lake_v.1.f90``
(R. Castedo, UPM, 2011). Erosion de acantilados de arcillas (tills) en
entornos sin marea.

La traduccion es deliberadamente literal para reproducir los mismos
resultados que el ejecutable Fortran:

* Los vectores/matrices se indexan desde 1 (la posicion 0 se ignora), igual
  que en Fortran, para conservar los indices del original.
* Se respeta la semantica de los bucles DO de Fortran (valor de la variable
  de control al salir del bucle) y de los GO TO.
* Se respeta la aritmetica de Fortran: division entera truncada hacia cero
  (p.ej. ``tsim/Tdib``) y las constantes escritas sin ``d0`` en el original,
  que son de precision simple (``pi``, ``1.56``, frecuencias de oleaje,
  ``freshwater = 9.81``, ``9.81**0.5``, ``0.12``, ``0.0034``, ``512.61``,
  ``0.78``...) y se promocionan a doble precision.
* Las funciones elementales (exp, sqrt, tan...) son las de la libm de C, como
  en gfortran.
* Las operaciones en coma flotante siguen IEEE (divisiones por cero dan
  Inf/NaN en lugar de lanzar excepciones), igual que el ejecutable Fortran.

Uso:
    python great_lake.py [directorio_de_trabajo]

Lee ``data.dat`` del directorio de trabajo y escribe en el mismo directorio
``timerecesion.out``, ``Matriz_dibujo.out``, ``Matriz_roturas.out``,
``Matriz_recesion.out`` y ``Matriz_seepage.out``.
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


def set_waves(season):
    """Devuelve Waves, H, T, L, WF, Wn (vectores con indice desde 1)."""
    Hs, Ts, WFs = SEASONS[season]
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


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def run(workdir=".", verbose=True):
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

    # Perfil inicial: Pzy[k] = (cota, distancia, etiqueta), k = 1..NG
    NG = len(d["Pzy"])
    Pzy = np.zeros((4, NG + 1), dtype=F64)
    nlabel = 0
    for k, (a, b, c) in enumerate(d["Pzy"], start=1):
        Pzy[1, k], Pzy[2, k], Pzy[3, k] = a, b, c
        if c == 1.0:
            nlabel += 1

    if Utsim == 2:
        tsim = tsim * 12
    if UTdib == 2:
        Tdib = Tdib * 12

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
    controllevel = TL

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
    changelakelevel = 0
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
    nslides = F64(5)
    dimmatrizrotura = 0; dimmatrizrecesion = 0; dimmatrizHwdraw = 0

    # Variables que en Fortran conservan su valor entre iteraciones
    WBD = F64(0); HB = F64(0); Uo = F64(0); LL = F64(0)
    Playa = 0; DC = 0; BD = 0
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
    Waves, H, T, L, WF, Wn = set_waves("winter")

    fout2 = open(os.path.join(workdir, "timerecesion.out"), "w")

    def record(jt):
        """Guarda el perfil en DibMatriz si jt es un instante de salida."""
        nonlocal jdibujo
        if jt == jdibujo * Tdib and jdibujo * Tdib <= tsim:
            DibMatriz[1, jdibujo + 2] = jdibujo * Tdib
            DibMatriz[2:N + 2, jdibujo + 2] = P[2, 1:N + 1]
            jdibujo += 1

    def change_level_exact(t):
        """Cambio interanual del nivel del lago (changelakelevel == 1)."""
        nonlocal controllevel, ciclo
        base = 5 * 12 + ciclo
        for n in range(1, 21):
            if t == base + 9 * n:
                controllevel += -1 if (n <= 5 or n >= 16) else 1
                if n == 20:
                    ciclo += 180
                return

    def change_level_range(t0, t1):
        nonlocal controllevel, ciclo
        base = 5 * 12 + ciclo
        for n in range(1, 21):
            if t0 <= base + 9 * n <= t1:
                controllevel += -1 if (n <= 5 or n >= 16) else 1
                if n == 20:
                    ciclo += 180
                return

    Tglob = 1
    end_sim = False

    # ======================================================================
    # Bucle de erosion - recesion (etiqueta 500 en Fortran)
    # ======================================================================
    while Tglob <= tsim:
        if verbose:
            print(f"Time steep: {Tglob:8d}")
        goto = 1001

        if invierno == 1:
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
            Waves, H, T, L, WF, Wn = set_waves("spring")

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
        else:
            # ------------------------------------------------------------
            # EROSION - Trenhaile (2009)
            # ------------------------------------------------------------
            BEACH = np.zeros(Waves + 1, dtype=F64)
            TER = np.zeros(N + 1, dtype=F64)
            TRunup = np.zeros(N + 1, dtype=F64)
            aux5 = P[2].copy()
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

            # Nuevas posiciones tras la erosion
            for j in range(1, N + 1):
                P[2, j] = P[2, j] - TER[j] - TRunup[j]
            P[2, TL] = P[2, TL] - Ey

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
            while True:
                DistPcirc2 = P[2, roturanotch] - retroceso
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
                    break  # -> 41

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

                def lowy(X):
                    disc = sqrt(4 * CYc ** 2 + 4 * (radio ** 2 - CXc ** 2 - X ** 2 + 2 * CXc * X - CYc ** 2))
                    y1 = (2 * CYc + disc) / 2.0
                    y2 = (2 * CYc - disc) / 2.0
                    return y1 if y1 < y2 else y2

                for i in range(1, 6):
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

                FS = Num / Deno

                if FS >= 1:
                    retroceso = retroceso + F64(0.5)
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

                    # Matriz de niveles freaticos
                    if groundwater == 1:
                        groundwater_curve()
                        dimmatrizHwdraw = len(MatrizHwdraw) - 1
                        h1 = np.zeros(1002, dtype=F64); h2 = np.zeros(1002, dtype=F64)
                        h1[1] = Tglob; h2[1] = Tglob
                        h1[2:1002] = Hwdraw[1:1001]; h2[2:1002] = distancia[1:1001]
                        MatrizHwdraw.append(h1)
                        MatrizHwdraw.append(h2)
                    continue  # GO TO 1
            # 41
            if parametrorotura == 1:
                volacumulado2 = F64(0)
                volacumulado = volacumulado1
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
                        for ii in range(Ptocorte, roturanotch + 1):
                            P[2, ii] = vvr[ii]
                        volacumulado = volacumulado2
                    elif Sigmatensilemax < Sigmatensilerock:
                        volacumulado = F64(0)
            # 54
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
            if Pmonth == 0:
                if changelakelevel == 1:
                    change_level_exact(Tglob)
                if Tglob == 3 * salto:
                    Waves, H, T, L, WF, Wn = set_waves("summer")
                    goto = 1001
                elif Tglob == 3 * (salto + 1):
                    Waves, H, T, L, WF, Wn = set_waves("fall")
                    goto = 1001
                elif Tglob == 3 * (salto + 2):
                    invierno = 1
                    counter = counter + 3
                    Waves, H, T, L, WF, Wn = set_waves("winter")
                    goto = 1002
            elif Pmonth != 0:
                Draw = Tglob
                Tglob = Tglob + Pmonth
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
                    Waves, H, T, L, WF, Wn = set_waves("summer")
                    goto = 1001
                elif 3 * (salto + 1) <= Tglob < 3 * (salto + 2):
                    Waves, H, T, L, WF, Wn = set_waves("fall")
                    goto = 1001
                elif 3 * (salto + 2) <= Tglob < 3 * (salto + 3):
                    invierno = 1
                    counter = counter + 3
                    Waves, H, T, L, WF, Wn = set_waves("winter")
                    goto = 1002
                elif Tglob == 3 * (salto + 3):
                    counter = 1
                    salto = salto + 4
                    Waves, H, T, L, WF, Wn = set_waves("spring")
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
    }


if __name__ == "__main__":
    wd = sys.argv[1] if len(sys.argv) > 1 else "."
    if not os.path.isfile(os.path.join(wd, "data.dat")):
        sys.exit(f"No se encuentra data.dat en «{os.path.abspath(wd)}».\n"
                 "Uso: python great_lake.py CARPETA_CON_DATA_DAT\n"
                 "Ejemplo: python great_lake.py ejemplos/caso_simple")
    run(wd)
