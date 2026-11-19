#!/usr/bin/env python3
"""
Great Lake - herramienta grafica para el modelo de recesion de acantilados.

Interfaz de escritorio (PySide6 + matplotlib) sobre ``great_lake_motor.py``
(el motor de la herramienta: ``great_lake.py`` mas opciones adicionales; con
las opciones por defecto da los mismos resultados que el original):

* Carga y guarda ficheros ``data.dat``.
* Edicion del perfil inicial, de los parametros y de los estratos, con
  visualizacion inmediata del perfil, el nivel del lago, el nivel freatico y
  los puntos singulares.
* Ejecucion del modelo en segundo plano, con barra de progreso.
* Visualizacion de los perfiles en cada instante de salida, de las roturas
  (antes/despues) y de graficas de retroceso y volumen.
* Tabla de oleaje por estacion (editable) y carga de series de nivel/mareas.

Uso:
    python great_lake_gui.py [ruta/a/data.dat]

Requiere: numpy, PySide6, matplotlib.
"""

import copy
import csv
import json
import os
import sys
import time

import numpy as np
from PySide6.QtCore import QEvent, QObject, QPointF, QRectF, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QActionGroup, QBrush, QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView, QAbstractScrollArea, QAbstractSpinBox, QApplication, QCheckBox, QComboBox, QDockWidget,
    QDoubleSpinBox, QFileDialog, QFormLayout, QFrame, QGroupBox, QHBoxLayout,
    QHeaderView, QLabel, QMainWindow, QMenu, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QSpinBox, QTableWidget, QTableWidgetItem,
    QTabWidget, QToolButton, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
import matplotlib.cm  # noqa: E402
import matplotlib.colors  # noqa: E402
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas  # noqa: E402
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import great_lake_motor as gl  # noqa: E402
from ayuda_parametros import AYUDA  # noqa: E402

APP_NAME = "Great Lake - Recesión de acantilados"
ACCENT = "#B5651D"
ACCENT_DARK = "#8F4E14"
CANVAS_BG = "#262b33"
CANVAS_FG = "#d9dde3"

ETIQUETAS = {0: "—", 1: "Punto de estudio", 2: "Cresta del talud",
             3: "Punto de rezume", 4: "Base del talud"}
ESTACIONES = [("spring", "Primavera"), ("summer", "Verano"),
              ("fall", "Otoño"), ("winter", "Invierno")]

STYLE = f"""
QMainWindow, QWidget {{ font-size: 10pt; }}
QMainWindow {{ background: #f3f3f3; }}
#topbar {{ background: #fafafa; border-bottom: 2px solid {ACCENT}; }}
#brand {{ color: {ACCENT}; font-weight: bold; font-size: 11pt; border: none;
          padding: 2px 8px; }}
#brand::menu-indicator {{ subcontrol-position: right center; }}
QPushButton {{ background: #ffffff; border: 1px solid #c9c9c9; border-radius: 3px;
               padding: 5px 12px; }}
QPushButton:hover {{ border-color: {ACCENT}; }}
QPushButton:disabled {{ color: #a0a0a0; }}
QPushButton#accent {{ background: {ACCENT}; color: white; border: 1px solid {ACCENT_DARK};
                      font-weight: bold; }}
QPushButton#accent:hover {{ background: {ACCENT_DARK}; }}
QPushButton#info {{ background: #1e88e5; color: white; border-radius: 11px; min-width: 22px;
                    max-width: 22px; min-height: 22px; max-height: 22px; padding: 0;
                    font-weight: bold; }}
#valuelabel {{ font-weight: bold; padding: 0 6px; }}
QTabWidget::pane {{ border: 1px solid #d0d0d0; background: white; }}
QTabBar::tab {{ background: #ececec; border: 1px solid #d0d0d0; padding: 6px 14px;
                color: #555; }}
QTabBar::tab:selected {{ background: white; color: #222; font-weight: bold;
                         border-top: 2px solid {ACCENT}; }}
QHeaderView::section {{ background: #f0f0f0; border: none; border-right: 1px solid #ddd;
                        border-bottom: 1px solid #ddd; padding: 4px; }}
QTableWidget {{ gridline-color: #e3e3e3; selection-background-color: #f6e3d1;
                selection-color: black; }}
QGroupBox {{ font-weight: bold; border: 1px solid #ddd; border-radius: 4px; margin-top: 10px;
             padding-top: 6px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 8px; padding: 0 3px; color: {ACCENT_DARK}; }}
QDockWidget::title {{ background: #ececec; padding: 5px; }}
QProgressBar {{ border: 1px solid #c9c9c9; border-radius: 3px; text-align: center;
                background: white; max-height: 16px; }}
QProgressBar::chunk {{ background: {ACCENT}; }}
QLabel#note {{ color: #666; font-size: 9pt; }}
"""


# ---------------------------------------------------------------------------
# Datos por defecto (caso nuevo)
# ---------------------------------------------------------------------------
def caso_por_defecto():
    """Caso con el que arranca la herramienta: el caso simple realista.

    Acantilado de 10 m a 45 grados (estable al inicio), 20 anos, playa pequena,
    nivel freatico y dos estratos, con parametros dentro de los rangos del
    capitulo 10 de la tesis (igual que ejemplos/caso_simple/data.dat).
    """
    d = dict(rn=1, Utsim=2, tsim=20, UTdib=2, Tdib=1, Dz=0.1, lakelevel=0,
             Qb=5.0, Gs=0.001, groundwater=1, profagua=4.0, ndatosresis=2,
             SFmin=5.0, Kss=2.4e-6, p=1.0, WFmin=70.0, MMM=1e-8,
             Rsc=5.0, Cb=1e-7, KRs=1.0, Bmat=1.0, Mf=2.0)
    for nm, v in (("sigmac", (200.0, 200.0)), ("Shearstrength", (8.0, 8.0)),
                  ("Ecohesion", (15.0, 23.0)), ("Ephi", (25.0, 25.0)),
                  ("pesounit", (20.0, 21.0)), ("potestrato", (4.0, 100.0))):
        d[nm] = np.array([0.0, *v])
    d["Pzy"] = [[10.0000001, -500.0, 0.0], [10.0000001, -10.0, 2.0], [5.0, -5.0, 3.0],
                [2.0, -2.0, 1.0], [0.0, 0.0, 4.0], [-50.0, 1250.0, 0.0]]
    return d


def _fmt(x):
    """Texto de un numero sin perder precision (15 cifras significativas)."""
    return f"{float(x):.15g}"


def _f2(x):
    t = f"{float(x):.2f}"
    return "0.00" if t == "-0.00" else t


def curva_freatica(P, seep, profagua):
    """Nivel freatico (misma formula que el modelo, groundwater = 1)."""
    if seep <= 0:
        return None, None
    prof = abs(P[1, seep] - P[1, 1])
    lon = abs(P[2, seep] - P[2, 1])
    if lon == 0:
        return None, None
    y = P[2, seep] - np.arange(1000) * (lon / 1000.0)
    z = np.sqrt((np.abs(y - P[2, seep]) / lon) * (prof - profagua) ** 2) + P[1, seep]
    return y, z


# ---------------------------------------------------------------------------
# Calculo en segundo plano
# ---------------------------------------------------------------------------
NOMBRES_EROSION = {"trenhaile": "Trenhaile", "castedo": "Castedo"}
NOMBRES_ESTABILIDAD = {"fellenius": "Fellenius", "bishop": "Bishop", "spencer": "Spencer"}


class Worker(QObject):
    """Ejecuta en segundo plano uno o varios calculos (uno por combinacion)."""

    progreso = Signal(int, int, int, int)   # calculo i (desde 1), n calculos, mes, tsim
    terminado = Signal(object, float)
    error = Signal(str)

    def __init__(self, configs):
        super().__init__()
        self.configs = configs
        self._cancelar = False

    def cancelar(self):
        self._cancelar = True

    def run(self):
        t0 = time.time()
        resultados = []
        n = len(self.configs)
        for i, c in enumerate(self.configs, start=1):
            last = [0.0]

            def prog(t, tsim, i=i):
                now = time.time()
                if now - last[0] > 0.05 or t >= tsim:
                    last[0] = now
                    self.progreso.emit(i, n, t, tsim)
                return not self._cancelar

            t1 = time.time()
            try:
                res = gl.run(c["dir"], verbose=False, progress=prog, seasons=c["seasons"],
                             nivel_lago=c["nivel_lago"], estabilidad=c["estabilidad"],
                             mar=c["mar"], erosion=c["erosion"])
                res["nivel_lago"] = c["nivel_lago"]
            except gl.Cancelled:
                self.error.emit("Cálculo detenido por el usuario.")
                return
            except Exception as exc:  # noqa: BLE001
                self.error.emit(f"El cálculo «{c['nombre']}» ha fallado:\n"
                                f"{type(exc).__name__}: {exc}")
                return
            resultados.append({"nombre": c["nombre"], "erosion": c["clave_erosion"],
                               "estabilidad": c["clave_estabilidad"], "dir": c["dir"],
                               "res": res, "segundos": time.time() - t1})
        self.terminado.emit(resultados, time.time() - t0)


# ---------------------------------------------------------------------------
# Widgets auxiliares
# ---------------------------------------------------------------------------
def dspin(minimo, maximo, decimales, paso, sufijo=""):
    w = QDoubleSpinBox()
    w.setRange(minimo, maximo)
    w.setDecimals(decimales)
    w.setSingleStep(paso)
    if sufijo:
        w.setSuffix(" " + sufijo)
    w.setKeyboardTracking(False)
    return w


class SciEdit(QDoubleSpinBox):
    """Campo numerico que admite notacion cientifica (coeficientes pequenos)."""

    def __init__(self):
        super().__init__()
        self.setRange(-1e30, 1e30)
        self.setDecimals(15)
        self.setKeyboardTracking(False)
        self.setButtonSymbols(QDoubleSpinBox.NoButtons)

    def textFromValue(self, v):
        return f"{v:.6g}"

    def valueFromText(self, text):
        try:
            return float(text.replace(",", "."))
        except ValueError:
            return self.value()

    def validate(self, text, pos):
        from PySide6.QtGui import QValidator
        t = text.replace(",", ".")
        if t in ("", "-", "+", ".", "-.", "+.") or t.lower().endswith(("e", "e-", "e+")):
            return QValidator.Intermediate, text, pos
        try:
            float(t)
            return QValidator.Acceptable, text, pos
        except ValueError:
            return QValidator.Invalid, text, pos


class Lienzo(FigureCanvas):
    def __init__(self, oscuro=False):
        self.fig = Figure(figsize=(6, 4), dpi=100)
        super().__init__(self.fig)
        self.oscuro = oscuro
        self.ax = self.fig.add_subplot(111)
        self.estilo()

    def estilo(self):
        ax = self.ax
        if self.oscuro:
            self.fig.set_facecolor(CANVAS_BG)
            ax.set_facecolor(CANVAS_BG)
            for s in ax.spines.values():
                s.set_color("#5a6270")
            ax.tick_params(colors=CANVAS_FG, labelsize=8)
            ax.xaxis.label.set_color(CANVAS_FG)
            ax.yaxis.label.set_color(CANVAS_FG)
            ax.title.set_color(CANVAS_FG)
            ax.grid(True, color="#3a414d", lw=0.6)
        else:
            self.fig.set_facecolor("white")
            ax.tick_params(labelsize=8)
            ax.grid(True, color="#e6e6e6", lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)


class FiltroRueda(QObject):
    """Evita que la rueda del raton cambie desplegables y campos numericos.

    Al girar la rueda sobre un QComboBox o un campo numerico no se cambia su
    valor: el giro se pasa al panel desplazable que lo contiene (que se
    desplaza, como en cualquier otra parte). Para cambiar el valor hay que
    hacer clic, escribir o usar el teclado.
    """

    def eventFilter(self, obj, ev):
        if ev.type() == QEvent.Wheel and isinstance(obj, (QComboBox, QAbstractSpinBox)):
            padre = obj.parentWidget()
            while padre is not None and not isinstance(padre, QAbstractScrollArea):
                padre = padre.parentWidget()
            if padre is not None:
                QApplication.sendEvent(padre.viewport(), ev)
            return True
        return False


class InfoBoton(QToolButton):
    """Boton (i) que despliega una explicacion sencilla y, debajo, la detallada."""

    def __init__(self, clave, parent=None):
        super().__init__(parent)
        self.clave = clave
        self.setText("i")
        self.setCursor(Qt.PointingHandCursor)
        self.setToolTip("Qué es esto")
        self.setFixedSize(18, 18)
        self.setStyleSheet(
            "QToolButton { background: #e7eef8; color: #1e66c8; border: 1px solid #9dbbe6;"
            " border-radius: 9px; font-weight: bold; font-family: serif; font-size: 9pt;"
            " padding: 0; }"
            "QToolButton:hover { background: #1e88e5; color: white; border-color: #1e88e5; }")
        self.clicked.connect(self.mostrar)
        self._popup = None

    def mostrar(self):
        sencillo, detalle = AYUDA.get(self.clave, ("(Sin ayuda)", None))
        pop = QFrame(self, Qt.Popup)
        pop.setObjectName("infopopup")
        pop.setStyleSheet(
            "#infopopup { background: white; border: 1px solid #9dbbe6; border-radius: 6px; }"
            "QLabel { background: transparent; }")
        lay = QVBoxLayout(pop)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(6)
        l1 = QLabel(sencillo)
        l1.setWordWrap(True)
        l1.setStyleSheet("font-size: 10pt; font-weight: bold; color: #1d2a3a;")
        lay.addWidget(l1)
        if detalle:
            sep = QFrame()
            sep.setFrameShape(QFrame.HLine)
            sep.setStyleSheet("color: #dfe6ef;")
            lay.addWidget(sep)
            titulo = QLabel("En detalle")
            titulo.setStyleSheet(f"font-size: 8pt; font-weight: bold; color: {ACCENT_DARK};"
                                 " text-transform: uppercase;")
            lay.addWidget(titulo)
            l2 = QLabel(detalle)
            l2.setWordWrap(True)
            l2.setTextFormat(Qt.RichText)
            l2.setStyleSheet("font-size: 9pt; color: #3a4554;")
            lay.addWidget(l2)
        pop.setFixedWidth(420)
        pop.adjustSize()
        # debajo del boton, sin salirse de la pantalla
        pos = self.mapToGlobal(self.rect().bottomLeft())
        pantalla = self.screen().availableGeometry() if self.screen() else None
        x, y = pos.x() - 10, pos.y() + 4
        if pantalla is not None:
            x = max(pantalla.left() + 4, min(x, pantalla.right() - pop.width() - 4))
            if y + pop.height() > pantalla.bottom():
                y = self.mapToGlobal(self.rect().topLeft()).y() - pop.height() - 4
        pop.move(x, y)
        pop.show()
        self._popup = pop


def etiqueta_info(texto, clave):
    """Etiqueta de un formulario con su boton (i) a la derecha."""
    w = QWidget()
    h = QHBoxLayout(w)
    h.setContentsMargins(0, 0, 0, 0)
    h.setSpacing(5)
    h.addWidget(QLabel(texto))
    h.addWidget(InfoBoton(clave))
    h.addStretch(1)
    return w


COLOR_MES_CALCULADO = "#3a414d"
COLOR_INVIERNO = "#9cc9f0"      # azul hielo
COLOR_DERRUBIOS = "#d9a25b"     # arena
COLOR_ROTURA = "#ff4d4d"
COLOR_VUELCO = "#ffb020"
COLOR_CURSOR = "#00e5ff"


class LineaTemporal(QWidget):
    """Barra de tiempo mes a mes con marcas de roturas y meses sin calculo."""

    valor_cambiado = Signal(int)

    def __init__(self):
        super().__init__()
        self.tsim = 0
        self.valor = 0
        self.roturas = {}     # mes -> numero de deslizamientos circulares
        self.vuelcos = {}     # mes -> numero de vuelcos
        self.sin_calc = {}    # mes -> "invierno" | "derrubios"
        self.texto_invierno = "invierno: lago helado, sin oleaje ni erosión"
        self.setMinimumHeight(58)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setMouseTracking(True)

    def configurar(self, tsim, roturas, invierno, derrubios, vuelcos=()):
        self.tsim = int(tsim)
        self.roturas = {}
        for m in roturas:
            self.roturas[m] = self.roturas.get(m, 0) + 1
        self.vuelcos = {}
        for m in vuelcos:
            self.vuelcos[m] = self.vuelcos.get(m, 0) + 1
        self.sin_calc = {m: "invierno" for m in invierno}
        self.sin_calc.update({m: "derrubios" for m in derrubios})
        self.valor = 0
        self.update()

    def _x0x1(self):
        return 14.0, float(self.width() - 14)

    def _x(self, mes):
        x0, x1 = self._x0x1()
        return x0 + (x1 - x0) * (mes / max(1, self.tsim))

    def set_valor(self, v, emitir=True):
        v = max(0, min(self.tsim, int(v)))
        if v != self.valor:
            self.valor = v
            self.update()
            if emitir:
                self.valor_cambiado.emit(v)

    def _valor_en(self, x):
        x0, x1 = self._x0x1()
        return round((x - x0) / max(1.0, x1 - x0) * self.tsim)

    def mousePressEvent(self, ev):
        if self.tsim:
            self.set_valor(self._valor_en(ev.position().x()))

    def describir_mes(self, mes):
        if mes <= 0:
            return "Mes 0: perfil inicial"
        tipo = self.sin_calc.get(mes)
        if tipo == "invierno":
            estado = self.texto_invierno
        elif tipo == "derrubios":
            estado = "protegido por derrubios: sin erosión"
        else:
            estado = "calculado (oleaje y erosión)"
        n = self.roturas.get(mes, 0)
        nv = self.vuelcos.get(mes, 0)
        txt = f"Mes {mes} (año {mes / 12:.2f}): {estado}"
        if n:
            txt += f" · {n} deslizamiento{'s' if n > 1 else ''}"
        if nv:
            txt += f" · {nv} vuelco{'s' if nv > 1 else ''}"
        return txt

    def mouseMoveEvent(self, ev):
        if not self.tsim:
            return
        if ev.buttons() & Qt.LeftButton:
            self.set_valor(self._valor_en(ev.position().x()))
        from PySide6.QtWidgets import QToolTip
        mes = max(0, min(self.tsim, self._valor_en(ev.position().x())))
        QToolTip.showText(ev.globalPosition().toPoint(), self.describir_mes(mes), self)

    def wheelEvent(self, ev):
        if self.tsim:
            self.set_valor(self.valor + (1 if ev.angleDelta().y() < 0 else -1))

    def keyPressEvent(self, ev):
        pasos = {Qt.Key_Right: 1, Qt.Key_Left: -1, Qt.Key_PageDown: 12, Qt.Key_PageUp: -12}
        if ev.key() in pasos:
            self.set_valor(self.valor + pasos[ev.key()])
        elif ev.key() == Qt.Key_Home:
            self.set_valor(0)
        elif ev.key() == Qt.Key_End:
            self.set_valor(self.tsim)
        else:
            super().keyPressEvent(ev)

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.fillRect(self.rect(), QColor(CANVAS_BG))
        x0, x1 = self._x0x1()
        yb, hb = 18.0, 12.0
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(COLOR_MES_CALCULADO))
        p.drawRoundedRect(QRectF(x0, yb, x1 - x0, hb), 4, 4)
        if not self.tsim:
            p.setPen(QColor("#8a93a3"))
            p.drawText(self.rect(), Qt.AlignCenter, "Calcula para activar la línea temporal")
            return
        w = (x1 - x0) / self.tsim
        # meses sin calculo
        for m, tipo in self.sin_calc.items():
            c = QColor(COLOR_INVIERNO if tipo == "invierno" else COLOR_DERRUBIOS)
            p.fillRect(QRectF(self._x(m - 1), yb, max(1.0, w), hb), c)
        # parte recorrida: linea fina bajo la barra (no tapa los colores)
        p.fillRect(QRectF(x0, yb + hb + 1.5, max(0.0, self._x(self.valor) - x0), 2.5),
                   QColor(COLOR_CURSOR))
        # vuelcos (rayas cortas ambar) y deslizamientos (rayas rojas)
        for m, n in self.vuelcos.items():
            x = self._x(m - 0.5)
            p.setPen(QPen(QColor(COLOR_VUELCO), 1.6))
            p.drawLine(QPointF(x, yb - 4), QPointF(x, yb - 1))
        for m, n in self.roturas.items():
            x = self._x(m - 0.5)
            p.setPen(QPen(QColor(COLOR_ROTURA), 2 if n > 1 else 1.4))
            p.drawLine(QPointF(x, yb - 7 - min(6, 2 * (n - 1))), QPointF(x, yb - 1))
        # escala en anos
        anos = self.tsim / 12.0
        paso = 1
        for cand in (1, 2, 5, 10, 20, 25, 50, 100):
            if (x1 - x0) / max(1e-9, anos / cand) >= 42:
                paso = cand
                break
        f = p.font()
        f.setPointSizeF(7.5)
        p.setFont(f)
        a = 0
        while a * 12 <= self.tsim:
            x = self._x(a * 12)
            p.setPen(QPen(QColor("#8a93a3"), 1))
            p.drawLine(QPointF(x, yb + hb + 4), QPointF(x, yb + hb + 7))
            p.drawText(QRectF(x - 20, yb + hb + 7, 40, 14), Qt.AlignHCenter, f"{a}")
            a += paso
        ultimo = self._x((a - paso) * 12)
        if x1 - ultimo > 60:
            p.drawText(QRectF(x1 - 60, yb + hb + 7, 60, 14), Qt.AlignRight, "años")
        # cursor
        x = self._x(self.valor)
        p.setPen(QPen(QColor("#00e5ff"), 2))
        p.drawLine(QPointF(x, yb - 6), QPointF(x, yb + hb + 3))
        p.setBrush(QBrush(QColor("#00e5ff")))
        p.setPen(QPen(QColor("black"), 1))
        p.drawEllipse(QPointF(x, yb + hb / 2), 6, 6)
        if self.hasFocus():
            p.setPen(QPen(QColor("#00e5ff"), 1, Qt.DotLine))
            p.setBrush(Qt.NoBrush)
            p.drawRect(self.rect().adjusted(0, 0, -1, -1))


# ---------------------------------------------------------------------------
# Ventana principal
# ---------------------------------------------------------------------------
class Ventana(QMainWindow):
    @property
    def seasons(self):
        return self.seasons_por_modo[self.modo]

    @seasons.setter
    def seasons(self, valor):
        self.seasons_por_modo[self.modo] = valor

    def _seasons_defecto(self):
        return gl.SEASONS if self.modo == "lago" else gl.SEASONS_MAR

    def __init__(self, ruta=None):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1600, 950)
        self._filtro_rueda = FiltroRueda(self)
        QApplication.instance().installEventFilter(self._filtro_rueda)
        self.ruta = None
        self.res = None
        self.res_dir = None
        self.serie_nivel = None
        self.modo = "lago"
        self._eventos = []
        self.resultados = []
        self.seasons_por_modo = {"lago": copy.deepcopy(gl.SEASONS),
                                 "mar": copy.deepcopy(gl.SEASONS_MAR)}
        self._bloqueo = False
        self.hilo = None
        self.worker = None

        self._crear_barra()
        self._crear_panel_datos()
        self._crear_centro()
        self._crear_panel_derecho()
        self.statusBar().showMessage("Listo")

        if ruta:
            self.abrir(ruta)
        else:
            self.cargar_datos(caso_por_defecto())
            self.dibujar_grafica()

    # ------------------------------------------------------------------ UI
    def _crear_barra(self):
        barra = QWidget(objectName="topbar")
        lay = QHBoxLayout(barra)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(8)

        marca = QToolButton(objectName="brand")
        marca.setText("LAGO")
        marca.setPopupMode(QToolButton.InstantPopup)
        self.marca = marca
        menu = QMenu(marca)
        grupo_modo = QActionGroup(self)
        self.acciones_modo = {}
        for clave, texto in (("lago", "Lago: sin marea (modelo Great_lake de la tesis)"),
                             ("mar", "Mar: con marea (adaptación de la herramienta)")):
            a = QAction(texto, self, checkable=True)
            a.setChecked(clave == "lago")
            a.triggered.connect(lambda _=False, c=clave: self.cambiar_modo(c))
            grupo_modo.addAction(a)
            menu.addAction(a)
            self.acciones_modo[clave] = a
        menu.addSeparator()
        for texto, fn in (("Nuevo caso", self.nuevo), ("Abrir data.dat…", self.dialogo_abrir),
                          ("Guardar data.dat", self.guardar),
                          ("Guardar como…", self.guardar_como), (None, None),
                          ("Abrir carpeta de resultados…", self.cargar_resultados_dialogo),
                          (None, None), ("Salir", self.close)):
            if texto is None:
                menu.addSeparator()
            else:
                a = QAction(texto, self)
                a.triggered.connect(fn)
                menu.addAction(a)
        marca.setMenu(menu)
        lay.addWidget(marca)
        lay.addWidget(InfoBoton("modo"))

        b_abrir = QPushButton("Abrir data.dat…")
        b_abrir.clicked.connect(self.dialogo_abrir)
        lay.addWidget(b_abrir)
        self.b_calc = QPushButton("Calcular", objectName="accent")
        self.b_calc.clicked.connect(self.calcular)
        lay.addWidget(self.b_calc)
        b_guardar = QPushButton("Guardar")
        b_guardar.clicked.connect(self.guardar)
        lay.addWidget(b_guardar)
        lay.addWidget(self._sep())

        lay.addWidget(QLabel("Duración:"))
        self.sp_tsim = QSpinBox()
        self.sp_tsim.setRange(1, 100000)
        self.cb_utsim = QComboBox()
        self.cb_utsim.addItems(["meses", "años"])
        lay.addWidget(self.sp_tsim)
        lay.addWidget(self.cb_utsim)
        lay.addWidget(QLabel("Salida cada:"))
        self.sp_tdib = QSpinBox()
        self.sp_tdib.setRange(1, 100000)
        self.cb_utdib = QComboBox()
        self.cb_utdib.addItems(["meses", "años"])
        lay.addWidget(self.sp_tdib)
        lay.addWidget(self.cb_utdib)
        lay.addWidget(self._sep())

        self.lb_nivel_barra = QLabel("Nivel del lago:")
        lay.addWidget(self.lb_nivel_barra)
        self.sp_lake = QSpinBox()
        self.sp_lake.setRange(-100, 100)
        self.sp_lake.setSuffix(" m")
        lay.addWidget(self.sp_lake)
        lay.addWidget(QLabel("Δz:"))
        self.sp_dz = dspin(0.01, 5, 3, 0.05, "m")
        lay.addWidget(self.sp_dz)
        self.lb_npts = QLabel("", objectName="valuelabel")
        lay.addWidget(self.lb_npts)
        lay.addWidget(self._sep())

        lay.addWidget(QLabel("Vista:"))
        self.cb_vista = QComboBox()
        self.cb_vista.addItems(["Acantilado", "Perfil completo"])
        self.cb_vista.currentIndexChanged.connect(self.dibujar_perfil)
        lay.addWidget(self.cb_vista)
        b_info = QPushButton("i", objectName="info")
        b_info.clicked.connect(self.acerca_de)
        lay.addWidget(b_info)
        self.ck_leyenda = QCheckBox("Leyenda")
        self.ck_leyenda.setChecked(True)
        self.ck_leyenda.toggled.connect(self.dibujar_perfil)
        lay.addWidget(self.ck_leyenda)
        lay.addStretch(1)
        self.progreso = QProgressBar()
        self.progreso.setFixedWidth(220)
        self.progreso.setVisible(False)
        lay.addWidget(self.progreso)

        for w in (self.sp_tsim, self.sp_tdib, self.sp_lake):
            w.valueChanged.connect(self._cambio_datos)
        self.sp_dz.valueChanged.connect(self._cambio_datos)
        for w in (self.cb_utsim, self.cb_utdib):
            w.currentIndexChanged.connect(self._cambio_datos)

        self.setMenuWidget(barra)

    @staticmethod
    def _sep():
        f = QFrame()
        f.setFrameShape(QFrame.VLine)
        f.setStyleSheet("color: #d0d0d0;")
        return f

    def _crear_panel_datos(self):
        dock = QDockWidget("Datos", self)
        dock.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable)
        tabs = QTabWidget()

        # --- Perfil
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        b_add = QPushButton("Añadir punto")
        b_add.clicked.connect(self.anadir_punto)
        b_del = QPushButton("Eliminar")
        b_del.clicked.connect(self.eliminar_punto)
        b_ord = QPushButton("Ordenar por cota")
        b_ord.clicked.connect(self.ordenar_perfil)
        fila.addWidget(b_add)
        fila.addWidget(b_del)
        fila.addWidget(b_ord)
        fila.addStretch(1)
        fila.addWidget(etiqueta_info("Ayuda", "perfil"))
        v.addLayout(fila)
        self.tb_perfil = QTableWidget(0, 3)
        self.tb_perfil.setHorizontalHeaderLabels(["Cota z (m)", "Distancia y (m)", "Punto"])
        self.tb_perfil.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_perfil.verticalHeader().setDefaultSectionSize(28)
        self.tb_perfil.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tb_perfil.itemChanged.connect(self._cambio_datos)
        v.addWidget(self.tb_perfil)
        nota = QLabel("Puntos de arriba (tierra) a abajo (lago), con la cota decreciente. "
                      "Hacen falta la cresta del talud y la base del talud; el punto de "
                      "rezume, si hay nivel freático.", objectName="note")
        nota.setWordWrap(True)
        v.addWidget(nota)
        tabs.addTab(w, "Perfil")

        # --- Parametros
        cont = QWidget()
        v = QVBoxLayout(cont)
        self.campos = {}

        def grupo(titulo, filas):
            g = QGroupBox(titulo)
            f = QFormLayout(g)
            for clave, etiqueta, widget in filas:
                self.campos[clave] = widget
                f.addRow(etiqueta_info(etiqueta, clave) if clave in AYUDA else etiqueta, widget)
                if isinstance(widget, QComboBox):
                    widget.currentIndexChanged.connect(self._cambio_datos)
                else:
                    widget.valueChanged.connect(self._cambio_datos)
            v.addWidget(g)

        sp_rn = QSpinBox()
        sp_rn.setRange(0, 10 ** 6)
        grupo("Simulación", [("rn", "Número de ejecución", sp_rn)])
        grupo("Erosión por esfuerzo cortante",
              [("SFmin", "Umbral de erosión [Pa]", SciEdit()),
               ("Kss", "Coeficiente Kss", SciEdit()),
               ("p", "Coeficiente de calibración p", SciEdit())])
        grupo("Erosión por impacto del oleaje",
              [("WFmin", "Fuerza umbral [kg/m²]", SciEdit()),
               ("MMM", "Coeficiente MMM", SciEdit())])
        grupo("Run-up",
              [("Rsc", "Distancia umbral Rsc", SciEdit()),
               ("Cb", "Coeficiente Cb", SciEdit()),
               ("KRs", "Coeficiente KRs", SciEdit()),
               ("Bmat", "Factor material playa", SciEdit())])
        grupo("Playa",
              [("Qb", "Volumen de arena [m³/m]", SciEdit()),
               ("Gs", "Diámetro del grano [m]", SciEdit())])
        cb_gw = QComboBox()
        cb_gw.addItems(["0: sin nivel freático", "1: nivel freático", "2: saturado"])
        grupo("Agua subterránea",
              [("groundwater", "Tipo", cb_gw),
               ("profagua", "Profundidad tierra adentro [m]", SciEdit())])
        grupo("Derrubios", [("Mf", "Coef. de movilidad [m²/día]", SciEdit())])

        g = QGroupBox("Métodos de cálculo (marca uno o varios para compararlos)")
        f = QFormLayout(g)
        self.ck_ero = {}
        fila_e = QWidget()
        he = QHBoxLayout(fila_e)
        he.setContentsMargins(0, 0, 0, 0)
        for clave, texto, ayuda in (("trenhaile", "Trenhaile (original del lago)", "erosion_trenhaile"),
                                    ("castedo", "Castedo (marina, tesis)", "erosion_castedo")):
            ck = QCheckBox(texto)
            ck.toggled.connect(self._cambio_estabilidad)
            self.ck_ero[clave] = ck
            he.addWidget(ck)
            he.addWidget(InfoBoton(ayuda))
            he.addSpacing(8)
        he.addStretch(1)
        self.ck_ero["trenhaile"].setChecked(True)
        f.addRow(etiqueta_info("Erosión", "erosion"), fila_e)
        self.ck_est = {}
        fila_s = QWidget()
        hs = QHBoxLayout(fila_s)
        hs.setContentsMargins(0, 0, 0, 0)
        for clave, texto, ayuda in (("fellenius", "Fellenius (original)", "estab_fellenius"),
                                    ("bishop", "Bishop", "estab_bishop"),
                                    ("spencer", "Spencer", "estab_spencer")):
            ck = QCheckBox(texto)
            ck.toggled.connect(self._cambio_estabilidad)
            self.ck_est[clave] = ck
            hs.addWidget(ck)
            hs.addWidget(InfoBoton(ayuda))
            hs.addSpacing(8)
        hs.addStretch(1)
        self.ck_est["fellenius"].setChecked(True)
        f.addRow(etiqueta_info("Estabilidad", "metodo"), fila_s)
        self.sp_rebanadas = QSpinBox()
        self.sp_rebanadas.setRange(1, 200)
        self.sp_rebanadas.setValue(5)
        self.cb_criterio = QComboBox()
        for clave, nombre in (("todas", "Todas las inestables (original)"),
                              ("profunda", "La más profunda (tesis)"),
                              ("minimo", "La de menor FS")):
            self.cb_criterio.addItem(nombre, clave)
        f.addRow(etiqueta_info("Nº de rebanadas", "rebanadas"), self.sp_rebanadas)
        f.addRow(etiqueta_info("Criterio de rotura", "criterio"), self.cb_criterio)
        b = QPushButton("Valores del original")
        b.clicked.connect(self._estabilidad_original)
        f.addRow("", b)
        self.lb_aviso_estab = QLabel("", objectName="note")
        self.lb_aviso_estab.setWordWrap(True)
        f.addRow(self.lb_aviso_estab)
        self.cb_criterio.currentIndexChanged.connect(self._cambio_estabilidad)
        self.sp_rebanadas.valueChanged.connect(self._cambio_estabilidad)
        v.addWidget(g)

        g = QGroupBox("Erosión de Castedo (si está marcada)")
        f = QFormLayout(g)
        self.sp_K = dspin(0.01, 1e7, 2, 10, "")
        self.sp_K.setValue(gl.CASTEDO_DEFECTO["K"])
        self.sp_pmax = dspin(5, 89, 1, 5, "°")
        self.sp_pmax.setValue(gl.CASTEDO_DEFECTO["pendiente_max"])
        self.sp_fdesl = dspin(0, 100, 0, 5, "%")
        self.sp_fdesl.setValue(100 * gl.CASTEDO_DEFECTO["derrubios_deslizamiento"])
        self.sp_fvuel = dspin(0, 100, 0, 5, "%")
        self.sp_fvuel.setValue(100 * gl.CASTEDO_DEFECTO["derrubios_vuelco"])
        f.addRow(etiqueta_info("Constante de calibración K", "castedo_K"), self.sp_K)
        f.addRow(etiqueta_info("Pendiente máxima", "castedo_pendiente"), self.sp_pmax)
        f.addRow(etiqueta_info("Derrubios de un deslizamiento", "castedo_derrubios"), self.sp_fdesl)
        f.addRow(etiqueta_info("Derrubios de un vuelco", "castedo_derrubios"), self.sp_fvuel)
        b = QPushButton("Valores de la tesis (Holderness)")
        b.clicked.connect(self._castedo_tesis)
        f.addRow("", b)
        v.addWidget(g)
        v.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(cont)
        tabs.addTab(scroll, "Parámetros")

        # --- Estratos
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        b_add = QPushButton("Añadir estrato")
        b_add.clicked.connect(self.anadir_estrato)
        b_del = QPushButton("Eliminar")
        b_del.clicked.connect(self.eliminar_estrato)
        fila.addWidget(b_add)
        fila.addWidget(b_del)
        fila.addStretch(1)
        fila.addWidget(etiqueta_info("Ayuda", "estratos"))
        v.addLayout(fila)
        self.tb_estratos = QTableWidget(0, 6)
        self.tb_estratos.setHorizontalHeaderLabels(
            ["UCS\n(kPa)", "Resist. corte\n(kPa)", "Cohesión\n(kPa)", "Fricción\n(°)",
             "Peso esp.\n(kN/m³)", "Espesor\n(m)"])
        self.tb_estratos.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_estratos.itemChanged.connect(self._cambio_datos)
        v.addWidget(self.tb_estratos)
        nota = QLabel("De arriba abajo, empezando en la cota del primer punto del perfil.",
                      objectName="note")
        nota.setWordWrap(True)
        v.addWidget(nota)
        tabs.addTab(w, "Estratos")

        # --- Resultados (perfiles por instante de salida)
        w = QWidget()
        v = QVBoxLayout(w)
        self.fila_resultado = QWidget()
        fr = QHBoxLayout(self.fila_resultado)
        fr.setContentsMargins(0, 0, 0, 0)
        fr.addWidget(QLabel("Cálculo mostrado:"))
        self.cb_resultado = QComboBox()
        self.cb_resultado.currentIndexChanged.connect(self._elegir_resultado)
        fr.addWidget(self.cb_resultado, 1)
        self.fila_resultado.setVisible(False)
        v.addWidget(self.fila_resultado)
        fila = QHBoxLayout()
        b_todos = QPushButton("Ver todos")
        b_todos.clicked.connect(lambda: self._marcar_salidas(True))
        b_ninguno = QPushButton("Ninguno")
        b_ninguno.clicked.connect(lambda: self._marcar_salidas(False))
        fila.addWidget(b_todos)
        fila.addWidget(b_ninguno)
        fila.addStretch(1)
        v.addLayout(fila)
        self.tb_salidas = QTableWidget(0, 5)
        self.tb_salidas.setHorizontalHeaderLabels(
            ["Ver", "Mes", "Año", "Retroceso\nnivel lago (m)", "Retroceso\ncresta (m)"])
        h = self.tb_salidas.horizontalHeader()
        h.setSectionResizeMode(QHeaderView.Stretch)
        h.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.tb_salidas.verticalHeader().setVisible(False)
        self.tb_salidas.itemChanged.connect(self._cambio_salidas)
        v.addWidget(self.tb_salidas)
        self.ck_primer = QCheckBox("Mostrar el perfil tras el primer mes calculado "
                                   "(roturas iniciales)")
        self.ck_primer.setChecked(True)
        self.ck_primer.toggled.connect(self.dibujar_perfil)
        v.addWidget(self.ck_primer)
        nota = QLabel("La línea blanca es el perfil inicial de los datos (año 0). Si el "
                      "perfil inicial es inestable con los parámetros del terreno, en el "
                      "primer mes calculado (mes 4; enero-marzo no se calculan) se producen "
                      "varias roturas seguidas: véase la pestaña «Roturas».",
                      objectName="note")
        nota.setWordWrap(True)
        v.addWidget(nota)
        self.lb_resdir = QLabel("Sin resultados. Pulsa «Calcular».", objectName="note")
        self.lb_resdir.setWordWrap(True)
        v.addWidget(self.lb_resdir)
        tabs.addTab(w, "Resultados")
        self.tabs_datos = tabs

        dock.setWidget(tabs)
        dock.setMinimumWidth(470)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def _crear_centro(self):
        w = QWidget()
        w.setStyleSheet(f"background: {CANVAS_BG};")
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        self.lienzo = Lienzo(oscuro=True)
        self._bg = None
        self._linea_mes = None
        self._modo_activo = False
        self.lienzo.mpl_connect("draw_event", self._on_draw)
        nav = NavigationToolbar2QT(self.lienzo, w)
        nav.setStyleSheet("background: #e9e9e9; border: none;")
        v.addWidget(nav)
        v.addWidget(self.lienzo, 1)

        barra = QWidget()
        barra.setStyleSheet(f"background: {CANVAS_BG}; color: {CANVAS_FG};")
        bl = QVBoxLayout(barra)
        bl.setContentsMargins(8, 2, 8, 6)
        bl.setSpacing(2)
        fila = QHBoxLayout()
        estilo_b = ("QPushButton { background: #3a414d; color: #e6e9ee; border: 1px solid #4b5462;"
                    " border-radius: 3px; padding: 3px 8px; min-width: 26px; }"
                    "QPushButton:hover { border-color: %s; }"
                    "QPushButton:disabled { color: #6b7380; }" % ACCENT)
        self.b_inicio = QPushButton("⏮")
        self.b_atras = QPushButton("▕◀")
        self.b_play = QPushButton("▶")
        self.b_adelante = QPushButton("▶▏")
        self.b_fin = QPushButton("⏭")
        for b, tip in ((self.b_inicio, "Ir al inicio (perfil inicial)"),
                       (self.b_atras, "Mes anterior (←)"), (self.b_play, "Reproducir / pausa"),
                       (self.b_adelante, "Mes siguiente (→)"), (self.b_fin, "Ir al final")):
            b.setToolTip(tip)
        for b, fn in ((self.b_inicio, lambda: self.linea.set_valor(0)),
                      (self.b_atras, lambda: self.linea.set_valor(self.linea.valor - 1)),
                      (self.b_play, self.reproducir),
                      (self.b_adelante, lambda: self.linea.set_valor(self.linea.valor + 1)),
                      (self.b_fin, lambda: self.linea.set_valor(self.linea.tsim))):
            b.setStyleSheet(estilo_b)
            b.clicked.connect(fn)
            fila.addWidget(b)
        fila.addSpacing(8)
        fila.addWidget(QLabel("Velocidad:"))
        self.cb_vel = QComboBox()
        self.cb_vel.addItems(["Lenta", "Normal", "Rápida", "Muy rápida"])
        self.cb_vel.setCurrentIndex(1)
        self.cb_vel.setStyleSheet("background: #3a414d; color: #e6e9ee;")
        self.cb_vel.currentIndexChanged.connect(self._velocidad)
        fila.addWidget(self.cb_vel)
        fila.addStretch(1)
        leyenda = QLabel(self._texto_leyenda())
        self.lb_leyenda_linea = leyenda
        leyenda.setStyleSheet("color: #c3cad5; font-size: 8pt;")
        fila.addWidget(leyenda)
        fila.addWidget(InfoBoton("linea_temporal"))
        bl.addLayout(fila)
        self.linea = LineaTemporal()
        self.linea.valor_cambiado.connect(self._cambio_mes)
        bl.addWidget(self.linea)
        self.lb_mes = QLabel("")
        self.lb_mes.setStyleSheet(f"color: {CANVAS_FG}; font-weight: bold;")
        bl.addWidget(self.lb_mes)
        v.addWidget(barra)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._tick)
        self._velocidad()
        self._activar_linea(False)
        self.setCentralWidget(w)

    def _crear_panel_derecho(self):
        dock = QDockWidget("", self)
        dock.setTitleBarWidget(QWidget())
        dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        tabs = QTabWidget()

        # --- Recesion
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        fila.addWidget(QLabel("Gráfica:"))
        self.cb_grafica = QComboBox()
        self.cb_grafica.addItems(["Retroceso en el nivel del lago", "Retroceso de la cresta",
                                  "Retroceso en los puntos de estudio",
                                  "Volumen de roturas acumulado (deslizamientos + vuelcos)"])
        self.cb_grafica.currentIndexChanged.connect(self.dibujar_grafica)
        fila.addWidget(self.cb_grafica, 1)
        v.addLayout(fila)
        self.graf = Lienzo()
        v.addWidget(self.graf, 1)
        fila = QHBoxLayout()
        fila.addStretch(1)
        b = QPushButton("Guardar figura…")
        b.clicked.connect(lambda: self.guardar_figura(self.graf.fig))
        fila.addWidget(b)
        b = QPushButton("Exportar CSV…")
        b.clicked.connect(self.exportar_csv)
        fila.addWidget(b)
        v.addLayout(fila)
        tabs.addTab(w, "Recesión")

        # --- Comparacion de metodos
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        fila.addWidget(etiqueta_info("Comparación de los cálculos", "comparacion"))
        fila.addStretch(1)
        v.addLayout(fila)
        self.tb_comp = QTableWidget(0, 7)
        self.tb_comp.setHorizontalHeaderLabels(
            ["Erosión", "Estabilidad", "Desliz.", "Vuelcos", "Cresta\n(m/año)",
             "Nivel agua\n(m/año)", "Tiempo\n(s)"])
        self.tb_comp.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_comp.verticalHeader().setVisible(False)
        self.tb_comp.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tb_comp.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_comp.setMaximumHeight(190)
        self.tb_comp.itemSelectionChanged.connect(self._fila_comparacion)
        v.addWidget(self.tb_comp)
        fila = QHBoxLayout()
        fila.addWidget(QLabel("Gráfica:"))
        self.cb_comp = QComboBox()
        self.cb_comp.addItems(["Retroceso de la cresta", "Retroceso en el nivel del agua",
                               "Volumen de roturas acumulado", "Perfil final"])
        self.cb_comp.currentIndexChanged.connect(self.dibujar_comparacion)
        fila.addWidget(self.cb_comp, 1)
        v.addLayout(fila)
        self.graf_comp = Lienzo()
        v.addWidget(self.graf_comp, 1)
        fila = QHBoxLayout()
        fila.addStretch(1)
        b = QPushButton("Guardar figura…")
        b.clicked.connect(lambda: self.guardar_figura(self.graf_comp.fig))
        fila.addWidget(b)
        b = QPushButton("Exportar tabla CSV…")
        b.clicked.connect(self.exportar_comparacion)
        fila.addWidget(b)
        v.addLayout(fila)
        tabs.addTab(w, "Comparación")

        # --- Roturas
        w = QWidget()
        v = QVBoxLayout(w)
        self.tb_roturas = QTableWidget(0, 5)
        self.tb_roturas.setHorizontalHeaderLabels(["Mes", "Año", "Tipo", "Volumen (m³/m)", "FS"])
        self.tb_roturas.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_roturas.verticalHeader().setVisible(False)
        self.tb_roturas.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tb_roturas.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_roturas.itemSelectionChanged.connect(self.dibujar_perfil)
        v.addWidget(self.tb_roturas)
        nota = QLabel("Selecciona una o varias roturas para ver el perfil antes (discontinuo) "
                      "y después en la vista central.", objectName="note")
        nota.setWordWrap(True)
        fila_n = QHBoxLayout()
        fila_n.addWidget(nota, 1)
        fila_n.addWidget(etiqueta_info("Tipos de rotura", "vuelco"))
        v.addLayout(fila_n)
        b = QPushButton("Quitar selección")
        b.clicked.connect(self.tb_roturas.clearSelection)
        v.addWidget(b)
        tabs.addTab(w, "Roturas")

        # --- Oleaje
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        fila.addWidget(etiqueta_info("Estación:", "oleaje"))
        self.cb_estacion = QComboBox()
        for _, nombre in ESTACIONES:
            self.cb_estacion.addItem(nombre)
        self.cb_estacion.currentIndexChanged.connect(self.mostrar_oleaje)
        fila.addWidget(self.cb_estacion, 1)
        v.addLayout(fila)
        fila = QHBoxLayout()
        b = QPushButton("Añadir ola")
        b.clicked.connect(self.anadir_ola)
        fila.addWidget(b)
        b = QPushButton("Eliminar")
        b.clicked.connect(self.eliminar_ola)
        fila.addWidget(b)
        b = QPushButton("Restablecer originales")
        b.clicked.connect(self.restablecer_oleaje)
        fila.addWidget(b)
        self.b_temporales = QPushButton("Con temporales (ilustrativo)")
        self.b_temporales.setToolTip("Distribución de olas con temporales cuya altura media "
                                     "coincide con la de Holderness; útil para la erosión de "
                                     "Trenhaile, que solo actúa por encima de unos umbrales.")
        self.b_temporales.clicked.connect(self._oleaje_temporales)
        self.b_temporales.setVisible(False)
        fila.addWidget(self.b_temporales)
        v.addLayout(fila)
        self.tb_oleaje = QTableWidget(0, 3)
        self.tb_oleaje.setHorizontalHeaderLabels(["H (m)", "T (s)", "Frecuencia (%)"])
        self.tb_oleaje.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_oleaje.itemChanged.connect(self._cambio_oleaje)
        v.addWidget(self.tb_oleaje)
        self.graf_oleaje = Lienzo()
        self.graf_oleaje.setMinimumHeight(220)
        v.addWidget(self.graf_oleaje)
        self.lb_oleaje = QLabel("", objectName="note")
        self.lb_oleaje.setWordWrap(True)
        v.addWidget(self.lb_oleaje)
        tabs.addTab(w, "Oleaje")

        # --- Nivel / mareas
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        fila.addWidget(InfoBoton("serie_nivel"))
        b = QPushButton("Cargar serie (marea.dat)…")
        b.clicked.connect(self.cargar_serie_nivel)
        fila.addWidget(b)
        b = QPushButton("Quitar")
        b.clicked.connect(self.quitar_serie_nivel)
        fila.addWidget(b)
        fila.addStretch(1)
        v.addLayout(fila)

        g = QGroupBox("Variación del nivel del lago en el cálculo")
        self.grp_nivel_lago = g
        gf = QVBoxLayout(g)
        self.ck_estacional = QCheckBox("Variación estacional (subida respecto al nivel del lago):")
        fe0 = QHBoxLayout()
        fe0.addWidget(self.ck_estacional)
        fe0.addWidget(InfoBoton("nivel_estacional"))
        fe0.addStretch(1)
        gf.addLayout(fe0)
        fe = QHBoxLayout()
        self.sp_estac = {}
        for clave, nombre in ESTACIONES:
            fe.addWidget(QLabel(nombre))
            sp = dspin(-5, 5, 2, 0.1, "m")
            sp.setValue(gl.NIVEL_ESTACIONAL_TESIS[clave])
            sp.valueChanged.connect(self._cambio_opciones_nivel)
            self.sp_estac[clave] = sp
            fe.addWidget(sp)
        fe.addStretch(1)
        gf.addLayout(fe)
        self.ck_interanual = QCheckBox("Ciclo interanual de 15 años (±5·Δz, cambios cada 9 meses "
                                       "desde el mes 69)")
        fi0 = QHBoxLayout()
        fi0.addWidget(self.ck_interanual)
        fi0.addWidget(InfoBoton("nivel_interanual"))
        fi0.addStretch(1)
        gf.addLayout(fi0)
        fb = QHBoxLayout()
        b = QPushButton("Valores de la tesis")
        b.clicked.connect(self._nivel_tesis)
        fb.addWidget(b)
        fb.addStretch(1)
        gf.addLayout(fb)
        self.lb_aviso_nivel = QLabel("", objectName="note")
        self.lb_aviso_nivel.setWordWrap(True)
        gf.addWidget(self.lb_aviso_nivel)
        for ck in (self.ck_estacional, self.ck_interanual):
            ck.toggled.connect(self._cambio_opciones_nivel)
        v.addWidget(g)

        gm_ = QGroupBox("Marea y nivel del mar (modo mar)")
        self.grp_marea = gm_
        fm = QFormLayout(gm_)
        self.sp_amp = dspin(0, 10, 2, 0.1, "m")
        self.sp_amp.setValue(gl.MAR_TESIS["amplitud"])
        self.sp_subida = dspin(-50, 50, 1, 0.5, "mm/año")
        self.sp_subida.setValue(gl.MAR_TESIS["subida_mm_ano"])
        self.sp_nmarea = QSpinBox()
        self.sp_nmarea.setRange(1, 81)
        self.sp_nmarea.setValue(gl.MAR_TESIS["niveles"])
        fm.addRow(etiqueta_info("Amplitud de la marea", "mar_amplitud"), self.sp_amp)
        fm.addRow(etiqueta_info("Subida del nivel del mar", "mar_subida"), self.sp_subida)
        fm.addRow(etiqueta_info("Niveles de marea", "mar_niveles"), self.sp_nmarea)
        b = QPushButton("Valores por defecto (Holderness, tesis cap. 9)")
        b.clicked.connect(self._mar_defecto)
        fm.addRow("", b)
        self.lb_aviso_mar = QLabel(
            "El modo mar es una adaptación de la herramienta: no existe un Fortran con el "
            "que comparar. Con amplitud 0, sin subida y el oleaje del lago reproduce "
            "exactamente el modo lago.", objectName="note")
        self.lb_aviso_mar.setWordWrap(True)
        fm.addRow(self.lb_aviso_mar)
        for w_ in (self.sp_amp, self.sp_subida):
            w_.valueChanged.connect(self._cambio_opciones_nivel)
        self.sp_nmarea.valueChanged.connect(self._cambio_opciones_nivel)
        v.addWidget(gm_)
        gm_.setVisible(False)

        self.graf_nivel = Lienzo()
        v.addWidget(self.graf_nivel, 1)
        self.lb_nivel = QLabel("", objectName="note")
        self.lb_nivel.setWordWrap(True)
        v.addWidget(self.lb_nivel)
        tabs.addTab(w, "Nivel")
        tabs.setTabToolTip(tabs.count() - 1, "Nivel del lago / mareas")

        self.tabs_der = tabs
        dock.setWidget(tabs)
        dock.setMinimumWidth(560)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.mostrar_oleaje()
        self._cambio_opciones_nivel()
        self._cambio_estabilidad()

    # ------------------------------------------------------- datos <-> UI
    def cargar_datos(self, d):
        self._bloqueo = True
        try:
            self.sp_tsim.setValue(int(d["tsim"]))
            self.cb_utsim.setCurrentIndex(1 if int(d["Utsim"]) == 2 else 0)
            self.sp_tdib.setValue(int(d["Tdib"]))
            self.cb_utdib.setCurrentIndex(1 if int(d["UTdib"]) == 2 else 0)
            self.sp_lake.setValue(int(d["lakelevel"]))
            self.sp_dz.setValue(float(d["Dz"]))
            for clave, w in self.campos.items():
                if clave == "groundwater":
                    w.setCurrentIndex(int(d["groundwater"]))
                elif clave == "rn":
                    w.setValue(int(d["rn"]))
                else:
                    w.setValue(float(d[clave]))
            self.tb_perfil.setRowCount(0)
            for z, y, e in d["Pzy"]:
                self._fila_perfil(float(z), float(y), int(e))
            self.tb_estratos.setRowCount(0)
            names = ["sigmac", "Shearstrength", "Ecohesion", "Ephi", "pesounit", "potestrato"]
            for i in range(1, int(d["ndatosresis"]) + 1):
                self._fila_estrato([float(d[nm][i]) for nm in names])
        finally:
            self._bloqueo = False
        self._cambio_datos()

    def _fila_perfil(self, z, y, e, fila=None):
        t = self.tb_perfil
        if fila is None:
            fila = t.rowCount()
        t.insertRow(fila)
        t.setItem(fila, 0, QTableWidgetItem(_fmt(z)))
        t.setItem(fila, 1, QTableWidgetItem(_fmt(y)))
        cb = QComboBox()
        for k in sorted(ETIQUETAS):
            cb.addItem(ETIQUETAS[k], k)
        cb.setCurrentIndex(e if e in ETIQUETAS else 0)
        cb.currentIndexChanged.connect(self._cambio_datos)
        t.setCellWidget(fila, 2, cb)

    def _fila_estrato(self, vals, fila=None):
        t = self.tb_estratos
        if fila is None:
            fila = t.rowCount()
        t.insertRow(fila)
        for k, x in enumerate(vals):
            t.setItem(fila, k, QTableWidgetItem(_fmt(x)))

    @staticmethod
    def _num(item):
        if item is None:
            raise ValueError("celda vacía")
        return float(item.text().replace(",", "."))

    def leer_datos(self):
        """Datos actuales de la interfaz como diccionario (formato read_data)."""
        d = dict(rn=self.campos["rn"].value(), tsim=self.sp_tsim.value(),
                 Utsim=2 if self.cb_utsim.currentIndex() == 1 else 1,
                 Tdib=self.sp_tdib.value(),
                 UTdib=2 if self.cb_utdib.currentIndex() == 1 else 1,
                 Dz=gl.F64(self.sp_dz.value()), lakelevel=self.sp_lake.value(),
                 groundwater=self.campos["groundwater"].currentIndex())
        for clave in ("SFmin", "Kss", "p", "WFmin", "MMM", "Rsc", "Cb", "KRs", "Bmat",
                      "Qb", "Gs", "profagua", "Mf"):
            d[clave] = gl.F64(self.campos[clave].value())
        pts = []
        for r in range(self.tb_perfil.rowCount()):
            z = self._num(self.tb_perfil.item(r, 0))
            y = self._num(self.tb_perfil.item(r, 1))
            e = self.tb_perfil.cellWidget(r, 2).currentData()
            pts.append([gl.F64(z), gl.F64(y), gl.F64(e)])
        d["Pzy"] = pts
        n = self.tb_estratos.rowCount()
        d["ndatosresis"] = n
        names = ["sigmac", "Shearstrength", "Ecohesion", "Ephi", "pesounit", "potestrato"]
        for k, nm in enumerate(names):
            arr = np.zeros(n + 1)
            for r in range(n):
                arr[r + 1] = self._num(self.tb_estratos.item(r, k))
            d[nm] = arr
        return d

    def validar(self, d):
        errores = []
        pts = d["Pzy"]
        if len(pts) < 3:
            errores.append("El perfil necesita al menos 3 puntos.")
        zs = [p[0] for p in pts]
        if any(zs[i + 1] > zs[i] for i in range(len(zs) - 1)):
            errores.append("Las cotas del perfil deben ir de mayor a menor.")
        etiquetas = [int(p[2]) for p in pts[1:]]
        if 2 not in etiquetas:
            errores.append("Falta la cresta del talud (punto 2).")
        if 4 not in etiquetas:
            errores.append("Falta la base del talud (punto 4).")
        if d["groundwater"] == 1 and 3 not in etiquetas:
            errores.append("Con nivel freático (tipo 1) hace falta el punto de rezume (3).")
        if etiquetas.count(1) > 10:
            errores.append("Como máximo 10 puntos de estudio.")
        if d["ndatosresis"] < 1:
            errores.append("Hace falta al menos un estrato.")
        if pts and d["ndatosresis"] >= 1:
            altura = pts[0][0] - pts[-1][0]
            if sum(d["potestrato"][1:]) < altura:
                errores.append("Los estratos no cubren toda la altura del perfil "
                               f"({altura:.1f} m).")
        if pts:
            Lmax = max(gl.f32(1.56) * gl.f32(t) ** 2 for s in self.seasons.values()
                       for t in s[1])
            fondo = pts[-1][0]
            amp = self.sp_amp.value() if self.modo == "mar" else 0.0
            necesario = Lmax / 2 + 1 + amp
            if d["lakelevel"] - fondo < necesario:
                que = ("del nivel medio del mar (bajamar + media longitud de onda de la ola "
                       "más larga)" if self.modo == "mar" else
                       "del nivel del lago (media longitud de onda de la ola más larga)")
                errores.append(
                    f"El perfil debe llegar al menos a {necesario:.1f} m por debajo {que}; "
                    f"ahora llega a {d['lakelevel'] - fondo:.1f} m.")
        if self.sp_tdib.value() * (12 if d["UTdib"] == 2 else 1) > \
                self.sp_tsim.value() * (12 if d["Utsim"] == 2 else 1):
            errores.append("El intervalo de salida es mayor que la duración.")
        return errores

    def _cambio_datos(self, *_):
        if self._bloqueo:
            return
        self.dibujar_perfil()
        self.dibujar_nivel()

    # --------------------------------------------------------------- dibujo
    def _perfil_preliminar(self):
        try:
            d = self.leer_datos()
        except ValueError:
            return None, None
        if len(d["Pzy"]) < 2 or d["Dz"] <= 0:
            return d, None
        try:
            g = gl.prepare_profile(d["Pzy"], d["Dz"], d["lakelevel"])
        except Exception:  # noqa: BLE001
            return d, None
        self.lb_npts.setText(f"{g['N']} puntos")
        return d, g

    def dibujar_perfil(self, *_):
        lz = self.lienzo
        lz.fig.clf()
        lz.ax = lz.fig.add_subplot(111)
        ax = lz.ax
        lz.estilo()
        d, g = self._perfil_preliminar()
        if d is None or g is None:
            ax.text(0.5, 0.5, "Perfil no válido", transform=ax.transAxes, ha="center",
                    color=CANVAS_FG)
            lz.draw_idle()
            return
        N = g["N"]
        P = g["P"]
        if self.res is not None and self.res["meta"]["N"] == N:
            P = self.res["profile_inicial"]
        z = P[1, 1:N + 1]
        y0 = P[2, 1:N + 1]
        zlago = P[1, g["TL"]] if 1 <= g["TL"] <= N else d["lakelevel"]

        ymin, ymax = float(np.min(y0)), float(np.max(y0))
        zmin, zmax = float(np.min(z)), float(np.max(z))
        # Agua y terreno
        ax.fill_betweenx(z, y0, ymax + 1e6, where=z <= zlago + 1e-9, color="#2e6f9e",
                         alpha=0.35, lw=0, label="_agua")
        modo_mes = (self.res is not None and self.res["meta"]["N"] == N
                    and self.linea.valor > 0)
        self._modo_activo = modo_mes
        terreno = y0
        if self.res is not None and self.res["meta"]["N"] == N and not modo_mes:
            vis = [r for r in range(self.tb_salidas.rowCount())
                   if self.tb_salidas.item(r, 0) is not None
                   and self.tb_salidas.item(r, 0).checkState() == Qt.Checked]
            if vis:
                terreno = self.res["DibMatriz"][2:N + 2, vis[-1] + 3]
        self._terreno_base = terreno
        self._fill_terreno = ax.fill_betweenx(z, ymin - 1e6, terreno, color="#6b5a45",
                                              alpha=0.55, lw=0, zorder=0)
        self._lineas_tenues = []
        self._art_mes = []
        self._z_actual = z
        self._ymin_actual = ymin
        self._y0_actual = y0
        self._zlago_actual = zlago
        self._fill_erosion = []
        if self.modo == "mar":
            amp = self.sp_amp.value()
            ax.axhline(zlago, color="#4fa3e0", lw=1.2, ls="-",
                       label=f"Nivel medio del mar ({zlago:.2f} m)")
            if amp > 0:
                ax.axhspan(zlago - amp, zlago + amp, color="#4fa3e0", alpha=0.10, lw=0)
                ax.axhline(zlago + amp, color="#4fa3e0", lw=0.9, ls=":",
                           label=f"Pleamar / bajamar (±{amp:.1f} m)")
                ax.axhline(zlago - amp, color="#4fa3e0", lw=0.9, ls=":")
        else:
            ax.axhline(zlago, color="#4fa3e0", lw=1.2, ls="-",
                       label=f"Nivel del lago ({zlago:.2f} m)")
        ax.plot(y0, z, color="#f0f0f0", lw=1.8, zorder=4, label="Perfil inicial (año 0, datos)")
        ax.plot([p[1] for p in d["Pzy"]], [p[0] for p in d["Pzy"]], "o", ms=3.5,
                color="#bbbbbb", alpha=0.8, label="_vertices")

        # Nivel freatico
        if d["groundwater"] == 1:
            yf, zf = curva_freatica(P, g["seepagelevel"], d["profagua"])
            if yf is not None:
                ax.plot(yf, zf, ls="--", color="#5fd1e6", lw=1.2, label="Nivel freático")
        elif d["groundwater"] == 2:
            ax.text(0.01, 0.02, "Terreno saturado", transform=ax.transAxes, color="#5fd1e6",
                    fontsize=8)

        # Resultados
        hay_res = self.res is not None and self.res["meta"]["N"] == N
        if hay_res:
            Dib = self.res["DibMatriz"]
            sel = sorted({i.row() for i in self.tb_roturas.selectedItems()})
            cols = [r for r in range(self.tb_salidas.rowCount())
                    if self.tb_salidas.item(r, 0) is not None
                    and self.tb_salidas.item(r, 0).checkState() == Qt.Checked]
            cmap = matplotlib.colormaps["plasma"]
            anos = Dib[1, 3:][Dib[1, 3:] > 0] / 12.0
            amin = float(anos.min()) if anos.size else 0.0
            amax = float(anos.max()) if anos.size else 1.0
            if amax <= amin:
                amin, amax = amin - 0.5, amax + 0.5
            norm = matplotlib.colors.Normalize(amin, amax)
            pp = self.res.get("primer_paso")
            if pp is not None and self.ck_primer.isChecked():
                mes_pp, perfil_pp = pp
                ln, = ax.plot(perfil_pp, z, color="#9fb3c8", lw=1.3, ls=(0, (4, 2)),
                              alpha=0.22 if modo_mes else (0.4 if sel else 1.0),
                              label=f"Tras el 1.er mes calculado (mes {mes_pp})")
                self._lineas_tenues.append(ln)
            for r in cols:
                k = r + 3
                ano = Dib[1, k] / 12.0
                if ano <= 0:
                    continue
                ln, = ax.plot(Dib[2:N + 2, k], z, color=cmap(norm(ano)), lw=1.2,
                              alpha=0.22 if modo_mes else (0.3 if sel else 1.0))
                self._lineas_tenues.append(ln)
            if cols:
                sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
                cb = lz.fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.04)
                cb.set_label("Perfiles de salida: año", color=CANVAS_FG, fontsize=8)
                cb.ax.tick_params(colors=CANVAS_FG, labelsize=7)
                cb.outline.set_edgecolor("#5a6270")
            for r in sel:
                if r >= len(self._eventos):
                    continue
                mes, tipo, vol, antes, despues, _ = self._eventos[r]
                c1, c2 = (("#ff6b6b", "#ff3030") if tipo == "Deslizamiento"
                          else ("#ffc966", COLOR_VUELCO))
                ax.plot(antes, z, color=c1, lw=1.1, ls="--",
                        label=f"Antes: {tipo.lower()} mes {mes}")
                ax.plot(despues, z, color=c2, lw=1.6,
                        label=f"Después: {tipo.lower()} mes {mes} ({vol:.1f} m³/m)")

        # Puntos singulares
        marcas = [(g["crestatalud"], "^", "#ffcc33", "Cresta del talud"),
                  (g["roturanotch"], "s", "#ff9f43", "Base del talud")]
        if d["groundwater"] == 1:
            marcas.append((g["seepagelevel"], "v", "#5fd1e6", "Punto de rezume"))
        for k, m, c, nombre in marcas:
            if 1 <= k <= N:
                ax.plot(P[2, k], P[1, k], m, color=c, ms=8, mec="black", mew=0.6, label=nombre)
        for i, k in enumerate(g["label"][1:g["nlabel"] + 1]):
            if 1 <= k <= N:
                ax.plot(P[2, k], P[1, k], "o", color="#7bed9f", ms=7, mec="black", mew=0.6,
                        label="Punto de estudio" if i == 0 else "_")

        ax.set_xlabel("Distancia horizontal y (m)")
        ax.set_ylabel("Cota z (m)")
        if self.cb_vista.currentIndex() == 0:
            kc = max(1, g["crestatalud"])
            kb = min(N, max(1, g["roturanotch"]))
            kt = min(N, max(kb, g["TL"]))
            H = max(1.0, abs(P[1, kc] - P[1, kb]))
            tramo = slice(kc, kt + 1)
            ys = [float(np.min(P[2, tramo]))]
            # mismo encuadre para todos los métodos calculados, para que al
            # cambiar de resultado no cambie la escala horizontal
            for r_ in ([x["res"] for x in self.resultados] if hay_res else []):
                if r_["meta"]["N"] == N and r_["DibMatriz"].shape[1] > 3:
                    ys.append(float(np.min(r_["DibMatriz"][kc + 1:kt + 2, 3:])))
            x0 = min(ys) - 0.5 * H
            x1 = float(np.max(P[2, kc:kb + 1])) + 1.5 * H
            ax.set_xlim(x0, x1)
            ax.set_ylim(zlago - 0.5 * H - (self.sp_amp.value() if self.modo == "mar" else 0),
                        zmax + 0.15 * H)
        else:
            ax.set_xlim(ymin, ymax)
            ax.set_ylim(zmin, zmax + 0.05 * (zmax - zmin))
        if self.ck_leyenda.isChecked():
            leg = ax.legend(loc="best", fontsize=7.5, facecolor="#2f3540",
                            edgecolor="#5a6270", labelcolor=CANVAS_FG)
            leg.set_draggable(True)
        lz.fig.tight_layout()
        if hay_res:
            self._linea_mes, = ax.plot(y0, z, color="#00e5ff", lw=2.6, zorder=6, visible=False,
                                       animated=True)
            self._txt_mes = ax.text(0.015, 0.975, "", transform=ax.transAxes, va="top",
                                    ha="left", fontsize=10, fontweight="bold", color="#00e5ff",
                                    zorder=7, animated=True,
                                    bbox=dict(boxstyle="round,pad=0.35",
                                              fc="#1d2128", ec="#00e5ff", lw=0.8))
            self._txt_mes.set_visible(False)
            self._linea_nivel_mes = ax.axhline(zlago, color="#7cc4ff", lw=1.6, ls="--",
                                               zorder=6, visible=False, animated=True)
            self._actualizar_mes(redibujar=False)
        else:
            self._linea_mes = None
        lz.draw_idle()

    # --------------------------------------------------------- resultados
    def _series(self):
        """Series temporales a partir de los perfiles de salida."""
        r = self.res
        N = r["meta"]["N"]
        Dib = r["DibMatriz"]
        P0 = r["profile_inicial"]
        meses = Dib[1, 3:]
        validos = meses > 0
        cols = np.arange(3, Dib.shape[1])[validos]
        meses = meses[validos]
        TL = r["meta"]["TL"]
        kc = r["meta"]["crestatalud"]
        lago = np.array([P0[2, TL] - Dib[TL + 1, c] for c in cols])
        cresta = np.array([P0[2, kc] - Dib[kc + 1, c] for c in cols])
        estudio = [np.array([P0[2, k] - Dib[k + 1, c] for c in cols]) for k in r["meta"]["label"]]
        return meses, lago, cresta, estudio, N

    def mostrar_resultados(self):
        r = self.res
        meses, lago, cresta, _, _ = self._series()
        t = self.tb_salidas
        t.blockSignals(True)
        t.setRowCount(0)
        for i, mes in enumerate(meses):
            t.insertRow(i)
            it = QTableWidgetItem()
            it.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            it.setCheckState(Qt.Checked)
            t.setItem(i, 0, it)
            for k, val in enumerate((f"{mes:.0f}", f"{mes / 12:.1f}", _f2(lago[i]),
                                     _f2(cresta[i])), start=1):
                c = QTableWidgetItem(val)
                c.setFlags(Qt.ItemIsEnabled)
                c.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                t.setItem(i, k, c)
        t.blockSignals(False)

        # deslizamientos circulares y vuelcos, ordenados por mes
        fss = r.get("roturas_fs", [])
        ev = [(m, 0, "Deslizamiento", v, a, b, fss[i] if i < len(fss) else None)
              for i, (m, v, a, b) in enumerate(r["roturas"])]
        ev += [(m, 1, "Vuelco", v, a, b, None) for (m, v, a, b) in r.get("vuelcos", [])]
        ev.sort(key=lambda e: (e[0], e[1]))
        self._eventos = [(m, t, v, a, b, f) for m, _, t, v, a, b, f in ev]
        tr = self.tb_roturas
        tr.setRowCount(0)
        for i, (mes, tipo, vol, _, _, fs) in enumerate(self._eventos):
            tr.insertRow(i)
            fs_txt = f"{fs:.3f}" if fs is not None else "—"
            for k, val in enumerate((f"{mes}", f"{mes / 12:.1f}", tipo, f"{vol:.3f}", fs_txt)):
                c = QTableWidgetItem(val)
                c.setTextAlignment((Qt.AlignLeft if k == 2 else Qt.AlignRight) | Qt.AlignVCenter)
                if tipo == "Vuelco":
                    c.setForeground(QColor("#9a6400"))
                tr.setItem(i, k, c)
        n_d = len(r["roturas"])
        n_v = len(r.get("vuelcos", []))
        self.tabs_der.setTabText(2, f"Roturas ({n_d}+{n_v})")
        self.tabs_der.setTabToolTip(2, f"{n_d} deslizamientos y {n_v} vuelcos")
        self._set_invierno = set(r["mensual"]["invierno"])
        self._set_derrubios = set(r["mensual"]["derrubios"])
        self._roturas_mes = {}
        for e in self._eventos:
            self._roturas_mes.setdefault(e[0], []).append(e)
        self.linea.texto_invierno = ("invierno: lago helado, sin oleaje ni erosión"
                                     if r.get("modo") != "mar" else "estación sin oleaje")
        self.linea.configurar(r["meta"]["tsim"], [m for m, *_ in r["roturas"]],
                              r["mensual"]["invierno"], r["mensual"]["derrubios"],
                              [m for m, *_ in r.get("vuelcos", [])])
        self._activar_linea(True)
        self.dibujar_perfil()
        self.dibujar_grafica()
        self.dibujar_nivel()

    # ------------------------------------------------------ linea temporal
    def _activar_linea(self, activa):
        for b in (self.b_inicio, self.b_atras, self.b_play, self.b_adelante, self.b_fin,
                  self.cb_vel):
            b.setEnabled(activa)
        if not activa:
            self.timer.stop()
            self.b_play.setText("▶")
            self.linea.configurar(0, [], [], [])
            self.lb_mes.setText("")

    def _velocidad(self, *_):
        self.timer.setInterval([300, 120, 50, 16][self.cb_vel.currentIndex()])

    def reproducir(self):
        if self.timer.isActive():
            self.timer.stop()
            self.b_play.setText("▶")
            return
        if self.linea.valor >= self.linea.tsim:
            self.linea.set_valor(0)
        self.b_play.setText("⏸")
        self.timer.start()

    def _tick(self):
        if self.linea.valor >= self.linea.tsim:
            self.timer.stop()
            self.b_play.setText("▶")
            return
        self.linea.set_valor(self.linea.valor + 1)

    def _cambio_mes(self, _):
        self._actualizar_mes()

    def _perfil_mes(self, mes):
        """Perfil (posiciones y) al final del mes indicado; mes 0 = inicial."""
        r = self.res
        men = r["mensual"]
        i = int(np.searchsorted(men["meses"], mes, side="right")) - 1
        if mes <= 0 or i < 0:
            return r["profile_inicial"][2, 1:r["meta"]["N"] + 1]
        return men["perfiles"][i]

    def _animados(self):
        return (self._fill_erosion + self._art_mes
                + [self._linea_nivel_mes, self._linea_mes, self._txt_mes])

    def _on_draw(self, _evento):
        """Tras un dibujo completo guarda el fondo y pinta los elementos del mes."""
        if self._linea_mes is None:
            self._bg = None
            return
        lz = self.lienzo
        self._bg = lz.copy_from_bbox(lz.fig.bbox)
        for a in self._animados():
            if a.get_visible():
                lz.ax.draw_artist(a)

    def _actualizar_mes(self, redibujar=True):
        """Actualiza el perfil del mes seleccionado (dibujo rapido por blitting)."""
        if self.res is None or self._linea_mes is None:
            return
        mes = self.linea.valor
        activo = mes > 0
        if redibujar and activo != self._modo_activo:
            # cambia el fondo (perfiles atenuados / terreno): dibujo completo
            self.dibujar_perfil()
            return
        ax = self.lienzo.ax
        for a in self._art_mes + self._fill_erosion:
            a.remove()
        self._art_mes = []
        self._fill_erosion = []
        self._linea_mes.set_visible(activo)
        self._txt_mes.set_visible(activo)
        z = self._z_actual
        meta = self.res["meta"]
        nivel_var = (self.res.get("nivel_lago") is not None
                     or (self.res.get("modo") == "mar"
                         and (self.res.get("mar") or {}).get("subida_mm_ano", 0) != 0))
        self._linea_nivel_mes.set_visible(activo and nivel_var)
        txt_nivel = ""
        if activo and nivel_var:
            men = self.res["mensual"]
            i = int(np.searchsorted(men["meses"], mes, side="right")) - 1
            if i >= 0:
                zl = float(men["nivel"][i])
                self._linea_nivel_mes.set_ydata([zl, zl])
                txt_nivel = (f" · nivel medio del mar {zl:+.2f} m" if self.res.get("modo") == "mar"
                             else f" · nivel del lago {zl:+.2f} m")
        if activo:
            perfil = self._perfil_mes(mes)
            self._linea_mes.set_xdata(perfil)
            # zona erosionada hasta este mes: se "vacia" (aire o agua)
            y0 = self._y0_actual
            erosion = perfil < y0 - 1e-9
            for donde, color in ((erosion & (z > self._zlago_actual), CANVAS_BG),
                                 (erosion & (z <= self._zlago_actual), "#29435a")):
                if donde.any():
                    self._fill_erosion.append(ax.fill_betweenx(
                        z, perfil, y0, where=donde, color=color, alpha=0.9, lw=0,
                        zorder=1, animated=True))
            rot = self._roturas_mes.get(mes, [])
            for e in rot:
                ln, = ax.plot(e[3], z, color="#ff6b6b" if e[1] == "Deslizamiento" else "#ffc966",
                              lw=1.0, ls="--", zorder=5, animated=True)
                self._art_mes.append(ln)
            P0 = self.res["profile_inicial"]
            TL, kc = meta["TL"], meta["crestatalud"]
            r_lago = P0[2, TL] - perfil[TL - 1]
            r_cresta = P0[2, kc] - perfil[kc - 1]
            if mes in self._set_invierno:
                estado = self.linea.texto_invierno
            elif mes in self._set_derrubios:
                estado = "protegido por derrubios: sin erosión"
            else:
                estado = "calculado"
            txt = f"Mes {mes}  ·  año {mes / 12:.2f}"
            nd = sum(1 for e in rot if e[1] == "Deslizamiento")
            nv = len(rot) - nd
            partes = []
            if nd:
                partes.append(f"{nd} deslizamiento{'s' if nd > 1 else ''}")
            if nv:
                partes.append(f"{nv} vuelco{'s' if nv > 1 else ''}")
            if rot:
                txt += "\n" + ", ".join(partes) + f" ({sum(e[2] for e in rot):.1f} m³/m)"
            self._txt_mes.set_text(txt)
            nom_nivel = "nivel del mar" if self.res.get("modo") == "mar" else "nivel del lago"
            self.lb_mes.setText(
                f"Mes {mes} de {meta['tsim']} (año {mes / 12:.2f}) · {estado} · "
                f"retroceso: {nom_nivel} {_f2(r_lago)} m, cresta {_f2(r_cresta)} m"
                + (" · " + ", ".join(partes) if rot else "") + txt_nivel)
        else:
            self.lb_mes.setText(f"Mes 0 (perfil inicial) · arrastra la barra para ver la "
                                f"evolución mes a mes hasta el mes {meta['tsim']}")
        if not redibujar:
            return
        lz = self.lienzo
        if self._bg is None:
            lz.draw_idle()
            return
        lz.restore_region(self._bg)
        for a in self._animados():
            if a.get_visible():
                ax.draw_artist(a)
        lz.blit(lz.fig.bbox)

    def _marcar_salidas(self, estado):
        t = self.tb_salidas
        t.blockSignals(True)
        for r in range(t.rowCount()):
            t.item(r, 0).setCheckState(Qt.Checked if estado else Qt.Unchecked)
        t.blockSignals(False)
        self.dibujar_perfil()

    def _cambio_salidas(self, item):
        if item.column() == 0:
            self.dibujar_perfil()

    def dibujar_grafica(self, *_):
        lz = self.graf
        ax = lz.ax
        ax.clear()
        lz.estilo()
        if self.res is None:
            ax.text(0.5, 0.5, "Sin resultados: pulsa «Calcular»", transform=ax.transAxes,
                    ha="center", color="#888")
            ax.set_axis_off()
            lz.draw_idle()
            return
        meses, lago, cresta, estudio, _ = self._series()
        anos = meses / 12.0
        modo = self.cb_grafica.currentIndex()
        if modo == 0:
            ax.plot(anos, lago, "-o", color=ACCENT, ms=3.5)
            ax.set_ylabel("Retroceso en el nivel del lago (m)")
        elif modo == 1:
            ax.plot(anos, cresta, "-o", color=ACCENT, ms=3.5)
            ax.set_ylabel("Retroceso de la cresta (m)")
        elif modo == 2:
            if not estudio:
                ax.text(0.5, 0.5, "No hay puntos de estudio en el perfil",
                        transform=ax.transAxes, ha="center", color="#888")
            for i, s in enumerate(estudio):
                ax.plot(anos, s, "-o", ms=3, label=f"Punto de estudio {i + 1}")
            if estudio:
                ax.legend(fontsize=8)
            ax.set_ylabel("Retroceso (m)")
        else:
            ev = getattr(self, "_eventos", [])
            if ev:
                t = np.array([e[0] for e in ev]) / 12.0
                v = np.cumsum([e[2] for e in ev])
                ax.step(t, v, where="post", color=ACCENT, label="Total")
                ax.plot(t, v, "o", color=ACCENT, ms=3)
                for tipo, c in (("Deslizamiento", "#d0312d"), ("Vuelco", "#c98a00")):
                    sub = [e for e in ev if e[1] == tipo]
                    if sub:
                        ax.step(np.array([e[0] for e in sub]) / 12.0,
                                np.cumsum([e[2] for e in sub]), where="post", color=c,
                                lw=1, ls="--", label=tipo + "s")
                ax.legend(fontsize=8)
            else:
                ax.text(0.5, 0.5, "No se han producido roturas",
                        transform=ax.transAxes, ha="center", color="#888")
            ax.set_ylabel("Volumen acumulado (m³/m)")
        ax.set_xlabel("Tiempo (años)")
        lz.fig.tight_layout()
        lz.draw_idle()

    # ------------------------------------------------------------- oleaje
    def _clave_estacion(self):
        return ESTACIONES[self.cb_estacion.currentIndex()][0]

    def mostrar_oleaje(self, *_):
        H, T, WF = self.seasons[self._clave_estacion()]
        t = self.tb_oleaje
        t.blockSignals(True)
        t.setRowCount(0)
        for i, (h, p, f) in enumerate(zip(H, T, WF)):
            t.insertRow(i)
            for k, x in enumerate((h, p, f)):
                t.setItem(i, k, QTableWidgetItem(f"{x:g}"))
        t.blockSignals(False)
        self._dibujar_oleaje()

    def _cambio_oleaje(self, *_):
        t = self.tb_oleaje
        try:
            filas = [[self._num(t.item(r, k)) for k in range(3)] for r in range(t.rowCount())]
        except ValueError:
            return
        self.seasons[self._clave_estacion()] = tuple(list(c) for c in zip(*filas)) \
            if filas else ([0], [0], [0])
        self._dibujar_oleaje()

    def anadir_ola(self):
        t = self.tb_oleaje
        r = t.rowCount()
        t.insertRow(r)
        for k, x in enumerate((1.0, 4.0, 1.0)):
            t.setItem(r, k, QTableWidgetItem(f"{x:g}"))

    def eliminar_ola(self):
        filas = sorted({i.row() for i in self.tb_oleaje.selectedItems()}, reverse=True)
        for r in filas:
            self.tb_oleaje.removeRow(r)
        self._cambio_oleaje()

    def _oleaje_temporales(self):
        self.seasons = copy.deepcopy(gl.SEASONS_MAR_TEMPORALES)
        self.mostrar_oleaje()

    def restablecer_oleaje(self):
        self.seasons = copy.deepcopy(self._seasons_defecto())
        self.mostrar_oleaje()

    def _dibujar_oleaje(self):
        lz = self.graf_oleaje
        ax = lz.ax
        ax.clear()
        lz.estilo()
        H, T, WF = self.seasons[self._clave_estacion()]
        if len(H) and any(T):
            ax.bar(range(len(H)), WF, color=ACCENT, alpha=0.85)
            ax.set_xticks(range(len(H)), [f"H={h:g} m\nT={p:g} s" for h, p in zip(H, T)],
                          fontsize=7)
            ax.set_ylabel("Frecuencia (%)")
        else:
            ax.text(0.5, 0.5, "Sin oleaje en esta estación" + (" (el invierno no erosiona)"
                                                                if self.modo == "lago" else ""),
                    transform=ax.transAxes, ha="center", color="#888")
        lz.fig.tight_layout()
        lz.draw_idle()
        modificado = self.seasons != self._seasons_defecto()
        if self.modo == "lago":
            base = ("Oleaje por estación usado por el modelo (Trenhaile, 2009), lago Erie. "
                    "Primavera: meses 4-6, verano: 7-9, otoño: 10-12; los meses de invierno "
                    "(1-3) no se calculan (hielo).")
            aviso = ("\n⚠ Valores modificados respecto al código original: los resultados ya "
                     "no coincidirán con el Fortran.")
        else:
            base = ("Modo mar: se calculan los 12 meses con el oleaje de su estación "
                    "(invierno: ene-mar). Por defecto, la ola media de Holderness que usa la "
                    "tesis (0,98 m, 7,66 s), adecuada para la erosión de Castedo. La de "
                    "Trenhaile solo erosiona con olas grandes: para ella usa «Con temporales» "
                    "o, mejor, datos de una boya cercana.")
            aviso = "\n⚠ Valores modificados respecto a los de por defecto del modo mar."
        self.lb_oleaje.setText(base + (aviso if modificado else ""))
        self.lb_oleaje.setStyleSheet("color: #b0451d;" if modificado else "")

    # -------------------------------------------------------- nivel/mareas
    def cargar_serie_nivel(self):
        ruta, _ = QFileDialog.getOpenFileName(self, "Serie de nivel / mareas",
                                              self._carpeta(), "Datos (*.dat *.txt *.csv);;Todos (*)")
        if not ruta:
            return
        t, z = [], []
        with open(ruta, encoding="latin-1") as fh:
            for linea in fh:
                toks = linea.replace(",", " ").replace(";", " ").split()
                try:
                    vals = [float(x) for x in toks]
                except ValueError:
                    continue
                if len(vals) >= 2:
                    t.append(vals[0])
                    z.append(vals[1])
                elif len(vals) == 1:
                    t.append(len(t) + 1)
                    z.append(vals[0])
        if not t:
            QMessageBox.information(self, "Serie de nivel",
                                    f"«{os.path.basename(ruta)}» está vacío o no contiene datos "
                                    "numéricos.\nFormato esperado: una o dos columnas "
                                    "(tiempo, nivel).")
            return
        self.serie_nivel = (os.path.basename(ruta), np.array(t), np.array(z))
        self.dibujar_nivel()

    def quitar_serie_nivel(self):
        self.serie_nivel = None
        self.dibujar_nivel()

    def metodos_elegidos(self):
        ero = [k for k, ck in self.ck_ero.items() if ck.isChecked()]
        est = [k for k, ck in self.ck_est.items() if ck.isChecked()]
        return ero, est

    def opciones_estabilidad(self, metodo="fellenius"):
        """Opciones de estabilidad para gl.run (None = como el original)."""
        op = {"metodo": metodo, "rebanadas": self.sp_rebanadas.value(),
              "criterio": self.cb_criterio.currentData()}
        if op == {"metodo": "fellenius", "rebanadas": 5, "criterio": "todas"}:
            return None
        return op

    def opciones_erosion(self, metodo):
        """Opciones de erosion para gl.run (None = Trenhaile, como el original)."""
        if metodo != "castedo":
            return None
        return {"metodo": "castedo", "K": self.sp_K.value(),
                "pendiente_max": self.sp_pmax.value(),
                "derrubios_deslizamiento": self.sp_fdesl.value() / 100.0,
                "derrubios_vuelco": self.sp_fvuel.value() / 100.0}

    def _castedo_tesis(self):
        self.sp_K.setValue(gl.CASTEDO_DEFECTO["K"])
        self.sp_pmax.setValue(gl.CASTEDO_DEFECTO["pendiente_max"])
        self.sp_fdesl.setValue(100 * gl.CASTEDO_DEFECTO["derrubios_deslizamiento"])
        self.sp_fvuel.setValue(100 * gl.CASTEDO_DEFECTO["derrubios_vuelco"])

    def _estabilidad_original(self):
        for k, ck in self.ck_est.items():
            ck.setChecked(k == "fellenius")
        self.sp_rebanadas.setValue(5)
        self.cb_criterio.setCurrentIndex(0)

    def _cambio_estabilidad(self, *_):
        if not all(hasattr(self, x) for x in ("ck_est", "sp_rebanadas", "cb_criterio",
                                               "lb_aviso_estab")):
            return   # todavia se esta construyendo la ventana
        ero, est = self.metodos_elegidos()
        n = len(ero) * len(est)
        if not ero or not est:
            self.lb_aviso_estab.setText("⚠ Marca al menos un método de erosión y uno de "
                                        "estabilidad.")
            self.lb_aviso_estab.setStyleSheet("color: #b0451d;")
            return
        partes = []
        if n > 1:
            partes.append(f"Se harán {n} cálculos, uno por combinación; se comparan en la "
                          f"pestaña «Comparación».")
        original = (ero == ["trenhaile"] and est == ["fellenius"]
                    and self.sp_rebanadas.value() == 5
                    and self.cb_criterio.currentData() == "todas")
        if original and self.modo == "lago":
            partes.append("Como el modelo original: resultados idénticos al Fortran.")
        if "castedo" in ero:
            partes.append("Castedo usa la resistencia a compresión (UCS) de los estratos y no "
                          "usa los parámetros de Trenhaile (umbrales, run-up, playa).")
        if "spencer" in est:
            partes.append("Spencer tarda unas 1,5 veces más que Bishop.")
        if not original and self.modo == "lago":
            partes.append("⚠ Con métodos distintos del original los resultados no coinciden "
                          "con el Fortran.")
        self.lb_aviso_estab.setText(" ".join(partes))
        self.lb_aviso_estab.setStyleSheet("color: #b0451d;" if "⚠" in " ".join(partes) else "")

    def opciones_mar(self):
        """Opciones del modo mar para gl.run (None en modo lago)."""
        if self.modo != "mar":
            return None
        return {"amplitud": self.sp_amp.value(), "subida_mm_ano": self.sp_subida.value(),
                "niveles": self.sp_nmarea.value()}

    def _mar_defecto(self):
        self.sp_amp.setValue(gl.MAR_TESIS["amplitud"])
        self.sp_subida.setValue(gl.MAR_TESIS["subida_mm_ano"])
        self.sp_nmarea.setValue(gl.MAR_TESIS["niveles"])

    def cambiar_modo(self, modo, limpiar=True):
        """Cambia entre lago (modelo original) y mar (adaptacion con marea)."""
        if modo == self.modo and not limpiar:
            return
        cambia = modo != self.modo
        self.modo = modo
        self.acciones_modo[modo].setChecked(True)
        if cambia:
            # metodo de erosion por defecto de cada modo
            self.ck_ero["trenhaile"].setChecked(modo == "lago")
            self.ck_ero["castedo"].setChecked(modo == "mar")
        self.b_temporales.setVisible(modo == "mar")
        self.marca.setText("LAGO" if modo == "lago" else "MAR")
        self.lb_nivel_barra.setText("Nivel del lago:" if modo == "lago" else "Nivel medio del mar:")
        self.grp_nivel_lago.setVisible(modo == "lago")
        self.grp_marea.setVisible(modo == "mar")
        self.lb_leyenda_linea.setText(self._texto_leyenda())
        self._actualizar_titulo()
        self._cambio_estabilidad()
        self.tb_salidas.setHorizontalHeaderLabels(
            ["Ver", "Mes", "Año",
             "Retroceso\nnivel lago (m)" if modo == "lago" else "Retroceso nivel\nmedio mar (m)",
             "Retroceso\ncresta (m)"])
        self.mostrar_oleaje()
        if limpiar:
            self._limpiar_resultados()
        self.dibujar_perfil()
        self._cambio_opciones_nivel()

    def _texto_leyenda(self):
        inv = "invierno: lago helado" if getattr(self, "modo", "lago") == "lago" else "sin oleaje"
        return (f'<span style="color:{COLOR_ROTURA}">┃</span> deslizamiento &nbsp; '
                f'<span style="color:{COLOR_VUELCO}">╻</span> vuelco &nbsp;&nbsp; '
                f'<span style="color:{COLOR_MES_CALCULADO}; background:{COLOR_MES_CALCULADO}">'
                f'&nbsp;&nbsp;&nbsp;</span> mes con erosión &nbsp;&nbsp; '
                f'<span style="color:{COLOR_INVIERNO}">█</span> {inv} &nbsp;&nbsp; '
                f'<span style="color:{COLOR_DERRUBIOS}">█</span> pie protegido por derrubios')

    def _actualizar_titulo(self):
        modo = "lago" if self.modo == "lago" else "mar"
        t = f"{APP_NAME} [{modo}]"
        if self.ruta:
            t += f" — {self.ruta}"
        self.setWindowTitle(t)

    def opciones_nivel(self):
        """Opciones de nivel variable para gl.run (None = constante, como el original)."""
        if not hasattr(self, "ck_estacional") or self.modo != "lago":
            return None
        est = self.ck_estacional.isChecked()
        inter = self.ck_interanual.isChecked()
        if not (est or inter):
            return None
        return {"estacional": ({k: sp.value() for k, sp in self.sp_estac.items()}
                               if est else None),
                "interanual": inter}

    def _nivel_tesis(self):
        for k, sp in self.sp_estac.items():
            sp.setValue(gl.NIVEL_ESTACIONAL_TESIS[k])
        self.ck_estacional.setChecked(True)
        self.ck_interanual.setChecked(True)

    def _cambio_opciones_nivel(self, *_):
        for sp in self.sp_estac.values():
            sp.setEnabled(self.ck_estacional.isChecked())
        self.dibujar_nivel()

    def dibujar_nivel(self):
        lz = self.graf_nivel
        ax = lz.ax
        ax.clear()
        lz.estilo()
        lago = self.sp_lake.value() if hasattr(self, "sp_lake") else 0
        hay = False
        res = self.res
        if res is not None and len(res["mensual"]["meses"]):
            men = res["mensual"]
            es_mar = res.get("modo") == "mar"
            ax.step(men["meses"] / 12.0, men["nivel"], where="post", color=ACCENT, lw=1.3,
                    label="Nivel medio del mar" if es_mar else "Nivel usado en el cálculo")
            if es_mar and res.get("mar") and res["mar"]["amplitud"] > 0:
                a_ = res["mar"]["amplitud"]
                ax.fill_between(men["meses"] / 12.0, men["nivel"] - a_, men["nivel"] + a_,
                                step="post", color="#4fa3e0", alpha=0.15,
                                label=f"Rango de marea (±{a_:.1f} m)")
            hay = True
        else:
            ax.axhline(lago, color=ACCENT, ls="--", lw=1, label="Nivel del lago (constante)")
        if self.serie_nivel is not None:
            nombre, t, z = self.serie_nivel
            ax.plot(t, z, color="#1e88e5", lw=1.0, label=f"Serie cargada: {nombre}")
            hay = True
        ax.set_xlabel("Tiempo (años)" if self.serie_nivel is None else "Tiempo")
        ax.set_ylabel("Nivel (m)")
        ax.ticklabel_format(axis="y", useOffset=False)
        lo, hi = ax.get_ylim()
        if hi - lo < 0.2:          # nivel (casi) constante: margen visible
            c = 0.5 * (lo + hi)
            ax.set_ylim(c - 0.5, c + 0.5)
        ax.legend(fontsize=8)
        lz.fig.tight_layout()
        lz.draw_idle()
        if not hay:
            pass
        op = self.opciones_nivel() if hasattr(self, "ck_estacional") else None
        if self.modo == "mar":
            pass
        elif op is None:
            self.lb_aviso_nivel.setText(
                "Desactivado: nivel constante, igual que el Fortran original (resultados "
                "idénticos).")
            self.lb_aviso_nivel.setStyleSheet("")
        else:
            self.lb_aviso_nivel.setText(
                "⚠ Activado: el original tiene esta variación comentada; con ella los "
                "resultados ya no coinciden con el Fortran. Los desniveles se redondean a "
                "múltiplos de Δz.")
            self.lb_aviso_nivel.setStyleSheet("color: #b0451d;")
        usado = ""
        if res is not None:
            nl = res.get("nivel_lago")
            usado = (" Último cálculo: nivel " +
                     ("constante." if nl is None else
                      "variable (" + ", ".join(x for x, ok in
                                               (("estacional", nl.get("estacional")),
                                                ("interanual", nl.get("interanual"))) if ok)
                      + ")."))
        if self.modo == "mar":
            if res is not None and res.get("modo") == "mar":
                usado = (f" Último cálculo: marea ±{res['mar']['amplitud']:.2f} m en "
                         f"{len(res['mar']['bandas'])} niveles, subida "
                         f"{res['mar']['subida_mm_ano']:g} mm/año.")
            else:
                usado = ""
            self.lb_nivel.setText(
                "La serie cargada (por ejemplo, marea.dat) solo se visualiza. En el cálculo "
                "interviene la marea sinusoidal definida arriba y la subida del nivel medio "
                "del mar." + usado)
        else:
            self.lb_nivel.setText(
                "La serie de marea.dat solo se visualiza: el modelo no la lee (en el Fortran "
                "original esa lectura está comentada). La variación que sí interviene en el "
                "cálculo es la de arriba: estacional (Tabla 10.iii de la tesis: +0,40 m en "
                "primavera y +0,30 m en verano) y ciclo interanual." + usado)

    # ------------------------------------------------------ edicion tablas
    def anadir_punto(self):
        fila = self.tb_perfil.currentRow()
        fila = self.tb_perfil.rowCount() if fila < 0 else fila + 1
        self._fila_perfil(0.0, 0.0, 0, fila)
        self._cambio_datos()

    def eliminar_punto(self):
        for r in sorted({i.row() for i in self.tb_perfil.selectedIndexes()}, reverse=True):
            self.tb_perfil.removeRow(r)
        self._cambio_datos()

    def ordenar_perfil(self):
        try:
            d = self.leer_datos()
        except ValueError as exc:
            QMessageBox.warning(self, APP_NAME, f"Dato no válido: {exc}")
            return
        pts = sorted(d["Pzy"], key=lambda p: (-p[0], p[1]))
        self._bloqueo = True
        self.tb_perfil.setRowCount(0)
        for z, y, e in pts:
            self._fila_perfil(float(z), float(y), int(e))
        self._bloqueo = False
        self._cambio_datos()

    def anadir_estrato(self):
        self._fila_estrato([200, 8, 23, 25, 23.4, 10])
        self._cambio_datos()

    def eliminar_estrato(self):
        for r in sorted({i.row() for i in self.tb_estratos.selectedIndexes()}, reverse=True):
            self.tb_estratos.removeRow(r)
        self._cambio_datos()

    # ------------------------------------------------------------ ficheros
    def _carpeta(self):
        return os.path.dirname(self.ruta) if self.ruta else os.getcwd()

    def nuevo(self):
        self.ruta = None
        self.res = None
        self._actualizar_titulo()
        self.cargar_datos(caso_por_defecto())
        self._limpiar_resultados()

    def dialogo_abrir(self):
        ruta, _ = QFileDialog.getOpenFileName(self, "Abrir data.dat", self._carpeta(),
                                              "Datos (*.dat);;Todos (*)")
        if ruta:
            self.abrir(ruta)

    def abrir(self, ruta):
        try:
            d = gl.read_data(ruta)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, APP_NAME, f"No se ha podido leer {ruta}:\n{exc}")
            return
        self.ruta = os.path.abspath(ruta)
        self.res = None
        self.cargar_datos(d)
        extra = ""
        ro = self._ruta_opciones(self.ruta)
        if os.path.exists(ro):
            try:
                with open(ro, encoding="utf-8") as fh:
                    self.aplicar_opciones_herramienta(json.load(fh))
                extra = f" (con las opciones de {os.path.basename(ro)}: modo {self.modo})"
            except Exception as exc:  # noqa: BLE001
                extra = f" (no se han podido leer las opciones: {exc})"
        self._actualizar_titulo()
        self._limpiar_resultados()
        self.statusBar().showMessage(f"Abierto {self.ruta}{extra}")

    def guardar(self):
        if not self.ruta:
            return self.guardar_como()
        return self._guardar_en(self.ruta)

    def guardar_como(self):
        ruta, _ = QFileDialog.getSaveFileName(self, "Guardar data.dat",
                                              os.path.join(self._carpeta(), "data.dat"),
                                              "Datos (*.dat)")
        if not ruta:
            return False
        if self._guardar_en(ruta):
            self.ruta = os.path.abspath(ruta)
            self._actualizar_titulo()
            return True
        return False

    @staticmethod
    def _ruta_opciones(ruta):
        return os.path.splitext(ruta)[0] + ".herramienta.json"

    def opciones_herramienta(self):
        """Opciones de la herramienta que no caben en data.dat."""
        return {
            "modo": self.modo,
            "mar": {"amplitud": self.sp_amp.value(), "subida_mm_ano": self.sp_subida.value(),
                    "niveles": self.sp_nmarea.value()},
            "estabilidad": {"metodos": self.metodos_elegidos()[1],
                            "rebanadas": self.sp_rebanadas.value(),
                            "criterio": self.cb_criterio.currentData()},
            "erosion": {"metodos": self.metodos_elegidos()[0], "K": self.sp_K.value(),
                        "pendiente_max": self.sp_pmax.value(),
                        "derrubios_deslizamiento": self.sp_fdesl.value() / 100.0,
                        "derrubios_vuelco": self.sp_fvuel.value() / 100.0},
            "nivel_lago": {"estacional": self.ck_estacional.isChecked(),
                           "interanual": self.ck_interanual.isChecked(),
                           "subidas": {k: sp.value() for k, sp in self.sp_estac.items()}},
            "oleaje": self.seasons_por_modo,
        }

    def aplicar_opciones_herramienta(self, op):
        mar = op.get("mar", {})
        self.sp_amp.setValue(mar.get("amplitud", gl.MAR_TESIS["amplitud"]))
        self.sp_subida.setValue(mar.get("subida_mm_ano", gl.MAR_TESIS["subida_mm_ano"]))
        self.sp_nmarea.setValue(mar.get("niveles", gl.MAR_TESIS["niveles"]))
        est = op.get("estabilidad", {})
        metodos = est.get("metodos", [est.get("metodo", "fellenius")])
        for k, ck in self.ck_est.items():
            ck.setChecked(k in metodos)
        self.sp_rebanadas.setValue(int(est.get("rebanadas", 5)))
        i = self.cb_criterio.findData(est.get("criterio", "todas"))
        self.cb_criterio.setCurrentIndex(max(0, i))
        nl = op.get("nivel_lago", {})
        for k, v in nl.get("subidas", {}).items():
            if k in self.sp_estac:
                self.sp_estac[k].setValue(v)
        self.ck_estacional.setChecked(bool(nl.get("estacional")))
        self.ck_interanual.setChecked(bool(nl.get("interanual")))
        ole = op.get("oleaje")
        if ole:
            for m in ("lago", "mar"):
                if m in ole:
                    self.seasons_por_modo[m] = {k: tuple(list(x) for x in v)
                                                for k, v in ole[m].items()}
        self.cambiar_modo(op.get("modo", "lago"))
        ero = op.get("erosion", {})
        if ero:
            for k, ck in self.ck_ero.items():
                ck.setChecked(k in ero.get("metodos", ["trenhaile"]))
            self.sp_K.setValue(ero.get("K", gl.CASTEDO_DEFECTO["K"]))
            self.sp_pmax.setValue(ero.get("pendiente_max", gl.CASTEDO_DEFECTO["pendiente_max"]))
            self.sp_fdesl.setValue(100 * ero.get("derrubios_deslizamiento",
                                                 gl.CASTEDO_DEFECTO["derrubios_deslizamiento"]))
            self.sp_fvuel.setValue(100 * ero.get("derrubios_vuelco",
                                                 gl.CASTEDO_DEFECTO["derrubios_vuelco"]))

    def _guardar_en(self, ruta):
        try:
            gl.write_data(ruta, self.leer_datos())
            with open(self._ruta_opciones(ruta), "w", encoding="utf-8") as fh:
                json.dump(self.opciones_herramienta(), fh, indent=2, ensure_ascii=False)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, APP_NAME, f"No se ha podido guardar:\n{exc}")
            return False
        self.statusBar().showMessage(f"Guardado {ruta} (y opciones de la herramienta en "
                                     f"{os.path.basename(self._ruta_opciones(ruta))})")
        return True

    def guardar_figura(self, fig):
        ruta, _ = QFileDialog.getSaveFileName(self, "Guardar figura",
                                              os.path.join(self._carpeta(), "figura.png"),
                                              "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if ruta:
            fig.savefig(ruta, dpi=200, facecolor=fig.get_facecolor())

    def exportar_csv(self):
        if self.res is None:
            QMessageBox.information(self, APP_NAME, "Primero hay que calcular.")
            return
        ruta, _ = QFileDialog.getSaveFileName(self, "Exportar CSV",
                                              os.path.join(self._carpeta(), "recesion.csv"),
                                              "CSV (*.csv)")
        if not ruta:
            return
        meses, lago, cresta, estudio, _ = self._series()
        with open(ruta, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=";")
            w.writerow(["mes", "anio", "retroceso_nivel_lago_m", "retroceso_cresta_m"]
                       + [f"retroceso_punto_estudio_{i + 1}_m" for i in range(len(estudio))])
            for i, m in enumerate(meses):
                w.writerow([f"{m:.0f}", f"{m / 12:.4f}", f"{lago[i]:.6f}", f"{cresta[i]:.6f}"]
                           + [f"{s[i]:.6f}" for s in estudio])
        self.statusBar().showMessage(f"Exportado {ruta}")

    def cargar_resultados_dialogo(self):
        QMessageBox.information(
            self, APP_NAME,
            "Los resultados de cada cálculo se guardan en la carpeta «resultados» junto al "
            "data.dat (Matriz_dibujo.out, Matriz_roturas.out, …).\n"
            + (f"Último cálculo: {self.res_dir}" if self.res_dir else "Aún no hay cálculos."))

    def _limpiar_resultados(self):
        self.res = None
        self._eventos = []
        self.resultados = []
        if hasattr(self, "cb_resultado"):
            self.cb_resultado.blockSignals(True)
            self.cb_resultado.clear()
            self.cb_resultado.blockSignals(False)
            self.fila_resultado.setVisible(False)
            self.tb_comp.setRowCount(0)
            self.dibujar_comparacion()
        self._activar_linea(False)
        self.tb_salidas.setRowCount(0)
        self.tb_roturas.setRowCount(0)
        self.tabs_der.setTabText(2, "Roturas")
        self.lb_resdir.setText("Sin resultados. Pulsa «Calcular».")
        self.dibujar_perfil()
        self.dibujar_grafica()

    # -------------------------------------------------------------- calculo
    def calcular(self):
        if self.worker is not None:
            self.worker.cancelar()
            self.b_calc.setEnabled(False)
            return
        try:
            d = self.leer_datos()
        except ValueError as exc:
            QMessageBox.warning(self, APP_NAME, f"Hay un dato no válido en las tablas: {exc}")
            return
        errores = self.validar(d)
        ero, est = self.metodos_elegidos()
        if not ero or not est:
            errores.append("Marca al menos un método de erosión y uno de estabilidad "
                           "(pestaña «Parámetros»).")
        if "castedo" in ero and any(float(x) <= 0 for x in d["sigmac"][1:]):
            errores.append("La erosión de Castedo necesita la resistencia a compresión (UCS) "
                           "de todos los estratos.")
        if errores:
            QMessageBox.warning(self, APP_NAME, "Revisa los datos:\n\n• " + "\n• ".join(errores))
            return
        if self.ruta:
            self.res_dir = os.path.join(os.path.dirname(self.ruta), "resultados")
        else:
            self.res_dir = os.path.join(os.path.expanduser("~"), "GreatLake", "resultados")
        os.makedirs(self.res_dir, exist_ok=True)
        seasons = (None if self.seasons == self._seasons_defecto()
                   else copy.deepcopy(self.seasons))
        varios = len(ero) * len(est) > 1
        configs = []
        for e in ero:
            for m in est:
                nombre = f"{NOMBRES_EROSION[e]} + {NOMBRES_ESTABILIDAD[m]}"
                carpeta = os.path.join(self.res_dir, f"{e}_{m}") if varios else self.res_dir
                os.makedirs(carpeta, exist_ok=True)
                gl.write_data(os.path.join(carpeta, "data.dat"), d)
                configs.append({"nombre": nombre, "clave_erosion": e, "clave_estabilidad": m,
                                "dir": carpeta, "seasons": seasons,
                                "nivel_lago": self.opciones_nivel(),
                                "estabilidad": self.opciones_estabilidad(m),
                                "mar": self.opciones_mar(),
                                "erosion": self.opciones_erosion(e)})

        self.hilo = QThread()
        self.worker = Worker(configs)
        self.worker.moveToThread(self.hilo)
        self.hilo.started.connect(self.worker.run)
        self.worker.progreso.connect(self._progreso)
        self.worker.terminado.connect(self._terminado)
        self.worker.error.connect(self._error)
        self.progreso.setVisible(True)
        self.progreso.setValue(0)
        self.b_calc.setText("Detener")
        self.statusBar().showMessage(f"Calculando {len(configs)} combinación(es)…")
        self.hilo.start()

    def _progreso(self, i, n, t, tsim):
        self.progreso.setMaximum(n * tsim)
        self.progreso.setValue((i - 1) * tsim + min(t, tsim))
        self.progreso.setFormat((f"Cálculo {i}/{n} · " if n > 1 else "") +
                                f"mes {min(t, tsim)} / {tsim}")

    def _fin_hilo(self):
        self.hilo.quit()
        self.hilo.wait()
        self.hilo = None
        self.worker = None
        self.b_calc.setText("Calcular")
        self.b_calc.setEnabled(True)
        self.progreso.setVisible(False)

    def _terminado(self, resultados, segundos):
        self._fin_hilo()
        self.resultados = resultados
        self.cb_resultado.blockSignals(True)
        self.cb_resultado.clear()
        for r in resultados:
            self.cb_resultado.addItem(r["nombre"])
        self.cb_resultado.blockSignals(False)
        self.fila_resultado.setVisible(len(resultados) > 1)
        self.res = resultados[0]["res"]
        self.res_dir_actual = resultados[0]["dir"]
        self.lb_resdir.setText(f"Resultados guardados en: {self.res_dir}")
        self.mostrar_resultados()
        self.mostrar_comparacion()
        self.tabs_datos.setCurrentIndex(3)
        if len(resultados) > 1:
            self.tabs_der.setCurrentIndex(1)
        est = self.res.get("estabilidad", {})
        txt_est = (f"{est.get('rebanadas')} rebanadas, criterio «{est.get('criterio')}»")
        if len(resultados) == 1:
            r = resultados[0]
            self.statusBar().showMessage(
                f"Cálculo terminado en {segundos:.1f} s · modo {self.res.get('modo', 'lago')} · "
                f"{r['nombre']} ({txt_est}) · {len(self.res['roturas'])} deslizamientos y "
                f"{len(self.res.get('vuelcos', []))} vuelcos · salidas en {self.res_dir}")
        else:
            self.statusBar().showMessage(
                f"{len(resultados)} cálculos terminados en {segundos:.1f} s · modo "
                f"{self.res.get('modo', 'lago')} · {txt_est} · compáralos en la pestaña "
                f"«Comparación» · salidas en {self.res_dir}")

    def _elegir_resultado(self, i):
        if 0 <= i < len(getattr(self, "resultados", [])):
            self.res = self.resultados[i]["res"]
            self.tb_roturas.clearSelection()
            self.mostrar_resultados()
            self.dibujar_nivel()

    def _fila_comparacion(self):
        filas = sorted({it.row() for it in self.tb_comp.selectedItems()})
        if filas and filas[0] != self.cb_resultado.currentIndex():
            self.cb_resultado.setCurrentIndex(filas[0])

    @staticmethod
    def _series_de(res):
        D = res["DibMatriz"]
        P0 = res["profile_inicial"]
        meses = D[1, 3:]
        val = meses > 0
        cols = np.arange(3, D.shape[1])[val]
        meses = meses[val]
        TL = res["meta"]["TL"]
        kc = res["meta"]["crestatalud"]
        lago = np.array([P0[2, TL] - D[TL + 1, c] for c in cols])
        cresta = np.array([P0[2, kc] - D[kc + 1, c] for c in cols])
        return meses, lago, cresta

    def mostrar_comparacion(self):
        t = self.tb_comp
        t.setRowCount(0)
        for i, r in enumerate(getattr(self, "resultados", [])):
            res = r["res"]
            meses, lago, cresta = self._series_de(res)
            anos = res["meta"]["tsim"] / 12.0
            vals = (NOMBRES_EROSION[r["erosion"]], NOMBRES_ESTABILIDAD[r["estabilidad"]],
                    str(len(res["roturas"])), str(len(res.get("vuelcos", []))),
                    _f2(cresta[-1] / anos) if len(cresta) else "",
                    _f2(lago[-1] / anos) if len(lago) else "", f"{r['segundos']:.1f}")
            t.insertRow(i)
            for k, val in enumerate(vals):
                c = QTableWidgetItem(val)
                c.setTextAlignment((Qt.AlignLeft if k < 2 else Qt.AlignRight) | Qt.AlignVCenter)
                c.setForeground(QColor(matplotlib.colors.to_hex(self._color_calculo(i))))
                t.setItem(i, k, c)
        self.dibujar_comparacion()

    @staticmethod
    def _color_calculo(i):
        return matplotlib.colormaps["tab10"](i % 10)

    def dibujar_comparacion(self, *_):
        lz = self.graf_comp
        ax = lz.ax
        ax.clear()
        lz.estilo()
        resultados = getattr(self, "resultados", [])
        if not resultados:
            ax.text(0.5, 0.5, "Calcula para comparar (marca varios métodos en «Parámetros»)",
                    transform=ax.transAxes, ha="center", color="#888", fontsize=9)
            ax.set_axis_off()
            lz.draw_idle()
            return
        modo = self.cb_comp.currentIndex()
        for i, r in enumerate(resultados):
            res = r["res"]
            c = self._color_calculo(i)
            meses, lago, cresta = self._series_de(res)
            if modo == 0:
                ax.plot(meses / 12.0, cresta, "-o", ms=3, color=c, label=r["nombre"])
            elif modo == 1:
                ax.plot(meses / 12.0, lago, "-o", ms=3, color=c, label=r["nombre"])
            elif modo == 2:
                ev = sorted([(m, v) for m, v, *_ in res["roturas"]] +
                            [(m, v) for m, v, *_ in res.get("vuelcos", [])])
                if ev:
                    ax.step([0] + [m / 12.0 for m, _ in ev] + [res["meta"]["tsim"] / 12.0],
                            np.concatenate([[0], np.cumsum([v for _, v in ev]),
                                            [sum(v for _, v in ev)]]),
                            where="post", color=c, label=r["nombre"])
            else:
                N = res["meta"]["N"]
                P0 = res["profile_inicial"]
                if i == 0:
                    ax.plot(P0[2, 1:N + 1], P0[1, 1:N + 1], color="#999", lw=1, ls="--",
                            label="Perfil inicial")
                ax.plot(res["DibMatriz"][2:N + 2, -1], P0[1, 1:N + 1], color=c, lw=1.4,
                        label=r["nombre"])
        if modo == 3:
            r0 = resultados[0]["res"]
            P0 = r0["profile_inicial"]
            kc, kb = r0["meta"]["crestatalud"], max(1, r0["meta"]["TL"])
            H = abs(P0[1, kc] - P0[1, kb]) or 1.0
            ymin = min(float(np.min(r["res"]["DibMatriz"][kc + 1:kb + 2, -1]))
                       for r in resultados)
            ax.set_xlim(ymin - 0.4 * H, float(P0[2, kc]) + 2.5 * H)
            ax.set_ylim(float(P0[1, kb]) - 0.6 * H, float(P0[1, 1]) + 0.1 * H)
            ax.set_xlabel("Distancia horizontal y (m)")
            ax.set_ylabel("Cota z (m)")
        else:
            ax.set_xlabel("Tiempo (años)")
            ax.set_ylabel(["Retroceso de la cresta (m)", "Retroceso en el nivel del agua (m)",
                           "Volumen acumulado (m³/m)"][modo])
        ax.legend(fontsize=8)
        lz.fig.tight_layout()
        lz.draw_idle()

    def exportar_comparacion(self):
        if not getattr(self, "resultados", []):
            QMessageBox.information(self, APP_NAME, "Primero hay que calcular.")
            return
        ruta, _ = QFileDialog.getSaveFileName(self, "Exportar comparación",
                                              os.path.join(self._carpeta(), "comparacion.csv"),
                                              "CSV (*.csv)")
        if not ruta:
            return
        with open(ruta, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter=";")
            w.writerow(["erosion", "estabilidad", "deslizamientos", "vuelcos",
                        "retroceso_cresta_m_ano", "retroceso_nivel_agua_m_ano", "tiempo_s",
                        "carpeta"])
            for r in self.resultados:
                res = r["res"]
                _, lago, cresta = self._series_de(res)
                anos = res["meta"]["tsim"] / 12.0
                w.writerow([r["erosion"], r["estabilidad"], len(res["roturas"]),
                            len(res.get("vuelcos", [])),
                            f"{cresta[-1] / anos:.4f}" if len(cresta) else "",
                            f"{lago[-1] / anos:.4f}" if len(lago) else "",
                            f"{r['segundos']:.1f}", r["dir"]])
        self.statusBar().showMessage(f"Exportado {ruta}")

    def _error(self, msg):
        self._fin_hilo()
        self.statusBar().showMessage(msg.splitlines()[0])
        QMessageBox.warning(self, APP_NAME, msg)

    def acerca_de(self):
        QMessageBox.about(
            self, APP_NAME,
            "<b>Great Lake</b> — erosión y recesión de acantilados de tills en lagos sin "
            "marea (lago Erie).<br><br>Modelo: R. Castedo (UPM), Great_lake_v.1 (Fortran), "
            "traducido a Python con resultados idénticos.<br><br>"
            "<b>Puntos del perfil</b>: 1 punto de estudio · 2 cresta del talud · "
            "3 punto de rezume · 4 base del talud.<br>"
            "<b>Vista</b>: rueda/lupa de la barra de matplotlib para ampliar; "
            "«Acantilado» encuadra la zona activa.")

    def closeEvent(self, ev):
        if self.worker is not None:
            self.worker.cancelar()
            self.hilo.quit()
            self.hilo.wait(3000)
        super().closeEvent(ev)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setStyleSheet(STYLE)
    f = QFont()
    f.setPointSize(10)
    app.setFont(f)
    ruta = sys.argv[1] if len(sys.argv) > 1 else None
    w = Ventana(ruta)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
