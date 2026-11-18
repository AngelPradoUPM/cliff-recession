#!/usr/bin/env python3
"""
Great Lake - herramienta grafica para el modelo de recesion de acantilados.

Interfaz de escritorio (PySide6 + matplotlib) sobre ``great_lake.py``:

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
import os
import sys
import time

import numpy as np
from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtGui import QAction, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDockWidget,
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
import great_lake as gl  # noqa: E402

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
    d = dict(rn=1, Utsim=2, tsim=10, UTdib=2, Tdib=1, Dz=0.1, lakelevel=0,
             Qb=0.0, Gs=0.001, groundwater=0, profagua=1.0, ndatosresis=1,
             SFmin=5.0, Kss=2.4e-6, p=1.0, WFmin=70.0, MMM=1e-8,
             Rsc=5.0, Cb=1e-7, KRs=1.0, Bmat=1.0, Mf=2.9)
    for nm, v in (("sigmac", 200.0), ("Shearstrength", 8.0), ("Ecohesion", 23.0),
                  ("Ephi", 25.0), ("pesounit", 23.4), ("potestrato", 120.0)):
        d[nm] = np.array([0.0, v])
    d["Pzy"] = [[15.0000001, -500.0, 0.0], [15.0000001, -6.1, 3.0],
                [15.0000001, -6.0, 2.0], [0.0, 0.0, 4.0], [-60.0, 1500.0, 0.0]]
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
class Worker(QObject):
    progreso = Signal(int, int)
    terminado = Signal(object, float)
    error = Signal(str)

    def __init__(self, workdir, seasons):
        super().__init__()
        self.workdir = workdir
        self.seasons = seasons
        self._cancelar = False

    def cancelar(self):
        self._cancelar = True

    def run(self):
        t0 = time.time()
        last = [0.0]

        def prog(t, tsim):
            now = time.time()
            if now - last[0] > 0.05 or t >= tsim:
                last[0] = now
                self.progreso.emit(t, tsim)
            return not self._cancelar

        try:
            res = gl.run(self.workdir, verbose=False, progress=prog, seasons=self.seasons)
        except gl.Cancelled:
            self.error.emit("Cálculo detenido por el usuario.")
            return
        except Exception as exc:  # noqa: BLE001
            self.error.emit(f"El cálculo ha fallado:\n{type(exc).__name__}: {exc}")
            return
        self.terminado.emit(res, time.time() - t0)


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


# ---------------------------------------------------------------------------
# Ventana principal
# ---------------------------------------------------------------------------
class Ventana(QMainWindow):
    def __init__(self, ruta=None):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1600, 950)
        self.ruta = None
        self.res = None
        self.res_dir = None
        self.serie_nivel = None
        self.seasons = copy.deepcopy(gl.SEASONS)
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

    # ------------------------------------------------------------------ UI
    def _crear_barra(self):
        barra = QWidget(objectName="topbar")
        lay = QHBoxLayout(barra)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(8)

        marca = QToolButton(objectName="brand")
        marca.setText("GREAT LAKE")
        marca.setPopupMode(QToolButton.InstantPopup)
        menu = QMenu(marca)
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

        lay.addWidget(QLabel("Nivel del lago:"))
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
                f.addRow(etiqueta, widget)
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
        nav = NavigationToolbar2QT(self.lienzo, w)
        nav.setStyleSheet("background: #e9e9e9; border: none;")
        v.addWidget(nav)
        v.addWidget(self.lienzo, 1)
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
                                  "Volumen de roturas acumulado"])
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

        # --- Roturas
        w = QWidget()
        v = QVBoxLayout(w)
        self.tb_roturas = QTableWidget(0, 3)
        self.tb_roturas.setHorizontalHeaderLabels(["Mes", "Año", "Volumen (m³/m)"])
        self.tb_roturas.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tb_roturas.verticalHeader().setVisible(False)
        self.tb_roturas.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tb_roturas.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tb_roturas.itemSelectionChanged.connect(self.dibujar_perfil)
        v.addWidget(self.tb_roturas)
        nota = QLabel("Selecciona una o varias roturas para ver el perfil antes (discontinuo) "
                      "y después en la vista central.", objectName="note")
        nota.setWordWrap(True)
        v.addWidget(nota)
        b = QPushButton("Quitar selección")
        b.clicked.connect(self.tb_roturas.clearSelection)
        v.addWidget(b)
        tabs.addTab(w, "Roturas")

        # --- Oleaje
        w = QWidget()
        v = QVBoxLayout(w)
        fila = QHBoxLayout()
        fila.addWidget(QLabel("Estación:"))
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
        b = QPushButton("Cargar serie (marea.dat)…")
        b.clicked.connect(self.cargar_serie_nivel)
        fila.addWidget(b)
        b = QPushButton("Quitar")
        b.clicked.connect(self.quitar_serie_nivel)
        fila.addWidget(b)
        fila.addStretch(1)
        v.addLayout(fila)
        self.graf_nivel = Lienzo()
        v.addWidget(self.graf_nivel, 1)
        self.lb_nivel = QLabel("", objectName="note")
        self.lb_nivel.setWordWrap(True)
        v.addWidget(self.lb_nivel)
        tabs.addTab(w, "Nivel / mareas")

        self.tabs_der = tabs
        dock.setWidget(tabs)
        dock.setMinimumWidth(560)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.mostrar_oleaje()
        self.dibujar_nivel()

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
            if d["lakelevel"] - fondo < Lmax / 2 + 1:
                errores.append(
                    f"El perfil debe llegar al menos a {Lmax / 2 + 1:.1f} m por debajo del nivel "
                    f"del lago (media longitud de onda de la ola más larga); ahora llega a "
                    f"{d['lakelevel'] - fondo:.1f} m.")
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
        terreno = y0
        if self.res is not None and self.res["meta"]["N"] == N:
            vis = [r for r in range(self.tb_salidas.rowCount())
                   if self.tb_salidas.item(r, 0) is not None
                   and self.tb_salidas.item(r, 0).checkState() == Qt.Checked]
            if vis:
                terreno = self.res["DibMatriz"][2:N + 2, vis[-1] + 3]
        ax.fill_betweenx(z, ymin - 1e6, terreno, color="#6b5a45", alpha=0.55, lw=0)
        ax.axhline(zlago, color="#4fa3e0", lw=1.2, ls="-", label=f"Nivel del lago ({zlago:.2f} m)")
        ax.plot(y0, z, color="#f0f0f0", lw=1.8, label="Perfil inicial")
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
            tmax = max(1.0, float(np.max(Dib[1, 3:])) / 12.0)
            norm = matplotlib.colors.Normalize(0.0, tmax)
            for r in cols:
                k = r + 3
                ano = Dib[1, k] / 12.0
                if ano <= 0:
                    continue
                ax.plot(Dib[2:N + 2, k], z, color=cmap(norm(ano)), lw=1.2,
                        alpha=0.3 if sel else 1.0)
            if cols:
                sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
                cb = lz.fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.04)
                cb.set_label("Perfil en el año", color=CANVAS_FG, fontsize=8)
                cb.ax.tick_params(colors=CANVAS_FG, labelsize=7)
                cb.outline.set_edgecolor("#5a6270")
            for r in sel:
                mes, vol, antes, despues = self.res["roturas"][r]
                ax.plot(antes, z, color="#ff6b6b", lw=1.1, ls="--",
                        label=f"Antes rotura mes {mes}")
                ax.plot(despues, z, color="#ff3030", lw=1.6,
                        label=f"Después rotura mes {mes} ({vol:.1f} m³/m)")

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
            if hay_res:
                ys.append(float(np.min(self.res["DibMatriz"][kc + 1:kt + 2, 3:])))
            x0 = min(ys) - 0.5 * H
            x1 = float(np.max(P[2, tramo])) + 1.0 * H
            ax.set_xlim(x0, x1)
            ax.set_ylim(zlago - 0.5 * H, zmax + 0.15 * H)
        else:
            ax.set_xlim(ymin, ymax)
            ax.set_ylim(zmin, zmax + 0.05 * (zmax - zmin))
        if self.ck_leyenda.isChecked():
            leg = ax.legend(loc="best", fontsize=7.5, facecolor="#2f3540",
                            edgecolor="#5a6270", labelcolor=CANVAS_FG)
            leg.set_draggable(True)
        lz.fig.tight_layout()
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

        tr = self.tb_roturas
        tr.setRowCount(0)
        for i, (mes, vol, _, _) in enumerate(r["roturas"]):
            tr.insertRow(i)
            for k, val in enumerate((f"{mes}", f"{mes / 12:.1f}", f"{vol:.3f}")):
                c = QTableWidgetItem(val)
                c.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                tr.setItem(i, k, c)
        self.tabs_der.setTabText(1, f"Roturas ({len(r['roturas'])})")
        self.dibujar_perfil()
        self.dibujar_grafica()

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
            rot = self.res["roturas"]
            if rot:
                t = np.array([m for m, *_ in rot]) / 12.0
                v = np.cumsum([vv for _, vv, *_ in rot])
                ax.step(t, v, where="post", color=ACCENT)
                ax.plot(t, v, "o", color=ACCENT, ms=3)
            else:
                ax.text(0.5, 0.5, "No se han producido roturas circulares",
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

    def restablecer_oleaje(self):
        self.seasons = copy.deepcopy(gl.SEASONS)
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
            ax.text(0.5, 0.5, "Sin oleaje (el invierno no erosiona)", transform=ax.transAxes,
                    ha="center", color="#888")
        lz.fig.tight_layout()
        lz.draw_idle()
        modificado = self.seasons != gl.SEASONS
        self.lb_oleaje.setText(
            "Oleaje por estación usado por el modelo (Trenhaile, 2009). "
            "Primavera: meses 4-6, verano: 7-9, otoño: 10-12; los meses de invierno (1-3) "
            "no se calculan."
            + ("\n⚠ Valores modificados respecto al código original: los resultados ya no "
               "coincidirán con el Fortran." if modificado else ""))
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

    def dibujar_nivel(self):
        lz = self.graf_nivel
        ax = lz.ax
        ax.clear()
        lz.estilo()
        lago = self.sp_lake.value() if hasattr(self, "sp_lake") else 0
        if self.serie_nivel is not None:
            nombre, t, z = self.serie_nivel
            ax.plot(t, z, color="#1e88e5", lw=1.2, label=nombre)
            ax.axhline(lago, color=ACCENT, ls="--", lw=1, label="Nivel usado por el modelo")
            ax.set_xlabel("Tiempo")
            ax.legend(fontsize=8)
        else:
            ax.axhline(lago, color=ACCENT, ls="--", lw=1)
            ax.text(0.5, 0.6, "Sin serie cargada", transform=ax.transAxes, ha="center",
                    color="#888")
        ax.set_ylabel("Nivel (m)")
        lz.fig.tight_layout()
        lz.draw_idle()
        self.lb_nivel.setText(
            "El modelo usa un nivel del lago constante («Nivel del lago» en la barra "
            "superior). En el código Fortran original la lectura de marea.dat y la variación "
            "interanual del nivel están desactivadas, así que la serie cargada aquí solo se "
            "visualiza: no interviene en el cálculo.")

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
        self.setWindowTitle(APP_NAME)
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
        self.setWindowTitle(f"{APP_NAME} — {self.ruta}")
        self.cargar_datos(d)
        self._limpiar_resultados()
        self.statusBar().showMessage(f"Abierto {self.ruta}")

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
            self.setWindowTitle(f"{APP_NAME} — {self.ruta}")
            return True
        return False

    def _guardar_en(self, ruta):
        try:
            gl.write_data(ruta, self.leer_datos())
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, APP_NAME, f"No se ha podido guardar:\n{exc}")
            return False
        self.statusBar().showMessage(f"Guardado {ruta}")
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
        self.tb_salidas.setRowCount(0)
        self.tb_roturas.setRowCount(0)
        self.tabs_der.setTabText(1, "Roturas")
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
        if errores:
            QMessageBox.warning(self, APP_NAME, "Revisa los datos:\n\n• " + "\n• ".join(errores))
            return
        if self.ruta:
            self.res_dir = os.path.join(os.path.dirname(self.ruta), "resultados")
        else:
            self.res_dir = os.path.join(os.path.expanduser("~"), "GreatLake", "resultados")
        os.makedirs(self.res_dir, exist_ok=True)
        gl.write_data(os.path.join(self.res_dir, "data.dat"), d)

        self.hilo = QThread()
        seasons = None if self.seasons == gl.SEASONS else copy.deepcopy(self.seasons)
        self.worker = Worker(self.res_dir, seasons)
        self.worker.moveToThread(self.hilo)
        self.hilo.started.connect(self.worker.run)
        self.worker.progreso.connect(self._progreso)
        self.worker.terminado.connect(self._terminado)
        self.worker.error.connect(self._error)
        self.progreso.setVisible(True)
        self.progreso.setValue(0)
        self.b_calc.setText("Detener")
        self.statusBar().showMessage("Calculando…")
        self.hilo.start()

    def _progreso(self, t, tsim):
        self.progreso.setMaximum(tsim)
        self.progreso.setValue(min(t, tsim))
        self.progreso.setFormat(f"Mes {min(t, tsim)} / {tsim}")

    def _fin_hilo(self):
        self.hilo.quit()
        self.hilo.wait()
        self.hilo = None
        self.worker = None
        self.b_calc.setText("Calcular")
        self.b_calc.setEnabled(True)
        self.progreso.setVisible(False)

    def _terminado(self, res, segundos):
        self._fin_hilo()
        self.res = res
        self.lb_resdir.setText(f"Resultados guardados en: {self.res_dir}")
        self.mostrar_resultados()
        self.tabs_datos.setCurrentIndex(3)
        self.statusBar().showMessage(
            f"Cálculo terminado en {segundos:.1f} s · {len(res['roturas'])} roturas · "
            f"salidas en {self.res_dir}")

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
