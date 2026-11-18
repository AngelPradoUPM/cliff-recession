"""Textos de ayuda de los parametros de la herramienta Great Lake.

Cada entrada es (explicacion sencilla, explicacion detallada o None). La
detallada admite HTML basico. Los valores de referencia proceden del capitulo
10 de la tesis (R. Castedo) y del codigo Great_lake_v.1.
"""

AYUDA = {
    # ------------------------------------------------------------ simulacion
    "rn": (
        "Un número para identificar la simulación.",
        "Se escribe en la primera fila de los ficheros de salida "
        "(<i>Matriz_dibujo.out</i>, <i>Matriz_roturas.out</i>…). No influye en el cálculo.",
    ),
    # --------------------------------------------------- esfuerzo cortante
    "SFmin": (
        "Cuánto tienen que «arrastrar» las olas sobre el fondo para empezar a erosionar "
        "la arcilla. Más alto = fondo más resistente.",
        "Es el esfuerzo cortante crítico τ<sub>c</sub> (Pa). Bajo el nivel del lago, en "
        "cada punto del perfil y para cada tipo de ola, el modelo calcula la velocidad "
        "orbital en el fondo U y el esfuerzo τ = ½·f<sub>w</sub>·ρ·U² "
        "(f<sub>w</sub> = 0,0034). La erosión es "
        "E = N<sub>o</sub>·K<sub>ss</sub>·(τ − τ<sub>c</sub>)<sup>p</sup>, solo si τ &gt; τ<sub>c</sub>."
        "<br><br>Tesis: 5 Pa (arcillas blandas) y 10 Pa (más resistentes); rango habitual "
        "0,5–20 Pa.",
    ),
    "Kss": (
        "Cuánto fondo se erosiona por cada unidad de «arrastre» por encima del umbral.",
        "Coeficiente K<sub>ss</sub> (m³/N) de la ecuación de erosión por esfuerzo "
        "cortante. Tesis: 2,4·10⁻⁶ m³/N (justificado en Trenhaile, 2009).",
    ),
    "p": (
        "Exponente de calibración de la erosión del fondo. Normalmente 1 (relación lineal).",
        "Exponente de (τ − τ<sub>c</sub>)<sup>p</sup>. Con p = 1 la erosión es proporcional "
        "al exceso de esfuerzo sobre el umbral; con p &gt; 1 crece más deprisa en temporales.",
    ),
    # ------------------------------------------------------ impacto del oleaje
    "WFmin": (
        "Fuerza mínima del golpe de las olas rotas para arrancar arcilla en el nivel del "
        "lago. Más alto = acantilado más resistente al oleaje.",
        "Umbral SF<sub>cr</sub> (kg/m²) de la erosión por impacto. La fuerza del oleaje "
        "roto es SF = 512,61·(H<sub>b</sub>/0,78)·e<sup>−0,05·S<sub>w</sub></sup>, con "
        "H<sub>b</sub> la altura de rotura y S<sub>w</sub> la anchura de la zona de "
        "rompientes (cuanto más ancha la plataforma, más se amortigua el golpe). La "
        "erosión, que se aplica en el punto del nivel del lago, es "
        "E = N<sub>o</sub>·MMM·(SF − SF<sub>cr</sub>).<br><br>Tesis: 70 (débil) y 175 "
        "(resistente) kg/m². Nota: la tesis da 0,01 como coeficiente de amortiguación y el "
        "código usa 0,05. Con playa, el impacto se anula si el nivel del lago queda dentro "
        "de la playa o si se forma una berma.",
    ),
    "MMM": (
        "Cuánto acantilado arranca cada golpe de ola por encima del umbral.",
        "Coeficiente K<sub>bf</sub> (m³/kg) de la erosión por impacto. Tesis: 1·10⁻⁸ m³/kg.",
    ),
    # --------------------------------------------------------------- run-up
    "Rsc": (
        "Umbral para que el agua que sube por la playa (run-up) erosione el pie del "
        "acantilado. Más alto = solo lo erosionan las olas grandes.",
        "Solo actúa si hay playa. El modelo calcula hasta dónde sube cada ola por encima "
        "de la playa (R<sub>s</sub>) y la erosión en el pie del acantilado es "
        "E = N<sub>o</sub>·C<sub>b</sub>·B<sub>mat</sub>·(K<sub>Rs</sub>·R<sub>s</sub> − "
        "R<sub>sc</sub>).<br><br>Tesis: 5 Pa (materiales débiles, erosionados por casi "
        "todas las olas) y 15 Pa (resistentes, solo las olas mayores).",
    ),
    "Cb": (
        "Cuánto erosiona el run-up por unidad de exceso sobre el umbral.",
        "Coeficiente K<sub>up</sub> (m³/N). Tesis: 1·10⁻⁷ m³/N, que da tasas de erosión por "
        "run-up de 0,10–0,8 m/año, coherentes con los Grandes Lagos.",
    ),
    "KRs": (
        "Convierte la distancia que sube el agua en un esfuerzo sobre el acantilado.",
        "Coeficiente m (N/m³) que multiplica la distancia de run-up R<sub>s</sub> para "
        "compararla con el umbral R<sub>sc</sub>. Tesis: 1 N/m³.",
    ),
    "Bmat": (
        "Factor del material de la playa en la erosión por run-up (1 = neutro).",
        "Multiplica la erosión por run-up. Permite representar playas cuyo material es "
        "más o menos abrasivo. El código original y la tesis usan 1.",
    ),
    # ---------------------------------------------------------------- playa
    "Qb": (
        "Cantidad de arena disponible para formar una playa al pie del acantilado "
        "(m³ por metro de costa). 0 = sin playa.",
        "Con suficiente arena, el modelo forma una playa desde el punto de rotura de las "
        "olas hasta el acantilado, con la pendiente de Sunamura (1989). Si no basta, mueve "
        "el pie de la playa hacia tierra de 0,1 en 0,1 m hasta que cabe o hasta llegar al "
        "acantilado. La playa protege del impacto directo, pero el run-up puede erosionar "
        "el pie.<br><br>Tesis: 0, 5 o 15 m³/m. Las tasas más altas se dan sin playa.",
    ),
    "Gs": (
        "Tamaño del grano de la arena de la playa. Grano más grueso = playa más empinada.",
        "Diámetro D (m) usado en la pendiente de la playa: "
        "tan β = 0,12 / [H<sub>b</sub>/(g<sup>½</sup>·D<sup>½</sup>·T)]<sup>½</sup> "
        "(Sunamura, 1989). Tesis: 0,001 m (arena) y 0,003 m (gránulo).",
    ),
    # ------------------------------------------------------ agua subterranea
    "groundwater": (
        "Cuánta agua hay dentro del acantilado. El agua lo hace menos estable.",
        "<b>0</b>: seco, sin presión de agua.<br>"
        "<b>1</b>: nivel freático que baja desde el punto de rezume del perfil hasta la "
        "profundidad indicada tierra adentro.<br>"
        "<b>2</b>: saturado, con el agua hasta la superficie.<br><br>"
        "En el cálculo de las roturas circulares, la presión del agua en la base de cada "
        "rebanada (u = γ<sub>w</sub>·h<sub>w</sub>) reduce la resistencia por rozamiento. "
        "Según la tesis, el retroceso es mayor en materiales saturados que en drenados.",
    ),
    "profagua": (
        "Profundidad del agua subterránea lejos del acantilado, medida desde la superficie "
        "(solo con el tipo 1).",
        "El nivel freático va del punto de rezume (marcado en el perfil) hasta esta "
        "profundidad en el extremo de tierra del perfil, con forma de curva. Conviene que "
        "el punto de rezume esté en la cara del acantilado y por debajo de la coronación.",
    ),
    # ------------------------------------------------------------- derrubios
    "Mf": (
        "Rapidez con la que el oleaje se lleva el material caído al pie. Mientras queda "
        "material, protege el acantilado.",
        "Coeficiente de movilidad M<sub>f</sub> (m²/día). Tras una rotura, el material "
        "caído tarda volumen/M<sub>f</sub> días en desaparecer. En ese tiempo el pie está "
        "protegido: el modelo salta los meses completos sin erosión y reduce el oleaje del "
        "mes parcial.<br><br>Tesis: 2 m²/día (movilidad baja) y 6 m²/día (alta).",
    ),
    # ------------------------------------------------------------ estabilidad
    "metodo": (
        "Cómo se calcula si un trozo del acantilado está a punto de deslizar (factor de "
        "seguridad: si es menor que 1, rompe).",
        "Los tres dividen la masa que desliza en rebanadas verticales:<br>"
        "<b>Fellenius</b> (original y tesis): ignora las fuerzas entre rebanadas. Es el más "
        "sencillo y el más conservador (da el FS más bajo).<br>"
        "<b>Bishop simplificado</b>: tiene en cuenta las fuerzas horizontales entre "
        "rebanadas. Es el estándar en la práctica para roturas circulares y suele dar un "
        "FS algo mayor (0–15 %).<br>"
        "<b>Spencer</b>: riguroso, cumple todas las ecuaciones de equilibrio. En roturas "
        "circulares queda muy cerca de Bishop (±3 %); es más lento.<br><br>"
        "Los parámetros de la tesis se calibraron con Fellenius.",
    ),
    "rebanadas": (
        "En cuántos trozos se divide la masa que desliza para calcular su estabilidad.",
        "El original usa 5. Con más rebanadas el cálculo es más preciso; a partir de "
        "20–30 apenas cambia. En los ejemplos, pasar de 5 a 30 cambia las tasas de "
        "retroceso menos de un 10 %.",
    ),
    "criterio": (
        "Qué deslizamiento se produce cuando hay varios posibles.",
        "<b>Todas las inestables</b> (original): prueba círculos cada vez más hacia "
        "tierra (cada 0,5 m) y rompe cada uno con FS &lt; 1, por orden, sobre el perfil ya "
        "roto. Puede encadenar varias roturas en un mes.<br>"
        "<b>La más profunda</b>: rompe el círculo inestable que llega más lejos tierra "
        "adentro, que es lo que describe la tesis.<br>"
        "<b>La de menor FS</b>: rompe el más inestable, que es el criterio clásico.",
    ),
    # ---------------------------------------------------------------- oleaje
    "oleaje": (
        "Las olas que llegan en cada estación: su altura, su periodo y qué porcentaje "
        "del tiempo se dan.",
        "<b>H</b>: altura de ola en aguas profundas (m). <b>T</b>: periodo (s). "
        "<b>Frecuencia</b>: % del tiempo en que se da esa ola. No suman 100 %: el resto es "
        "calma.<br><br>"
        "Para cada tipo de ola, el modelo calcula la longitud de onda L<sub>0</sub> = "
        "1,56·T² y el número de olas al mes, N<sub>o</sub> = (3600/T)·(frecuencia/100)·24·30, "
        "menos los días en que los derrubios protegen el pie. Propaga la ola desde una "
        "profundidad de L<sub>0</sub>/2 hacia la costa (asomeramiento), la rompe cuando "
        "H ≥ h/1,28 y calcula con ella la erosión del fondo, el impacto y el run-up.<br><br>"
        "Estaciones: primavera (abr–jun), verano (jul–sep) y otoño (oct–dic). En invierno "
        "(ene–mar) el hielo impide el oleaje y no se calcula. Valores por defecto: datos "
        "NOAA del lago Erie (tabla 10.iii de la tesis).",
    ),
    # ----------------------------------------------------------- nivel / mareas
    "nivel_estacional": (
        "Sube el nivel del lago en primavera y verano, como ocurre en los Grandes Lagos.",
        "Cada cambio de estación, el nivel pasa a ser «nivel del lago + subida de esa "
        "estación». Valores de la tesis (tabla 10.iii, cotas IGLD 1985): invierno y otoño "
        "173,99 m; primavera 174,39 m (+0,40); verano 174,29 m (+0,30). Las subidas se "
        "redondean a múltiplos de Δz. El Fortran original tiene esta variación comentada "
        "(desactivada).",
    ),
    "nivel_interanual": (
        "Un ciclo lento de subidas y bajadas del nivel que dura 15 años.",
        "A partir del mes 69, cada 9 meses el nivel cambia en un Δz: 5 subidas, 10 bajadas "
        "y 5 subidas (±0,5 m con Δz = 0,1 m), y vuelve a empezar. Se suma a la variación "
        "estacional. En el original está comentado y, además, perdía los cambios que caían "
        "en invierno o durante la protección por derrubios, de modo que el nivel derivaba "
        "hacia arriba. Aquí está corregido.",
    ),
    "serie_nivel": (
        "Para ver una serie de niveles medidos (por ejemplo, marea.dat) junto al nivel "
        "que usa el modelo. Solo se dibuja: no interviene en el cálculo.",
        "Formato: una columna (nivel) o dos (tiempo y nivel), separadas por espacios, "
        "comas o punto y coma. Las líneas que no son numéricas, como las cabeceras, se "
        "ignoran.",
    ),
    # --------------------------------------------------------- linea temporal
    "linea_temporal": (
        "Arrastra la barra para ver el perfil mes a mes. Los colores indican qué pasa "
        "cada mes.",
        "<b>Gris oscuro</b>: mes calculado normalmente (oleaje y erosión).<br>"
        "<b>Azul hielo</b>: invierno (enero–marzo). El lago está helado, no hay oleaje y el "
        "modelo salta esos meses sin erosión.<br>"
        "<b>Arena</b>: protección por derrubios. Tras una rotura, el material caído queda al "
        "pie del acantilado y el oleaje tarda volumen/M<sub>f</sub> días en llevárselo; "
        "mientras tanto el pie está protegido y el modelo salta esos meses sin erosión.<br>"
        "<b>Rayas rojas</b>: meses con rotura (más altas si hay varias).<br>"
        "<b>Línea cian bajo la barra</b>: tramo ya recorrido.<br><br>"
        "Al pasar el ratón por la barra se ve qué ocurre en cada mes. Teclas: ← → (un "
        "mes), RePág/AvPág (un año), Inicio/Fin.",
    ),
    # -------------------------------------------------------- otras pestanas
    "perfil": (
        "Puntos que definen la forma inicial del acantilado, de tierra (arriba) a lago "
        "(abajo).",
        "Cota z y distancia horizontal y en metros, con la cota siempre decreciente. El "
        "modelo interpola el perfil cada Δz de cota. Tipos de punto:<br>"
        "<b>Cresta del talud</b>: borde superior del acantilado.<br>"
        "<b>Base del talud</b>: pie del acantilado, donde arrancan las roturas circulares.<br>"
        "<b>Punto de rezume</b>: donde el nivel freático sale a la cara (solo con el tipo "
        "de agua 1).<br>"
        "<b>Punto de estudio</b>: puntos cuyo retroceso se guarda en "
        "<i>timerecesion.out</i>.<br><br>"
        "El perfil debe llegar, por debajo del nivel del lago, al menos a media longitud "
        "de onda de la ola más larga (unos 45 m con el oleaje por defecto).",
    ),
    "estratos": (
        "Capas de terreno de arriba abajo, con su resistencia y su espesor.",
        "<b>UCS</b>: resistencia a compresión simple. El código la lee pero no la usa; la "
        "tesis la emplea para estimar la resistencia a tracción (UCS/25 a UCS/5).<br>"
        "<b>Resist. corte</b>: en el código se usa como <b>resistencia a tracción</b> en el "
        "vuelco de bloques descalzados (tesis: 8 kPa).<br>"
        "<b>Cohesión</b> y <b>fricción</b>: parámetros efectivos (c′, φ′) para los "
        "deslizamientos circulares (tesis: 15–23 kPa y 25–28°).<br>"
        "<b>Peso esp.</b>: peso específico del terreno (kN/m³).<br>"
        "<b>Espesor</b>: grosor de la capa, contado desde la cota del primer punto del "
        "perfil. La suma debe cubrir toda la altura del perfil.",
    ),
}
