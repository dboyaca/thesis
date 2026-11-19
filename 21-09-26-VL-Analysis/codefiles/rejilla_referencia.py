"""
Rejilla de referencia sintetica y reindexado canonico de centroides.
=====================================================================

Para que sirve
--------------
El experimento 21-09-26 tiene imagenes de 1024x1024 pero NO tiene una imagen de
referencia propia.  La referencia del experimento 23-11 (`hpref.tif`, 480x640,
paso 19.95 px) no vale: el cociente de pasos contra estas imagenes es 1.855, o
sea +86 %, mientras que la senal corneal real es del orden del 9 %.  Reescalarla
absorberia en el factor de escala justo la magnitud que se quiere medir.

La alternativa es generar la referencia, igual que se hizo en 22-11 (cuyo
`22-11-RefCentroids.txt` es una rejilla generada, no medida).

Que hace este script
--------------------
1. Carga los centroides medidos.
2. Descarta blobs de ruido aislados.
3. Ajusta la rejilla regular de mejor ajuste (centro, paso, rotacion).
4. Reasigna `pos` a partir de los indices de esa rejilla, no de un contador
   secuencial.  Esto corrige un fallo silencioso: si la deteccion pierde un
   spot, el contador corre todos los indices siguientes y el mapeo a mano de
   los notebooks de analisis asigna cada punto a la celda del vecino.
5. Genera la rejilla de referencia sintetica con el mismo paso, centro,
   rotacion y apertura.
6. Escribe los dos ficheros y un diagnostico.

QUE SE PIERDE Y QUE SE CONSERVA
-------------------------------
El paso de la rejilla sintetica se ajusta AL PROPIO PATRON medido.  Eso fuerza a
cero la diferencia media de magnificacion, o sea el termino de DESENFOQUE.  En
consecuencia:

  * NO se recupera la queratometria absoluta (K, radio corneal).
  * SI se conservan astigmatismo, coma, trefoil y el defecto irregular, que es
    donde vive el cono.

Para recuperar la queratometria absoluta hace falta la geometria del montaje
(paso fisico de la rejilla del proyector y magnificacion), y entonces se fija
PASO_REF a mano en vez de ajustarlo.

Centrar la referencia en el centro ajustado del patron (opcion por defecto)
elimina por construccion el termino de INCLINACION.  Eso es deseable: la
inclinacion es descentramiento entre capturas, no forma corneal, y en el
experimento 23-11 suponia 341 de los 361 um del mapa crudo.
"""
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
#  CONFIGURACION
# ---------------------------------------------------------------------------
ENTRADA = "./../data/21-09-26-od.txt"
SALIDA_TEST = "./../data/21-09-26-od-corregido.txt"
SALIDA_REF = "./../data/21-09-26-ref-sintetica.txt"

#  Donde se ancla la rejilla de referencia:
#    'ajustado'  centro de la rejilla de mejor ajuste.  Anula el desplazamiento
#                MEDIO (gauge de minimos cuadrados): no hay inclinacion global.
#    'vertice'   nodo central del patron medido, que en estas capturas es el
#                punto que se coloca a mano.  Anula el desplazamiento EN EL
#                CENTRO: el vertice corneal queda como origen de elevaciones.
#                Es la convencion de apice fijo de los topografos.
#    'imagen'    centro del lienzo (512, 512).
#    (x, y)      explicito.
CENTRO_REF = 'vertice'
PASO_REF = None           # None -> el paso ajustado.  Un numero -> se impone
                          #         (necesario para queratometria absoluta)
ROTAR_REF = True          # aplicar tambien la rotacion ajustada
TAM_IMAGEN = 1024

# Apertura canonica de 19x19 de los notebooks de analisis: 253 puntos.
APERTURA = {0: (9, 9),   1: (5, 13),  2: (4, 14),  3: (3, 15),  4: (2, 16),
            5: (1, 17),  6: (1, 17),  7: (1, 17),  8: (1, 17),  9: (0, 18),
            10: (1, 17), 11: (1, 17), 12: (1, 17), 13: (1, 17),
            14: (2, 16), 15: (3, 15), 16: (4, 14), 17: (5, 13), 18: (9, 9)}


def celdas_canonicas():
    """(fila, columna) en el orden que produce imgprocessing: columna a columna
       de izquierda a derecha, y ascendente dentro de cada columna."""
    return [(r, c) for c in range(19)
            for r in range(APERTURA[c][0], APERTURA[c][1] + 1)]


# ---------------------------------------------------------------------------
#  AJUSTE DE LA REJILLA REGULAR
# ---------------------------------------------------------------------------
def asignar_indices(xy, paso, centro, rot):
    """Celda entera (u, v) mas cercana para cada punto medido."""
    dx = xy[:, 0] - centro[0]
    dy = xy[:, 1] - centro[1]
    ct, st = np.cos(-rot), np.sin(-rot)
    return (np.rint((dx * ct - dy * st) / paso).astype(int),
            np.rint((dx * st + dy * ct) / paso).astype(int))


def modelo(u, v, centro, paso, rot):
    ct, st = np.cos(rot), np.sin(rot)
    return np.column_stack([centro[0] + paso * (u * ct - v * st),
                            centro[1] + paso * (u * st + v * ct)])


def ajustar_rejilla(xy, paso0, centro0, n_iter=15):
    """Ajusta (cx, cy, paso, rotacion) alternando asignacion de indices y lstsq.

    El modelo  x = cx + paso*(u*cos r - v*sin r)
               y = cy + paso*(u*sin r + v*cos r)
    reparametrizado con (a, b) = paso*(cos r, sin r) es LINEAL en (cx, cy, a, b),
    asi que cada iteracion es un unico problema de minimos cuadrados:

        | 1  0   u  -v |  | cx |     | x |
        | 0  1   v   u |  | cy |  =  | y |
                          | a  |
                          | b  |

    apiladas las N ecuaciones en x sobre las N en y.  Despues se recupera
    paso = hypot(a, b) y rotacion = atan2(b, a).
    """
    centro, paso, rot = np.array(centro0, float), float(paso0), 0.0
    for _ in range(n_iter):
        u, v = asignar_indices(xy, paso, centro, rot)
        N = len(xy)
        A = np.zeros((2 * N, 4))
        A[:N, 0] = 1.0; A[:N, 2] = u; A[:N, 3] = -v
        A[N:, 1] = 1.0; A[N:, 2] = v; A[N:, 3] = u
        b = np.concatenate([xy[:, 0], xy[:, 1]])
        (cx, cy, a_, b_), *_ = np.linalg.lstsq(A, b, rcond=None)
        centro, paso, rot = np.array([cx, cy]), np.hypot(a_, b_), np.arctan2(b_, a_)

    # Re-centrar el origen de indices: el ajuste puede converger desplazado un
    # numero entero de celdas.  Se traslada para que las ocupadas queden
    # centradas en (0,0) y se corrige el centro en consecuencia.
    u, v = asignar_indices(xy, paso, centro, rot)
    du = int(np.rint((u.min() + u.max()) / 2.0))
    dv = int(np.rint((v.min() + v.max()) / 2.0))
    if du or dv:
        ct, st = np.cos(rot), np.sin(rot)
        centro = centro + paso * np.array([du * ct - dv * st, du * st + dv * ct])
        u, v = asignar_indices(xy, paso, centro, rot)

    return dict(centro=centro, paso=paso, rot=rot, u=u, v=v,
                residuo=xy - modelo(u, v, centro, paso, rot))


def limpiar(xy, factor=2.0):
    """Descarta puntos aislados: sin ningun vecino a menos de `factor` pasos.

    Un solo blob de ruido lejos del patron desplaza media y extension, y con
    ellas el paso inicial y el centro, haciendo que el ajuste converja a una
    rejilla equivocada.  En este experimento habia uno en (1000.3, 0.3), a
    332 px de su vecino mas cercano.
    """
    from scipy.spatial import cKDTree
    d, _ = cKDTree(xy).query(xy, k=2)
    nn = d[:, 1]
    return nn < factor * np.median(nn), nn


# ---------------------------------------------------------------------------
#  PROGRAMA
# ---------------------------------------------------------------------------
def main():
    df = pd.read_csv(ENTRADA, sep=r'\s+', header=None)
    df.columns = ['pos', 'x', 'y']
    xy = df[['x', 'y']].values
    print("entrada: %s  (%d puntos)" % (ENTRADA, len(df)))

    bueno, nn = limpiar(xy)
    for i in np.where(~bueno)[0]:
        print("   descartado pos %-4d (%8.2f, %8.2f): vecino mas cercano a %.1f px"
              % (df.pos[i], xy[i, 0], xy[i, 1], nn[i]))

    aj = ajustar_rejilla(xy[bueno], np.median(nn[bueno]),
                         np.median(xy[bueno], axis=0))
    res = np.hypot(*aj['residuo'].T)
    print("\nREJILLA REGULAR DE MEJOR AJUSTE")
    print("   centro   : (%.2f, %.2f) px" % tuple(aj['centro']))
    print("   paso     : %.4f px" % aj['paso'])
    print("   rotacion : %+.4f grados" % np.degrees(aj['rot']))
    print("   celdas   : u %d..%d,  v %d..%d   colisiones: %d"
          % (aj['u'].min(), aj['u'].max(), aj['v'].min(), aj['v'].max(),
             len(aj['u']) - len(set(zip(aj['u'], aj['v'])))))
    print("   residuo  : rms %.2f px   max %.2f px   <- campo de aberracion" % (res.std(), res.max()))

    # --- reindexado canonico ------------------------------------------------
    canon = celdas_canonicas()
    orden = {celda: k for k, celda in enumerate(canon, start=1)}
    fil, col = aj['v'] + 9, aj['u'] + 9
    xy_ok = xy[bueno]

    filas, fuera = [], 0
    for f, c, (X, Y) in zip(fil, col, xy_ok):
        if (f, c) not in orden:
            fuera += 1
            continue
        filas.append((orden[(f, c)], X, Y))
    test = pd.DataFrame(filas, columns=['pos', 'x', 'y']).sort_values('pos')

    desfase = sum(1 for p_nuevo, p_viejo in zip(test.pos.values, df.pos[bueno].values)
                  if p_nuevo != p_viejo)
    print("\nREINDEXADO CANONICO")
    print("   puntos reindexados: %d   fuera de la apertura canonica: %d" % (len(test), fuera))
    print("   celdas canonicas sin punto medido: %s"
          % sorted(set(orden.values()) - set(test.pos.values)))
    print("   posiciones que CAMBIAN respecto al archivo original: %d de %d"
          % (desfase, len(test)))
    if desfase:
        print("   -> el archivo original tenia el mapeo pos->celda corrido; los")
        print("      notebooks de analisis habrian asignado puntos a celdas vecinas.")

    # --- referencia sintetica ----------------------------------------------
    paso = aj['paso'] if PASO_REF is None else float(PASO_REF)
    rot = aj['rot'] if ROTAR_REF else 0.0
    if CENTRO_REF == 'ajustado':
        centro = aj['centro']
    elif CENTRO_REF == 'vertice':
        #  Nodo central de la apertura canonica = celda (fila 9, columna 9).
        #  En estas capturas es el punto que se inserta a mano en imgprocessing.
        centro_pos = canon.index((9, 9)) + 1
        fila_c = test[test.pos == centro_pos]
        if len(fila_c) == 0:
            raise ValueError("no hay punto medido en la celda central (9,9): "
                             "no se puede anclar la referencia al vertice")
        centro = np.array([fila_c.x.values[0], fila_c.y.values[0]])
    elif CENTRO_REF == 'imagen':
        centro = np.array([TAM_IMAGEN / 2.0, TAM_IMAGEN / 2.0])
    else:
        centro = np.array(CENTRO_REF, float)

    ref = pd.DataFrame(
        [(k, *modelo(np.array([c - 9]), np.array([r - 9]), centro, paso, rot)[0])
         for k, (r, c) in enumerate(canon, start=1)],
        columns=['pos', 'x', 'y'])

    print("\nREFERENCIA SINTETICA")
    print("   centro %s : (%.2f, %.2f) px" % (CENTRO_REF, *centro))
    print("   paso %s : %.4f px" % ("ajustado" if PASO_REF is None else "impuesto", paso))
    print("   rotacion : %+.4f grados" % np.degrees(rot))
    print("   %d puntos,  x %.1f-%.1f   y %.1f-%.1f"
          % (len(ref), ref.x.min(), ref.x.max(), ref.y.min(), ref.y.max()))
    if not (0 <= ref.x.min() and ref.x.max() < TAM_IMAGEN
            and 0 <= ref.y.min() and ref.y.max() < TAM_IMAGEN):
        print("   AVISO: la referencia se sale del lienzo de %d px." % TAM_IMAGEN)

    # --- desplazamientos ----------------------------------------------------
    m = test.merge(ref, on='pos', suffixes=('_t', '_r'))
    dx, dy = m.x_t - m.x_r, m.y_t - m.y_r
    print("\nCAMPO DE DESPLAZAMIENTO test - referencia")
    print("   %d nodos en comun" % len(m))
    print("   modulo: media %.2f px   max %.2f px" % (np.hypot(dx, dy).mean(),
                                                      np.hypot(dx, dy).max()))
    print("   traslacion rigida media: (%+.3f, %+.3f) px  <- cero por construccion"
          % (dx.mean(), dy.mean()))

    for ruta, tabla in [(SALIDA_TEST, test), (SALIDA_REF, ref)]:
        with open(ruta, 'w') as fh:
            for p, X, Y in tabla.itertuples(index=False):
                fh.write("%d %.8f %.8f\n" % (p, X, Y))
        print("escrito: %s  (%d puntos)" % (ruta, len(tabla)))


if __name__ == '__main__':
    main()
