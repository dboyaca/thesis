# ¿Por qué la cruz del ápice no cae en el mismo nodo en los tres métodos?

## Reconstrucción de topografía corneal a partir de patrones de Hartmann

**Experimento:** `23-11-KeratoconusMC` · paciente `maconst01` (OD) · referencia `hpref.tif`
**Notebook:** `codefiles/KeratoconusAnalysis.ipynb`. Los números de celda son el índice 0-based del
JSON. Las secciones son "Mapas de nivel para BFTE" (107), "Zernike y elevación sobre el mejor ajuste"
(109), "Ajuste geométrico tipo Pentacam" (111) y "Enhanced Best Fit Sphere" (113).
**Convención:** el origen es el centroide de los spots de referencia (eje del instrumento) y +Y es
superior. El paso de la rejilla es p = 0.409 mm.

Todos los números de este documento salen de ejecutar el notebook en su estado actual. Los scripts
de verificación están descritos en la sección 8.

---

## 1. Lo que se observa

| Figura (celda)                                    | Southwell      | Campbell       | Malacara       | Coinciden |
| ------------------------------------------------- | -------------- | -------------- | -------------- | --------- |
| Curvas de nivel, fila BFTE (107) y Zernike (109)  | (−0.40, −1.64) | (−0.40, −2.04) | (−0.40, −1.64) | S = M     |
| Ápice irregular, EBFS (113)                       | (−0.40, −2.04) | (−0.81, −2.04) | (−0.81, −2.04) | C = M     |
| BFTE geométrico tipo Pentacam (111)               | (+0.01, −1.64) | (−0.40, −1.64) | (−2.85, −1.66) | ninguno   |

Coordenadas en mm. Tres hechos saltan a la vista:

1. Salvo Malacara en la celda 111, las cruces están **a un nodo de distancia** (0.41 mm), ya sea en
   vertical o en horizontal. Las tres caen en el mismo bloque de 2×2 o 2×3 nodos de la zona inferior.
2. **La pareja que coincide cambia de una figura a otra**: S = M en 107/109 y C = M en 113. Si dos
   métodos coincidieran por ser "mejores", lo harían siempre los mismos dos.
3. La única discrepancia grande (Malacara en 111, en el borde) es de otra naturaleza. Se trata en la
   sección 5.

---

## 2. Cómo se coloca la cruz

La función `apice_cono` (celda 87) hace lo siguiente:

1. Toma el residuo de elevación después de restar el mejor ajuste (tórico).
2. Lo suaviza con un promedio 3×3 sobre los nodos interiores.
3. Devuelve **el nodo** donde ese promedio es máximo.

Así que la cruz:

- **solo puede caer sobre un nodo de la rejilla**, y su resolución es de ±p/2 ≈ ±0.2 mm;
- es un **argmax**: si dos nodos tienen casi el mismo valor, cualquier perturbación mínima decide
  cuál gana.

---

## 3. Razón principal: el pico del cono es una meseta y la cruz resuelve un empate

Valor del residuo suavizado en los cuatro nodos más altos de cada método (celdas 107/109):

| Método    | 1.º nodo                 | 2.º nodo                 | Diferencia 1.º − 2.º |
| --------- | ------------------------ | ------------------------ | -------------------- |
| Southwell | (−0.40, −1.64): 11.45 µm | (−0.40, −2.04): 11.42 µm | **0.03 µm**          |
| Campbell  | (−0.40, −2.04): 10.76 µm | (−0.40, −1.64): 10.76 µm | **0.00 µm**          |
| Malacara  | (−0.40, −1.64): 9.46 µm  | (−0.81, −1.64): 9.39 µm  | **0.07 µm**          |

Hay que compararlo con dos escalas de error del propio experimento:

- **Error de cierre del campo de pendientes:** 0.67 µm rms por celda, con un máximo de 2.14 µm
  (sección 4.1).
- **Diferencia entre las superficies reconstruidas:** de 1.3 a 1.9 µm rms entre métodos.

> Lo que separa al primer nodo del segundo (de 0.00 a 0.07 µm) es **entre 10 y 100 veces menor** que
> la incertidumbre de la elevación. Con estos datos, el ápice no se puede localizar con precisión de un
> solo nodo. Lo que sí está bien determinado es la **región** del ápice: el bloque de nodos
> x ∈ [−0.81, +0.01] mm, y ∈ [−2.04, −1.64] mm, que aparece entre los cuatro primeros nodos en los
> tres métodos.

Para argumentar en la tesis: el cono es un rasgo ancho, de unos 2 mm, es decir de 4 a 5 nodos. Cerca
de su cima la superficie es casi plana en la escala de la rejilla. Por eso que la cruz salte un nodo
entre métodos es lo **esperable** y no indica un error de ninguno.

### Con un ápice subnodo, los tres métodos coinciden

Si en lugar del argmax se calcula el **centroide ponderado** de los nodos que están a menos de 1 µm
del máximo:

| Método    | Centroide del ápice (mm) |
| --------- | ------------------------ |
| Southwell | (−0.32, −1.85)           |
| Campbell  | (−0.44, −1.85)           |
| Malacara  | (−0.50, −1.76)           |

La dispersión entre métodos baja a **unos 0.2 mm, la mitad de un paso de rejilla**. Con un umbral de
0.5 o 2 µm las posiciones cambian menos de 0.1 mm. En la tesis puede reportarse que los tres métodos
sitúan el ápice del cono en (−0.4 ± 0.2, −1.8 ± 0.2) mm, es decir **unos 1.8 mm por debajo del eje,
en la zona inferior**, que es la localización típica del queratocono.

---

## 4. Por qué las tres superficies difieren un poco

La meseta explica por qué una diferencia pequeña mueve la cruz. Queda explicar de dónde salen esas
diferencias pequeñas. Hay tres causas, todas medidas.

### 4.1 El campo medido no es exactamente un gradiente

Un campo de pendientes que viene de una superficie real tiene rotacional nulo. La integral de las
pendientes alrededor de cualquier celda cerrada debe dar cero. Con los datos medidos, esa circulación
(regla del trapecio en las cuatro aristas de cada celda) vale **0.67 µm rms y hasta 2.14 µm**, en 216
celdas.

Las causas son el error de centroide (K-means sobre spots de pocos píxeles), la discretización de los
spots y cualquier desajuste del emparejamiento spot ↔ nodo. Esa parte no integrable **no tiene una
única integral**, y cada método la reparte de forma distinta:

| Método        | Cómo trata la parte no integrable                                                                                                                                   |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Southwell** | Mínimos cuadrados globales: busca la superficie cuyo gradiente está más cerca del medido y reparte el error por toda la apertura (Southwell, 1980).                  |
| **Campbell**  | Una pasada por cuadrantes desde el centro, promediando las llegadas a cada nodo: es un promedio de varios caminos que compensa el error parcialmente (Campbell, 2022). |
| **Malacara**  | Integración por un **único camino** (fila central en x y luego cada columna en y): la superficie depende del camino elegido.                                          |

**Prueba directa de la dependencia del camino.** Se integró Malacara con el camino al revés: primero
la columna central en y y luego cada fila en x. Resultado:

- La superficie cambia **2.7 µm rms, con un máximo de 7.6 µm**.
- **La cruz salta de (−0.40, −1.64) a (−0.40, −2.04)**, el mismo nodo que elige Campbell.

> Por tanto, que Southwell y Malacara coincidan en la celda 107/109 depende de **qué camino eligió
> Malacara**. No refleja que esos dos métodos sean más parecidos o más exactos.

### 4.2 Los errores de los métodos de camino crecen con la distancia al centro

Diferencia media (en valor absoluto) respecto a Southwell, según la distancia al nodo central:

| Distancia al centro | Campbell − Southwell | Malacara − Southwell |
| ------------------- | -------------------- | -------------------- |
| 0 a 3 nodos         | 0.49 µm              | 0.75 µm              |
| 3 a 6 nodos         | 0.78 µm              | 1.08 µm              |
| 6 a 10 nodos        | 1.22 µm              | 1.99 µm              |

Campbell y Malacara fijan la altura en el nodo central y propagan hacia afuera, así que el error se
acumula como un paseo aleatorio. Además, la diferencia de Malacara está **organizada por filas**: la
dispersión de las medias por fila es de 1.43 µm, frente a 0.19 µm por columna. Es la huella de integrar
cada columna en y partiendo de la fila central: el error acumulado varía de una fila a otra.

El ápice está a unos 4.5 nodos del centro. A esa distancia el error relativo entre métodos (≈ 1 µm) ya
es mucho mayor que el empate de la sección 3.

### 4.3 Discretizaciones distintas de la misma ecuación

- **Southwell y Malacara** relacionan nodos vecinos con la misma regla del trapecio:
  $z_{j+1}-z_j = \tfrac{p}{2}(s_j+s_{j+1})$. Southwell la impone en mínimos cuadrados sobre todas las
  aristas; Malacara solo sobre las aristas de su camino.
- **Campbell** usa pendientes **promediadas sobre las cuatro esquinas** de cada celda, más un término
  de derivada cruzada. Eso equivale a un filtro paso bajo de la pendiente, que **aplana y ensancha** el
  pico.

Esto encaja con que Campbell tenga la meseta más plana de los tres (empate exacto: 10.76 = 10.76 µm)
y la cima más baja que Southwell (10.76 frente a 11.45 µm).

Aun así, en conjunto **Southwell y Campbell son los más parecidos** (r = 0.977 y 1.27 µm rms, frente a
r = 0.946 y 1.89 µm rms entre Southwell y Malacara). La similitud de discretización entre Southwell y
Malacara **no basta** para explicar su coincidencia en la cruz: la explicación correcta es el empate
de la sección 3.

---

## 5. Caso aparte: Malacara en el ajuste geométrico tipo Pentacam (celda 111)

En esta figura la cruz de Malacara cae en (−2.85, −1.66) mm, **en el borde de la apertura**. No es un
desplazamiento de un nodo como en los demás casos. Los cuatro nodos más altos del residuo de Malacara
sobre el BFTE geométrico son:

| Nodo (mm)      | Residuo suavizado |
| -------------- | ----------------- |
| (−2.85, −1.66) | 6.38 µm           |
| (−1.24, +2.87) | 6.36 µm           |
| (−0.00, +3.29) | 6.16 µm           |
| (−1.67, +2.81) | 5.99 µm           |

Están **repartidos por la periferia, con valores casi iguales**, y ninguno está en la zona del cono.
La explicación es esta:

1. El BFTE geométrico tiene radios por meridiano, eje y ápice desplazable. Absorbe más forma que el
   polinomio tórico de la celda 107, y deja al cono con un contraste de solo unos 7 µm (Southwell:
   7.3 µm).
2. El error acumulado de Malacara en la periferia es de unos 2 µm en promedio, y de hasta 7.6 µm según
   el camino (sección 4). Es **comparable al contraste del cono**.
3. El máximo deja de ser el cono y pasa a ser el lugar donde más error se acumuló.

Southwell (+0.01, −1.64) y Campbell (−0.40, −1.64) siguen siendo nodos vecinos dentro de la misma
meseta inferior. Conclusión para la tesis: **con el ajuste geométrico, Malacara no resuelve el cono**,
y su cruz en esa figura es un artefacto de integración.

> **No confundir** con la línea "centro en máx. elevación BFS (−2.85, −1.66)" de la celda 113. Ahí los
> **tres** métodos coinciden en el borde porque esa magnitud es el máximo de la elevación sobre la
> **esfera**, que incluye el astigmatismo regular y siempre crece hacia la periferia. Es otra cantidad,
> no el ápice del cono.

---

## 6. Argumentos que **no** conviene usar

| Afirmación                                                         | Por qué no se sostiene                                                                                                                     |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| "Southwell y Malacara coinciden porque son los más precisos."      | Al invertir el camino de Malacara, su cruz pasa al nodo de Campbell (sección 4.1). Además, la pareja que coincide cambia entre figuras.      |
| "Campbell está mal porque se separa de los otros dos."             | Sus dos nodos más altos empatan a 10.76 µm. Su cruz cae en el otro nodo del empate, y es el método más parecido a Southwell en conjunto.     |
| "El cambio de paso de 0.4 a 0.409 mm movió las cruces."            | Se verificó con ambos pasos: las cruces caen en los mismos nodos. El paso solo escala la elevación un 2.2 %.                                |
| "El ápice está en (−0.40, −1.64) mm." (dado como un punto exacto)  | La resolución es de ±0.2 mm y el empate impide elegir un nodo. Conviene dar el centroide con su incertidumbre (sección 3).                  |

---

## 7. Texto sugerido para la tesis

> Los tres métodos de integración zonal localizan el ápice del cono en la misma región, en el cuadrante
> inferior y unos 1.8 mm por debajo del eje del instrumento. Las cruces de la Figura X difieren en a lo
> sumo un nodo de la rejilla (0.41 mm). Esa diferencia no indica una discrepancia real entre métodos:
> es consecuencia de que el máximo del residuo forma una meseta. Los dos nodos más altos de cada
> reconstrucción difieren en menos de 0.1 µm, mientras que el error de cierre del campo de pendientes
> es de 0.67 µm rms por celda y las superficies reconstruidas difieren entre sí de 1.3 a 1.9 µm rms.
> Bajo esas condiciones, la elección del nodo máximo queda dominada por diferencias de integración que
> son inherentes a cada método. Southwell minimiza el error en mínimos cuadrados, Campbell promedia
> varios caminos y Malacara depende del camino elegido. Invertir el orden de integración de Malacara
> basta para trasladar su máximo al nodo que elige Campbell. Calculando el ápice como centroide
> ponderado de la meseta, las tres posiciones coinciden dentro de 0.2 mm, la mitad del paso de
> muestreo. En el ajuste al elipsoide tórico geométrico, el contraste residual del cono (≈ 7 µm) es
> comparable al error acumulado de la integración por caminos en la periferia, y el máximo de Malacara
> se desplaza al borde de la apertura. Por eso, para localizar el ápice se toma como referencia la
> solución de mínimos cuadrados (Southwell).

---

## 8. Cómo se obtuvieron los números

- Se ejecutaron las celdas de definición del notebook hasta la 105, y después las 107, 109, 111 y 113.
  Se usaron las funciones del propio notebook (`campo_pendientes`, `recon_*`, `quitar_mejor_ajuste`,
  `apice_cono`) con p = `p_med` = 0.4089 mm. Para comparar se repitió con p = 0.4.
- **Meseta:** el mismo promedio 3×3 de `apice_cono`, ordenando los cuatro nodos más altos.
- **Error de cierre:** la circulación trapezoidal de las pendientes alrededor de cada celda con los
  cuatro nodos válidos.
- **Malacara con camino invertido:** `recon_malacara(Sy.T, Sx.T, M.T, p)` transpuesto. Integra primero
  la columna central y luego las filas.
- **Centroide subnodo:** el promedio de las posiciones de los nodos con residuo suavizado ≥ máximo − 1
  µm, ponderado por (valor − umbral).

**Salvedad pendiente.** Las figuras suponen que la foto no está invertida (+Y = superior, igual que la
imagen). Falta confirmar con el montaje si un espejo o un divisor de haz invierte la imagen en X, en Y o
en ambos. Hasta entonces, "inferior" significa "hacia abajo en la foto". Antes de escribir "nasal" o
"temporal" también hay que saber qué ojo es y cómo está orientada la cámara. Nada de esto afecta al
argumento central, que es la **concordancia entre métodos**.

---

## Referencias

- W. H. Southwell, "Wave-front estimation from wave-front slope measurements," *J. Opt. Soc. Am.*
  **70**(8), 998–1006 (1980).
- C. E. Campbell, "Zonal reconstruction of surfaces from gradient fields," *Opt. Eng.* **61**(12),
  121802 (2022).
- D. Malacara (ed.), *Optical Shop Testing*, 3rd ed., Wiley (2007), capítulos sobre pruebas de
  Hartmann e integración de pendientes.
