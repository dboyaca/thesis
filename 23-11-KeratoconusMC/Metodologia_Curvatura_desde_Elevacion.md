# De la elevación a la curvatura

## Cómo se obtienen los mapas sagital y tangencial a partir de un mapa de elevación

**Experimento:** `23-11-KeratoconusMC` · paciente `maconst01` (OD)
**Notebooks:** `KeratoconusAnalysis.ipynb` (celda 106) y `KeratoconusAnalysisBFS.ipynb` (celda 106)

---

## 1. La idea de fondo: la curvatura **es** una derivada segunda

La sorpresa desaparece en cuanto se recuerda que curvatura y elevación no son dos medidas
independientes: **la curvatura está contenida en la elevación**, igual que la velocidad está
contenida en la trayectoria.

Para una curva plana $z(x)$, la curvatura en un punto es

$$k = \frac{z''}{\left(1+z'^2\right)^{3/2}}$$

No hay información nueva: $k$ se calcula derivando dos veces. El denominador corrige el hecho
de que $x$ no es la longitud de arco; si la superficie fuera casi plana ($z'\approx 0$) quedaría
simplemente $k\approx z''$.

Lo mismo vale para una superficie $z(x,y)$, solo que entonces la curvatura depende de la
dirección en que se mida, y hay que elegir cuáles reportar.

---

## 2. Las dos curvaturas que usa la topografía corneal

En la córnea se toma el eje óptico como referencia y se trabaja en coordenadas polares
$(r,\theta)$ centradas en el vértice. En cada punto se definen dos direcciones:

- el **meridiano**, la dirección radial (hacia el vértice),
- la **dirección sagital**, perpendicular a ella (azimutal).

### Curvatura tangencial (o meridional, o instantánea)

Es la curvatura de la curva que se obtiene al cortar la superficie con el plano que contiene el
eje y el punto. Es la fórmula de curva plana aplicada al perfil radial:

$$\boxed{\;k_{\text{tangencial}} = \frac{\partial^2 z/\partial r^2}{\left[1+(\partial z/\partial r)^2\right]^{3/2}}\;}$$

Es **local**: refleja la forma de la superficie justo ahí. Por eso es la que mejor delimita un
cono y la que se usa para localizarlo.

### Curvatura sagital (o axial)

No es una derivada segunda, sino una construcción geométrica: se traza la normal a la superficie
en el punto y se mide la distancia hasta que corta el eje óptico. Ese es el radio sagital.

$$\boxed{\;k_{\text{sagital}} = \frac{\partial z/\partial r}{r\,\sqrt{1+(\partial z/\partial r)^2}}\;}$$

Solo necesita la **primera** derivada. Es un promedio del comportamiento entre el vértice y el
punto, así que sale más suave que la tangencial y difumina los rasgos locales. Es la que
tradicionalmente reportan los queratómetros.

En el vértice ambas coinciden, y el código lo impone explícitamente para evitar el $0/0$:

```python
k_s = np.where(r < 1e-9, k_t, k_s)
```

### Cómo se pasa de $(x,y)$ a la dirección radial

Las derivadas se calculan en la malla cartesiana y luego se proyectan sobre el meridiano con la
regla de la cadena. Con $\theta$ el ángulo polar, $c=\cos\theta$, $s=\sin\theta$:

$$z_r = z_x c + z_y s$$
$$z_{rr} = z_{xx}c^2 + 2z_{xy}cs + z_{yy}s^2$$

La segunda es la forma cuadrática del hessiano evaluada en la dirección radial. En el código:

```python
th = np.arctan2(y, x); ct, st = np.cos(th), np.sin(th)
z_r  = zx*ct + zy*st
z_rr = zxx*ct**2 + 2*zxy*ct*st + zyy*st**2
```

### Unidades

La curvatura sale en mm⁻¹. Se convierte a dioptrías queratométricas con

$$D = \frac{n-1}{R} = \frac{0.3375}{R\,[\text{m}]} = \frac{337.5}{R\,[\text{mm}]} = 337.5\,k$$

usando el índice queratométrico convencional $n = 1.3375$.

---

## 3. El paso que no se ve: **el mapa de elevación no es la superficie**

Aquí está lo que hace que el asunto no sea trivial, y es la parte que conviene tener clara.

Un mapa de elevación es una **diferencia**: cuánto se separa la córnea de una superficie de
referencia. La córnea completa es

$$z_{\text{córnea}}(x,y) = z_{\text{referencia}}(x,y) + \Delta h(x,y)$$

y la curvatura hay que calcularla sobre $z_{\text{córnea}}$, no sobre $\Delta h$. Si derivaras
solo el residuo obtendrías la *perturbación* de curvatura, no la curvatura.

Como la derivada es lineal, las contribuciones se suman:

$$z_x = \underbrace{\frac{\partial z_{\text{ref}}}{\partial x}}_{\text{analítica}} + \underbrace{\frac{\partial \Delta h}{\partial x}}_{\text{diferencias finitas}}$$

**Y ese reparto es la decisión de diseño importante del código.** La referencia es una esfera (o
un toro), cuyas derivadas se conocen en forma cerrada. Para una esfera de radio $R$ centrada en
el eje, con $q=\sqrt{R^2-x^2-y^2}$:

$$z_x = \frac{x}{q}, \qquad z_{xx} = \frac{R^2-y^2}{q^3}, \qquad z_{xy} = \frac{xy}{q^3}$$

En el código:

```python
q = np.sqrt(np.maximum(R_ref**2 - x**2 - y**2, 1e-12))
zx  = x/q                    + dx
zy  = y/q                    + dy
zxx = (R_ref**2 - y**2)/q**3 + dxx
zyy = (R_ref**2 - x**2)/q**3 + dyy
zxy = x*y/q**3               + dxy
```

Las tres primeras columnas son exactas; solo `dx, dy, dxx, dyy, dxy` vienen de diferencias
finitas. El comentario de la función lo resume: *«La esfera de referencia aporta derivadas
exactas; FD solo sobre dh»*.

---

## 4. Por qué ese reparto es esencial: derivar dos veces amplifica el ruido

Esta es la razón práctica de todo lo anterior, y conviene cuantificarla.

Una derivada segunda por diferencias centradas

$$z_{xx} \approx \frac{z_{i+1}-2z_i+z_{i-1}}{p^2}$$

combina tres valores con coeficientes $(1,-2,1)$ y divide por $p^2$. Si cada nodo tiene ruido
independiente de desviación $\sigma_h$, el ruido de la derivada segunda es

$$\sigma_{z_{xx}} = \frac{\sqrt{1^2+2^2+1^2}}{p^2}\,\sigma_h = \frac{\sqrt{6}}{p^2}\,\sigma_h$$

Con $p = 0.4$ mm eso da $1/p^2 = 6.25\ \text{mm}^{-2}$, y en dioptrías:

$$337.5 \times \frac{\sqrt{6}}{0.16} \times 10^{-3} = 5.17\ \text{D por cada }\mu\text{m de ruido}$$

Medido empíricamente sobre estos datos, inyectando ruido gaussiano en la elevación y
observando cuánto se mueve la curvatura tangencial:

| ruido en $\Delta h$ | cambio en $k_{\text{tangencial}}$ | factor |
|---|---|---|
| 0.10 µm | 0.45 D | ×4.5 |
| 0.25 µm | 1.15 D | ×4.6 |
| 0.50 µm | 2.23 D | ×4.5 |
| 1.00 µm | 4.65 D | ×4.7 |

La medida (4.5–4.7 D/µm) coincide con la predicción teórica (5.17 D/µm) dentro de lo esperable.

**Y aquí está la consecuencia que hay que declarar:** el ruido real de la reconstrucción en este
experimento es de 0.91 µm rms (el residuo que no explica Zernike hasta orden 6). Eso solo ya
implica **≈ 4.7 D de incertidumbre punto a punto** en el mapa de curvatura.

Por eso no se deriva la superficie completa. Si se hiciera, el ruido se aplicaría sobre el
término dominante (≈ 45 D de la esfera) en vez de sobre la perturbación. Separando la parte
analítica, el ruido solo afecta a la parte pequeña — pero, como muestra la tabla, sigue siendo
el factor limitante y conviene interpretar los mapas de curvatura como suavizados, no punto a
punto.

---

## 5. Las diferencias finitas: la función `_fd`

Calcula las cinco derivadas del residuo. Lo que la hace algo más larga de lo esperado es que la
apertura tiene borde irregular: no todo nodo tiene vecinos a ambos lados.

```python
def _fd(z, m, p):
    nr, nc = z.shape
    nan = lambda: np.full((nr, nc), np.nan)
    zx, zy, zxx, zyy, zxy = nan(), nan(), nan(), nan(), nan()
    ok = lambda i, j: 0 <= i < nr and 0 <= j < nc and m[i, j]
    for i in range(nr):
        for j in range(nc):
            if not m[i, j]:
                continue
            if ok(i, j-1) and ok(i, j+1):
                zx[i, j]  = (z[i, j+1] - z[i, j-1]) / (2*p)
                zxx[i, j] = (z[i, j+1] - 2*z[i, j] + z[i, j-1]) / p**2
            elif ok(i, j+1) and ok(i, j+2):
                zx[i, j]  = (-3*z[i, j] + 4*z[i, j+1] - z[i, j+2]) / (2*p)
                zxx[i, j] = (z[i, j] - 2*z[i, j+1] + z[i, j+2]) / p**2
            elif ok(i, j-1) and ok(i, j-2):
                zx[i, j]  = ( 3*z[i, j] - 4*z[i, j-1] + z[i, j-2]) / (2*p)
                zxx[i, j] = (z[i, j] - 2*z[i, j-1] + z[i, j-2]) / p**2
```

(idéntico para `zy`, `zyy` recorriendo filas; y luego el término cruzado `zxy`.)

### Los tres casos

1. **Nodo interior** — diferencias **centradas**, error $O(p^2)$:

   $$z_x = \frac{z_{j+1}-z_{j-1}}{2p}, \qquad z_{xx} = \frac{z_{j+1}-2z_j+z_{j-1}}{p^2}$$

2. **Borde con dos vecinos hacia adelante** — fórmula **descentrada** de tres puntos, también
   $O(p^2)$:

   $$z_x = \frac{-3z_j+4z_{j+1}-z_{j+2}}{2p}$$

   Los coeficientes $(-3,4,-1)/2p$ salen de derivar la parábola que pasa por los tres nodos.

3. **Borde con dos vecinos hacia atrás** — la simétrica, con los signos cambiados.

Si no se cumple ninguno, el valor queda en `NaN` y ese nodo no aparece en el mapa. Es preferible
un hueco a un número inventado.

### El término cruzado

$$z_{xy} \approx \frac{z_{i+1,j+1}-z_{i+1,j-1}-z_{i-1,j+1}+z_{i-1,j-1}}{4p^2}$$

necesita las cuatro diagonales. Cuando alguna falta, el código recurre a derivar `zx` en la
dirección $y$, que es la misma cantidad calculada por otro camino:

```python
if   up and dn: zxy[i, j] = (zx[i+1, j] - zx[i-1, j]) / (2*p)
elif up:        zxy[i, j] = (zx[i+1, j] - zx[i, j]) / p
elif dn:        zxy[i, j] = (zx[i, j] - zx[i-1, j]) / p
else:           zxy[i, j] = 0.0
```

---

## 6. La función completa

```python
def sag_tan_fd(dh, mask, p=0.4, R_ref=7.78, center=None):
    """Curvaturas sagital (axial) y tangencial. dh = elevacion diferencial [mm].
       La esfera de referencia aporta derivadas exactas; FD solo sobre dh."""
    nr, nc = dh.shape
    ci, cj = (nr//2, nc//2) if center is None else center
    jj, ii = np.meshgrid(np.arange(nc), np.arange(nr))
    x = (jj - cj)*p; y = (ii - ci)*p
    m = mask.astype(bool)
    dx, dy, dxx, dyy, dxy = _fd(np.where(m, dh, 0.0), m, p)
    q = np.sqrt(np.maximum(R_ref**2 - x**2 - y**2, 1e-12))
    zx  = x/q                    + dx
    zy  = y/q                    + dy
    zxx = (R_ref**2 - y**2)/q**3 + dxx
    zyy = (R_ref**2 - x**2)/q**3 + dyy
    zxy = x*y/q**3               + dxy
    r  = np.sqrt(x**2 + y**2)
    th = np.arctan2(y, x); ct, st = np.cos(th), np.sin(th)
    z_r  = zx*ct + zy*st
    z_rr = zxx*ct**2 + 2*zxy*ct*st + zyy*st**2
    k_t = z_rr / (1 + z_r**2)**1.5
    with np.errstate(divide='ignore', invalid='ignore'):
        k_s = z_r / (r*np.sqrt(1 + z_r**2))
    k_s = np.where(r < 1e-9, k_t, k_s)
    return (np.where(m, 1/k_t, np.nan), np.where(m, 1/k_s, np.nan),
            (x.min()-p/2, x.max()+p/2, y.min()-p/2, y.max()+p/2))
```

Se lee de arriba abajo como la teoría de las secciones 2 y 3:

| líneas | qué hacen |
|---|---|
| `jj, ii = meshgrid` … | malla física en mm, centrada en el nodo central |
| `_fd(...)` | las cinco derivadas del residuo (§5) |
| `q = ...` y las cinco sumas | referencia analítica + residuo numérico (§3) |
| `z_r`, `z_rr` | proyección sobre el meridiano (§2) |
| `k_t`, `k_s` | las dos fórmulas de curvatura (§2) |
| `return 1/k_t, 1/k_s` | **devuelve radios en mm, no curvaturas** |

---

## 7. Aclaración importante sobre el notebook del BFTE

Conviene no confundirse con el título de este documento, porque el código no hace exactamente lo
que sugiere la pregunta.

En `KeratoconusAnalysis.ipynb` (el del elipsoide tórico), la llamada es:

```python
Z, mk = fn(Sx, Sy, M, p)
Rt, Rs, ext = sag_tan_fd(Z, mk, p, R)
```

Es decir: **a la curvatura se le pasa `Z`, la reconstrucción cruda, y el radio nominal
$R = 7.78$ mm.** No se le pasa el residuo sobre el BFTE ni el radio ajustado. El elipsoide
tórico se usa solo para los mapas de *elevación*; los mapas de curvatura de ese notebook se
calculan sobre la reconstrucción cruda más la esfera nominal.

En `KeratoconusAnalysisBFS.ipynb` sí se hizo la versión coherente: la referencia es la **esfera
ajustada a los datos**, con su radio $R_s$ y su ápice desplazado $(x_0,y_0)$, y lo que se deriva
numéricamente es el residuo sobre ella:

```python
u, v = x - x0, y - y0
q = np.sqrt(np.maximum(Rs**2 - u**2 - v**2, 1e-12))
dh = -np.where(m, np.nan_to_num(E), 0.0)/1000.0
dx, dy, dxx, dyy, dxy = _fd(dh, m, p)
zx  = rama*(-u/q)                 + dx
zxx = rama*(-(Rs**2 - v**2)/q**3) + dxx
```

y devuelve curvaturas en dioptrías en vez de radios en mm.

La diferencia numérica no es menor: con la esfera nominal la curvatura tangencial sale en el
rango **30.8 – 46.1 D**, y con la esfera ajustada **40.9 – 54.0 D**. El segundo rango es el
fisiológicamente creíble para esta córnea (queratometría media 45.8 D).

---

## 8. Limitaciones que deben declararse

1. **La curvatura no aporta información nueva.** Es una transformación de la elevación. Si la
   elevación está mal —por ejemplo por una referencia mal escalada— la curvatura hereda el error
   amplificado, no lo corrige.

2. **El ruido se amplifica ×4.7 D por µm.** Con 0.91 µm de ruido de reconstrucción, la
   incertidumbre punto a punto es de unas 4.7 D. Los mapas deben leerse como tendencias
   espaciales, no como valores nodo a nodo. Cualquier cifra que se cite debería promediarse
   sobre una zona.

3. **El borde es poco fiable.** Los nodos del anillo exterior usan fórmulas descentradas, que
   tienen más ruido, y algunos quedan en `NaN`. Por eso en los paneles se excluye ese anillo.

4. **La curvatura sagital usa $r$ medido desde el centro de la malla.** Si el vértice corneal
   real no coincide con ese centro, el mapa sagital queda sesgado — el tangencial, al ser local,
   es menos sensible a ello.

5. **$p = 0.4$ mm se supone constante.** El espaciado medido de la referencia en este
   experimento es 0.409 ± 0.014 mm. Un 2 % de error en $p$ es un 4 % en la derivada segunda, o
   sea cerca de 2 D sobre 45 D.
