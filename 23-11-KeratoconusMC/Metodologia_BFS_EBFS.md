# Superficies de referencia y ajuste por mínimos cuadrados

## Reconstrucción de topografía corneal a partir de patrones de Hartmann

**Experimento:** `23-11-KeratoconusMC` · paciente `maconst01` (OD) · referencia `hpref.tif`
**Notebook:** `codefiles/KeratoconusAnalysis.ipynb`, celdas 110 (BFS / BFTE) y 112 (EBFS)

---

## 1. Planteamiento: hay dos referencias encadenadas

Es la distinción que ordena todo lo demás, y confundirlas es la causa de que el cono
parezca invisible.

### Referencia 1 — física, experimental

El patrón de Hartmann de la superficie de calibración, `hpref.tif`, procesado a
`data/hpref-centroids-2.txt`. Todo el análisis arranca del campo de desplazamiento entre
los dos patrones de puntos:

```python
Dx = mat_R_X - mat_R_X_hpref      # paciente − referencia  [mm]
Dy = mat_R_Y - mat_R_Y_hpref
```

De ahí salen las pendientes `s = D/R₀`, que los tres métodos zonales (Southwell, Campbell,
Malacara) integran para obtener una elevación **diferencial** `Δh(x,y)`.

**El mapa crudo de `Δh` no sirve para diagnosticar.** Sobre estos datos vale 361 µm de pico
a valle, de los cuales 341 µm son inclinación pura: las dos capturas están descentradas
0.393 mm (≈ 20 px) una respecto de la otra. Integrar una pendiente constante de
0.393/7.78 = 0.0505 sobre una apertura de 7.37 mm produce 0.372 mm de rampa aunque la
córnea fuera una esfera perfecta. El cono, que vale ~12 µm, queda en el 3 % de la escala
de color.

### Referencia 2 — matemática, ajustada a posteriori

Una superficie geométrica ajustada por mínimos cuadrados a la córnea ya reconstruida. Es
lo que se resta para obtener el mapa de elevación clínico. Tres variantes:

| superficie                       | parámetros | qué absorbe                                  |
| -------------------------------- | ---------- | -------------------------------------------- |
| **BFS** — esfera de mejor ajuste | 4          | queratometría media y alineación             |
| **BFTE** — elipsoide tórico      | 6–8        | además el astigmatismo regular               |
| **EBFS** — esfera mejorada       | 4          | ídem BFS, pero excluyendo el cono del ajuste |

---

## 2. Paso 0 — Reconstruir la superficie corneal absoluta

Los métodos zonales devuelven `Δh`, no la córnea. Un ajuste geométrico necesita la
superficie real, así que primero se recompone sumando el casquete de referencia.

### El casquete de referencia es una ESFERA, no un toroide

$$z_{\text{ref}}(x,y) \;=\; R_0 \;-\; \sqrt{R_0^{\,2} - x^2 - y^2}, \qquad x^2+y^2 < R_0^{\,2}$$

con $R_0 = 7.78$ mm. Es la rama inferior de

$$x^2 + y^2 + (z-R_0)^2 = R_0^{\,2}$$

es decir, una esfera de radio $R_0$ centrada en $(0,0,R_0)$, con el vértice en el origen.
El _sag_ depende únicamente de $r = \sqrt{x^2+y^2}$ — sin dependencia azimutal — y por eso
es esférica. Vale 0 en el vértice y 0.9281 mm en el borde de la apertura ($r = 3.685$ mm).

El toroide **nunca** se supone: aparece solo como superficie _ajustada_ en el BFTE (§6).

### La superficie absoluta

$$z_{\text{córnea}}(x,y) \;=\; z_{\text{ref}}(x,y) \;-\; \sigma\,\Delta h(x,y), \qquad \sigma = +1$$

El signo $\sigma$ codifica la convención: $\Delta h > 0$ significa córnea más anterior, o sea
**menos** _sag_. Si el defecto apareciera hundido en vez de protuberante, esta constante
es la que hay que invertir.

```python
def superficie_absoluta(Z_m, zona):
    """sag de la cornea = casquete de referencia - elevacion diferencial."""
    xs, ys = Xc[zona], Yc[zona]
    sag_ref = R0 - np.sqrt(np.maximum(R0**2 - xs**2 - ys**2, 1e-12))
    return xs, ys, sag_ref - SIGNO * Z_m[zona]
```

> **Por qué importa este paso.** Sobre esta apertura $r/R_0 = 0.47$. El _sag_ de una esfera
> y el de un paraboloide equivalente difieren **55 µm** en el borde — cuatro veces el rasgo
> buscado. Ajustar una esfera _como esfera_, y no como parábola de segundo orden, no es un
> refinamiento cosmético.

---

## 3. Mínimos cuadrados lineales: el método matricial

Todo el ajuste algebraico de §5 se reduce a resolver un sistema **sobredeterminado**:
$N = 253$ ecuaciones, $p = 4$ incógnitas.

$$\mathbf{A}\,\mathbf{u} \;=\; \mathbf{f}, \qquad \mathbf{A}\in\mathbb{R}^{253\times 4},\;\; \mathbf{u}\in\mathbb{R}^{4}$$

No tiene solución exacta. Se busca la que minimiza la norma del residuo:

$$\hat{\mathbf{u}} \;=\; \arg\min_{\mathbf{u}} \;\lVert \mathbf{A}\mathbf{u}-\mathbf{f}\rVert_2^{2}$$

### 3.1 Ecuaciones normales — la vía que NO se usa

Derivando e igualando a cero:

$$
\mathbf{A}^{\mathsf T}\mathbf{A}\,\hat{\mathbf{u}} \;=\; \mathbf{A}^{\mathsf T}\mathbf{f}
\qquad\Longrightarrow\qquad
\hat{\mathbf{u}} = (\mathbf{A}^{\mathsf T}\mathbf{A})^{-1}\mathbf{A}^{\mathsf T}\mathbf{f}
$$

Es correcto en aritmética exacta y desastroso en punto flotante: el número de condición se
**eleva al cuadrado**,

$$\kappa(\mathbf{A}^{\mathsf T}\mathbf{A}) = \kappa(\mathbf{A})^2$$

lo que duplica la pérdida de dígitos significativos. Medido sobre la matriz real de este
ajuste:

```
forma de A: (253, 4)    rango: 4
valores singulares: [58.492  58.359  22.302   6.388]
cond(A)     = s_max/s_min = 9.16
cond(A^T A) = cond(A)^2   = 83.85
```

Aquí $\kappa = 9.16$ es benigno y las tres vías coinciden a $10^{-15}$, pero la práctica
correcta es no formar nunca $\mathbf{A}^{\mathsf T}\mathbf{A}$.

### 3.2 Descomposición QR

$\mathbf{A} = \mathbf{Q}\mathbf{R}$ con $\mathbf{Q}$ ortogonal ($\mathbf{Q}^{\mathsf T}\mathbf{Q}=\mathbf{I}$)
y $\mathbf{R}$ triangular superior. Como las transformaciones ortogonales preservan la
norma euclídea:

$$\lVert \mathbf{A}\mathbf{u}-\mathbf{f}\rVert = \lVert \mathbf{Q}\mathbf{R}\mathbf{u}-\mathbf{f}\rVert = \lVert \mathbf{R}\mathbf{u}-\mathbf{Q}^{\mathsf T}\mathbf{f}\rVert$$

y el problema se resuelve por sustitución hacia atrás sobre $\mathbf{R}\mathbf{u} = \mathbf{Q}^{\mathsf T}\mathbf{f}$,
conservando $\kappa(\mathbf{A})$ en vez de elevarlo al cuadrado.

### 3.3 SVD y pseudoinversa — lo que hace `np.linalg.lstsq`

$$
\mathbf{A} = \mathbf{U}\boldsymbol\Sigma\mathbf{V}^{\mathsf T}, \qquad
\hat{\mathbf{u}} = \mathbf{V}\boldsymbol\Sigma^{+}\mathbf{U}^{\mathsf T}\mathbf{f} = \mathbf{A}^{+}\mathbf{f}
$$

donde $\boldsymbol\Sigma^{+}$ invierte los valores singulares no nulos y deja en cero los que
caen por debajo de `rcond`. Es la vía más robusta: funciona aunque $\mathbf{A}$ sea
deficiente de rango, y entonces devuelve la solución de **norma mínima**.

`numpy.linalg.lstsq` usa la rutina LAPACK `gelsd`, que es SVD con _divide and conquer_.
Verificación numérica sobre este ajuste:

```
solucion por SVD        vs lstsq : diferencia max 5.33e-15
solucion por ec. normales vs lstsq : diferencia max 7.99e-15
```

### 3.4 Interpretación geométrica

$\mathbf{A}\hat{\mathbf{u}}$ es la **proyección ortogonal** de $\mathbf{f}$ sobre el subespacio
generado por las columnas de $\mathbf{A}$. El residuo $\mathbf{f}-\mathbf{A}\hat{\mathbf{u}}$
es perpendicular a ese subespacio — que es exactamente lo que dicen las ecuaciones normales
$\mathbf{A}^{\mathsf T}(\mathbf{f}-\mathbf{A}\hat{\mathbf{u}})=\mathbf{0}$.

---

## 4. Mínimos cuadrados no lineales

Cuando el modelo no es lineal en sus parámetros (§5.2, §6) se minimiza

$$
S(\boldsymbol\theta) = \tfrac12\lVert \mathbf{r}(\boldsymbol\theta)\rVert_2^2, \qquad
r_i(\boldsymbol\theta) = \text{modelo}(\boldsymbol\theta; x_i,y_i) - z_i
$$

### 4.1 Gauss-Newton

Se linealiza el residuo alrededor del punto actual mediante el Jacobiano
$J_{ij} = \partial r_i/\partial \theta_j$:

$$\mathbf{r}(\boldsymbol\theta+\boldsymbol\delta) \approx \mathbf{r}(\boldsymbol\theta) + \mathbf{J}\boldsymbol\delta$$

lo que convierte cada iteración en un problema lineal:

$$\mathbf{J}^{\mathsf T}\mathbf{J}\,\boldsymbol\delta = -\,\mathbf{J}^{\mathsf T}\mathbf{r}$$

### 4.2 Levenberg-Marquardt

Gauss-Newton diverge si $\mathbf{J}^{\mathsf T}\mathbf{J}$ está mal condicionada.
Levenberg-Marquardt añade amortiguación:

$$\left(\mathbf{J}^{\mathsf T}\mathbf{J} + \lambda\,\operatorname{diag}(\mathbf{J}^{\mathsf T}\mathbf{J})\right)\boldsymbol\delta = -\,\mathbf{J}^{\mathsf T}\mathbf{r}$$

con $\lambda \to 0$ recuperando Gauss-Newton (convergencia cuadrática cerca del óptimo) y
$\lambda \to \infty$ degenerando en descenso de gradiente (robusto lejos). `scipy.optimize.least_squares`
usa por defecto `trf` (_trust region reflective_), que es la variante con región de confianza
y admite cotas — necesarias en §6 para mantener los radios en rango físico.

### 4.3 Incertidumbre de los parámetros

De la inversa del Hessiano aproximado:

$$
\operatorname{Cov}(\hat{\boldsymbol\theta}) \approx \hat\sigma^2\,(\mathbf{J}^{\mathsf T}\mathbf{J})^{-1},
\qquad \hat\sigma^2 = \frac{\lVert\mathbf{r}\rVert^2}{N-p}
$$

Las raíces de la diagonal son los errores estándar. Medido para el BFS:

```
Jacobiano J: forma (253, 4)    cond(J) = 49.3
error estandar de R: 0.0141 mm
```

---

## 5. Best Fit Sphere (BFS)

### 5.1 Linealización algebraica de Kåsa

Una esfera de centro $(x_0,y_0,z_0)$ y radio $R$ cumple

$$(x-x_0)^2 + (y-y_0)^2 + (z-z_0)^2 = R^2$$

El truco: al desarrollar los cuadrados, los términos $x_0^2, y_0^2, z_0^2, R^2$ aparecen
**solo agrupados**, y se absorben en una única incógnita auxiliar

$$c \;=\; R^2 - x_0^2 - y_0^2 - z_0^2$$

quedando

$$\boxed{\;x^2+y^2+z^2 \;=\; 2x_0\,x \;+\; 2y_0\,y \;+\; 2z_0\,z \;+\; c\;}$$

que es **lineal** en $(x_0,y_0,z_0,c)$. Con $f_i = x_i^2+y_i^2+z_i^2$:

$$
\mathbf{A} = \begin{bmatrix}
2x_1 & 2y_1 & 2z_1 & 1\\
2x_2 & 2y_2 & 2z_2 & 1\\
\vdots & \vdots & \vdots & \vdots\\
2x_N & 2y_N & 2z_N & 1
\end{bmatrix},
\qquad
\mathbf{u} = \begin{bmatrix} x_0\\ y_0\\ z_0\\ c\end{bmatrix},
\qquad
\mathbf{f} = \begin{bmatrix} x_1^2+y_1^2+z_1^2\\ \vdots\\ x_N^2+y_N^2+z_N^2\end{bmatrix}
$$

y el radio se recupera como

$$R = \sqrt{c + x_0^2+y_0^2+z_0^2}$$

**Ventaja decisiva:** solución cerrada, sin valor inicial y sin mínimos locales.
**Limitación:** minimiza el residuo _algebraico_ $\lVert\mathbf{A}\mathbf{u}-\mathbf{f}\rVert$,
no la distancia geométrica, lo que pondera de más los puntos alejados del centro.

### 5.2 Refinamiento geométrico

La definición estricta de esfera de mejor ajuste minimiza la distancia ortogonal real:

$$\min_{x_0,y_0,z_0,R} \;\sum_{i=1}^{N} \Big( \underbrace{\sqrt{(x_i-x_0)^2+(y_i-y_0)^2+(z_i-z_0)^2}}_{d_i} \;-\; R \Big)^{2}$$

Es no lineal, pero arrancando de Kåsa converge en **3 iteraciones**:

```
R:  7.389341  ->  7.395614 mm      (cambio 6.3 um)
coste 0.5*||r||^2 : 1.027182e-02
```

### 5.3 Código

```python
def ajustar_esfera(x, y, z):
    """Best Fit Sphere minimizando la distancia ortogonal.

    Arranque algebraico (Kasa): para puntos sobre una esfera
        x^2+y^2+z^2 = 2*x0*x + 2*y0*y + 2*z0*z + (Rs^2 - x0^2 - y0^2 - z0^2)
    es LINEAL en (x0, y0, z0, c), asi que sale de un lstsq sin necesidad de
    valor inicial ni de riesgo de minimo local.  Despues se refina con la
    distancia geometrica real, que es la definicion estricta de mejor ajuste.
    """
    A = np.column_stack([2*x, 2*y, 2*z, np.ones_like(x)])
    q, *_ = np.linalg.lstsq(A, x**2 + y**2 + z**2, rcond=None)
    x0, y0, z0 = q[:3]
    Rs = np.sqrt(max(q[3] + x0**2 + y0**2 + z0**2, 1e-9))
    sol = least_squares(
        lambda w: np.sqrt((x-w[0])**2 + (y-w[1])**2 + (z-w[2])**2) - w[3],
        [x0, y0, z0, Rs])
    return sol.x


def sag_esfera(par, x, y, z_ref):
    """Sag de la esfera ajustada, eligiendo la rama que corresponde al casquete."""
    x0, y0, z0, Rs = par
    raiz = np.sqrt(np.maximum(Rs**2 - (x-x0)**2 - (y-y0)**2, 1e-12))
    a, b = z0 - raiz, z0 + raiz
    return a if np.std(z_ref - a) < np.std(z_ref - b) else b
```

### 5.4 El mapa de elevación

Un mapa de elevación es una diferencia **axial**, no una distancia perpendicular. Se evalúa
el _sag_ de la esfera ajustada en el mismo $(x,y)$ de cada nodo:

$$z_{\text{BFS}}(x,y) = z_0 \pm \sqrt{R^2-(x-x_0)^2-(y-y_0)^2}$$

(la rama se elige por la que mejor sigue al casquete, comparando desviaciones típicas) y

$$E(x,y) \;=\; \big[\, z_{\text{BFS}}(x,y) - z_{\text{córnea}}(x,y) \,\big]\times 10^{3} \quad [\mu\text{m}]$$

Positivo = la córnea está por delante de la esfera = protuberancia. Es la convención clínica.

```python
def mapa_elevacion(par_sag, xs, ys, zs, zona):
    """Elevacion axial [um] sobre toda la zona, positiva hacia delante."""
    out = np.full(M.shape, np.nan)
    out[zona] = (par_sag(xs, ys, zs) - zs) * 1000
    return out
```

### 5.5 Resultados

| método    | $R_{\text{BFS}}$ | K       | elevación rms | pico-valle |
| --------- | ---------------- | ------- | ------------- | ---------- |
| Southwell | 7.3956 mm        | 45.64 D | 10.05 µm      | 64.7 µm    |
| Campbell  | 7.3969 mm        | 45.63 D | 9.52 µm       | 66.6 µm    |
| Malacara  | 7.3949 mm        | 45.64 D | 8.95 µm       | 65.7 µm    |

Los tres algoritmos de reconstrucción —independientes entre sí— coinciden hasta la cuarta
cifra decimal.

---

## 6. Best Fit Toric Ellipsoid (BFTE)

### 6.1 La superficie bicónica

$$z(x,y) \;=\; \frac{c_x x'^2 + c_y y'^2}{1+\sqrt{1-(1+Q_x)c_x^2x'^2-(1+Q_y)c_y^2y'^2}}$$

con los ejes girados un ángulo $\varphi$,

$$
x' = (x-x_0)\cos\varphi + (y-y_0)\sin\varphi, \qquad
y' = -(x-x_0)\sin\varphi + (y-y_0)\cos\varphi
$$

y $c_x = 1/R_x$, $c_y = 1/R_y$ las curvaturas apicales de los dos meridianos principales.
$Q$ es la asfericidad: $Q=0$ esfera, $Q<0$ elipse prolata, $Q>0$ oblata.

**Por qué esta parametrización y no los semiejes $a,b,c$ del elipsoide:** sobre un casquete
poco profundo los semiejes quedan casi indeterminados, mientras que los radios apicales
están bien definidos. Es además la parametrización estándar en óptica oftálmica.

El astigmatismo corneal sale directo:

$$\text{astig [D]} = \left|\frac{337.5}{R_x} - \frac{337.5}{R_y}\right|$$

### 6.2 Degeneración entre asfericidad y toricidad

Con $Q_x, Q_y$ libres el ajuste tiene 8 parámetros y se vuelve degenerado sobre 7.4 mm:
asfericidad y toricidad se compensan entre sí. Medido con bootstrap de 20 remuestreos:

| modo      | nº par. | astigmatismo      | rms elev. | $\kappa(\mathbf{J})$ |
| --------- | ------- | ----------------- | --------- | -------------------- |
| Q libre   | 8       | 0.55 ± 0.19 D     | 4.04 µm   | 2.6·10²              |
| Q común   | 7       | 1.84 ± 0.07 D     | 4.42 µm   | 2.0·10²              |
| **Q = 0** | **6**   | **1.88 ± 0.08 D** | 4.53 µm   | **5.5·10¹**          |

Las asfericidades libres compran 0.5 µm de ajuste sobre 4.5 y a cambio triplican el
astigmatismo estimado, duplican su incertidumbre y hacen que el ajuste de Malacara converja
al meridiano perpendicular. El **eje**, en cambio, sale 72–73° en los tres modos con
concentración bootstrap 1.00: lo mal determinado es la magnitud, no la orientación.

Pentacam sí ajusta asfericidad, pero sobre zona de 8–9 mm y con mucha mayor densidad de
puntos; ahí la degeneración se rompe. Con 253 nodos sobre 7.37 mm, no. **Por eso el valor
por defecto es `MODO_TORO = 'Q_cero'`.**

### 6.3 Código

```python
def sag_toroide(par, x, y, modo=None):
    """Sag de una superficie biconica (elipsoide torico) con el eje girado.

        sag = (cx*xp^2 + cy*yp^2)
              / (1 + sqrt(1 - (1+Qx)*cx^2*xp^2 - (1+Qy)*cy^2*yp^2))

    Rx, Ry son los radios apicales de los dos meridianos principales y Qx, Qy
    sus asfericidades (Q=0 esfera, Q<0 elipse prolata, Q>0 oblata).  Es la
    parametrizacion estandar en optica oftalmica y esta mucho mejor
    condicionada que dar los semiejes a, b, c del elipsoide: sobre un casquete
    poco profundo los semiejes quedan casi indeterminados, mientras que los
    radios apicales estan bien definidos.
    """
    if modo is not None:
        par = _desempaquetar(par, modo)
    x0, y0, z0, Rx, Ry, Qx, Qy, phi = par
    dx, dy = x - x0, y - y0
    xp = dx*np.cos(phi) + dy*np.sin(phi)
    yp = -dx*np.sin(phi) + dy*np.cos(phi)
    cx, cy = 1.0/Rx, 1.0/Ry
    num = cx*xp**2 + cy*yp**2
    rad = 1 - (1+Qx)*cx**2*xp**2 - (1+Qy)*cy**2*yp**2
    return z0 + num / (1 + np.sqrt(np.maximum(rad, 1e-9)))

def ajustar_toroide(x, y, z, modo=None):
    """Best Fit Toric Ellipsoid. Arranca de la esfera de referencia (Q=0).

    Devuelve los 8 parametros desempaquetados, con Rx = meridiano MAS CURVO.
    """
    modo = MODO_TORO if modo is None else modo
    p0, bnd = _ARRANQUE[modo]
    sol = least_squares(lambda w: sag_toroide(w, x, y, modo) - z, p0, bounds=bnd)
    par = list(_desempaquetar(sol.x, modo))
    if par[3] > par[4]:                 # Rx = meridiano MAS CURVO (radio menor)
        par[3], par[4] = par[4], par[3]
        par[5], par[6] = par[6], par[5]
        par[7] += np.pi/2
    return np.array(par)
```

### 6.4 Resultados

| método    | $R$ curvo           | $R$ plano           | astigmatismo | eje   | elev. rms |
| --------- | ------------------- | ------------------- | ------------ | ----- | --------- |
| Southwell | 7.2461 mm (46.58 D) | 7.5521 mm (44.69 D) | 1.89 D       | 72.8° | 4.53 µm   |
| Campbell  | 7.2520 mm (46.54 D) | 7.5482 mm (44.71 D) | 1.83 D       | 71.7° | 4.39 µm   |
| Malacara  | 7.2675 mm (46.44 D) | 7.5276 mm (44.83 D) | 1.60 D       | 82.7° | 4.75 µm   |

---

## 7. Enhanced Best Fit Sphere (EBFS)

### 7.1 El problema que resuelve

La BFS estándar se ajusta a **todos** los puntos, el cono incluido. El cono tira del ajuste
hacia sí mismo y se auto-oculta: cuanto más pronunciado, más se lleva la esfera consigo.

### 7.2 El procedimiento (Belin / Ambrósio)

1. Ajustar la BFS estándar sobre toda la zona (§5).
2. Localizar el punto de referencia y **excluir un disco de 3.5 mm de diámetro** centrado
   en él.
3. Reajustar la esfera **solo con los puntos de fuera** — la periferia sana.
4. Calcular la elevación de **toda** la apertura contra esa esfera mejorada.
5. El índice diagnóstico es la diferencia de elevación en el punto de referencia entre el
   mapa EBFS y el BFS.

La clave está en el paso 3: la referencia deja de contener el cono, así que el cono
sobresale de verdad.

### 7.3 Código

```python
def enhanced(Z_m, zona, centro, modelo=None, diam=None):
    """Reajusta excluyendo un disco de diametro `diam` centrado en `centro`."""
    modelo = MODELO_EBFS if modelo is None else modelo
    diam = DIAM_EXCLUSION if diam is None else diam
    xs, ys, zs = superficie_absoluta(Z_m, zona)
    fuera = np.hypot(xs - centro[0], ys - centro[1]) > diam/2.0
    if fuera.sum() < 30:
        raise ValueError("quedan %d puntos fuera de la exclusion: insuficiente" % fuera.sum())
    par, sag = _ajustar(xs[fuera], ys[fuera], zs[fuera], modelo)
    return par, mapa_elevacion(sag, xs, ys, zs, zona), fuera.sum()
```

La localización del ápice se hace sobre un promedio local 3×3 restringido a nodos
interiores. El cono mide ~2 mm, o sea varios nodos; buscar el máximo nodo a nodo devolvería
cualquier punto ruidoso del anillo exterior, que solo tiene dos vecinos:

```python
def apice_cono(residuo, mascara):
    """Devuelve (i, j) del apice del defecto, o None si no se puede determinar.

    Busca el maximo de un promedio local 3x3 sobre los nodos interiores.  El
    rasgo mide ~2 mm, o sea varios nodos: buscar el maximo nodo a nodo
    devolveria cualquier punto ruidoso del anillo exterior.
    """
    interior = mascara_interior(mascara)
    suave = np.full_like(residuo, np.nan)
    nr_p, nc_p = residuo.shape
    for i in range(1, nr_p - 1):
        for j in range(1, nc_p - 1):
            if not interior[i, j]:
                continue
            bloque = residuo[i-1:i+2, j-1:j+2]
            if np.isfinite(bloque).sum() >= 6:
                suave[i, j] = np.nanmean(bloque)
    if not np.isfinite(suave).any():
        return None
    return np.unravel_index(np.nanargmax(suave), suave.shape)
```

### 7.4 Resultados

```
                     centro = max elevacion BFS      centro = apice irregular
Southwell                 +11.4 um                        -1.3 um
Campbell                  +12.1 um                        +2.2 um
Malacara                   +9.7 um                        -0.9 um
```

---

## 8. Limitaciones que deben declararse

### 8.1 El radio de referencia $R_0$ es un factor de escala, no una medida

$R_0 = 7.78$ mm es una constante escrita a mano. En el notebook hace **tres** trabajos a la
vez: convierte desplazamiento en pendiente (`s = D/R₀`), construye el casquete de referencia
y sirve de línea base en `c = 1/R + S`. Los comentarios del propio código la describen unas
veces como "radio de referencia" y otras como "distancia al plano de detección", que son
magnitudes físicas distintas.

Como $z_{\text{ref}} \approx r^2/2R_0$ y $\Delta h \propto 1/R_0$, **toda la superficie
absoluta escala como $1/R_0$** y el radio ajustado sale proporcional a $R_0$:

| $R_0$ supuesto | $R_{\text{BFS}}$ | K           | $R_{\text{BFS}}/R_0$ |
| -------------- | ---------------- | ----------- | -------------------- |
| 7.00           | 6.6738           | 50.57 D     | 0.9534               |
| 7.40           | 7.0437           | 47.92 D     | 0.9518               |
| **7.78**       | **7.3956**       | **45.64 D** | **0.9506**           |
| 8.20           | 7.7851           | 43.35 D     | 0.9494               |
| 8.60           | 8.1564           | 41.38 D     | 0.9484               |

El cociente es constante. **K = 45.64 D no es una medida independiente**: es
$337.5/(0.951\,R_0)$. Lo que sí se mide es el cociente: la córnea es un 4.9 % más curva que
la superficie de referencia.

La forma del mapa, en cambio, es robusta:

```
correlacion del mapa de elevacion contra el caso R0 = 7.78
   R0=7.00  ->  r = +0.99626    diferencia rms = 1.40 um
   R0=7.40  ->  r = +0.99929    diferencia rms = 0.62 um
   R0=8.20  ->  r = +0.99941    diferencia rms = 0.58 um
   R0=8.60  ->  r = +0.99812    diferencia rms = 1.06 um
```

Un 10 % de error en $R_0$ mueve el mapa menos de 1.4 µm rms, frente a los ~12 µm del cono.

> **Formulación sugerida.** _El radio de referencia $R_0$ actúa como factor de escala global
> de la reconstrucción; la queratometría absoluta hereda su incertidumbre, mientras que la
> topografía de elevación relativa —y por tanto la detección y localización del cono— es
> insensible a él dentro de ±10 %._

### 8.2 El centrado del EBFS es circular

Belin/Ambrósio centra la exclusión en el punto más delgado medido por **paquimetría**, un
dato independiente del mapa de elevación. Un montaje de Hartmann no mide espesor, así que
el centro hay que derivarlo del propio mapa. Los tres métodos de reconstrucción coinciden
dentro de ±1.2 µm, pero el criterio de centrado mueve el índice de ~+11 µm a ~0 µm.

### 8.3 El índice depende de un parámetro libre

```
exclusion 2.5 mm  ->  +6.9 um
exclusion 3.0 mm  ->  +9.3 um
exclusion 3.5 mm  -> +11.4 um
exclusion 4.0 mm  -> +13.5 um
exclusion 4.5 mm  -> +17.5 um
```

Los umbrales publicados (>12 µm sospechoso, >15 µm patológico) se cruzan ambos moviendo solo
ese parámetro, sin tocar el dato. Además se definen sobre zona de 8 mm y esta apertura mide
7.37 mm.

### 8.4 La geometría de exclusión no es la prevista

Con el criterio del máximo de elevación, el disco queda centrado a 3.3 mm del eje sobre una
apertura de radio 3.69 mm — pegado al borde. Belin/Ambrósio supone un cono central o
paracentral. Por eso $R_{\text{EBFS}} = 7.3611$ mm sale más curvo que $R_{\text{BFS}} = 7.3956$ mm,
que es lo contrario de lo esperable al excluir un cono.

### 8.5 Ruido del anillo exterior

Los nodos del borde tienen solo dos vecinos, su elevación se apoya en muy pocas ecuaciones y
son con diferencia los más ruidosos: llegan a −44 µm frente a los ~12 µm del rasgo real. Si
se dejan dentro, ellos solos fijan la escala de color. Por eso se excluyen del cálculo del
límite de color y de la búsqueda del ápice (`mascara_interior`, `escala_robusta`).

---

## 9. Conclusión

La **EBFS funciona bien como herramienta de visualización** — el cono se distingue con
claridad contra esa referencia — pero el índice numérico no es citable como diagnóstico
graduado con estos datos.

Formulación defendible:

> Se implementó la metodología Belin/Ambrósio de esfera de mejor ajuste mejorada. Los tres
> algoritmos de reconstrucción zonal producen índices consistentes entre sí (±1.2 µm), pero
> el valor absoluto depende del criterio de centrado en ausencia de paquimetría, por lo que
> se reporta como indicador cualitativo de la presencia del cono y no como índice diagnóstico
> graduado.

El resultado sólido del trabajo es la **validación cruzada**: tres algoritmos zonales
independientes —iterativo con ecuaciones normales y SOR (Southwell), una sola pasada por
cuadrantes con media móvil (Campbell), e integración de camino (Malacara)— reconstruyen la
misma superficie con correlación $r = 0.999$, y localizan el mismo cono dentro de un nodo de
malla (0.4 mm) y 1.6 µm de altura. Esa dispersión entre métodos constituye una estimación
empírica de la incertidumbre de la reconstrucción.

---

## Apéndice — Constantes del análisis

| símbolo  | código           | valor      | significado                                             |
| -------- | ---------------- | ---------- | ------------------------------------------------------- |
| $R_0$    | `R0`             | 7.78 mm    | radio del casquete de referencia (§8.1)                 |
| $p$      | `p`              | 0.4 mm     | espaciado nominal de la rejilla (medido: 0.409 ± 0.014) |
| $\delta$ | `delta`          | 0.02 mm/px | escala de la imagen                                     |
| $\sigma$ | `SIGNO`          | +1         | convención de signo de $\Delta h$                       |
| —        | `DIAM_EXCLUSION` | 3.5 mm     | disco excluido en el EBFS                               |
| —        | `MODO_TORO`      | `'Q_cero'` | asfericidades fijadas a cero (§6.2)                     |
| —        | `ZONA_MM`        | `None`     | zona de ajuste; `None` = toda la apertura (7.37 mm)     |
