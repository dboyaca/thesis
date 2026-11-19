# Descomposición de Zernike de la superficie reconstruida

## Qué es, para qué se usa aquí, y qué no es

**Experimento:** `21-09-26-VL-Analysis` · ojo derecho `ojoderechosinmirilla2907.jpg` (1024×1024)
**Notebook:** `codefiles/KeratoconusAnalysis.ipynb` · funciones en la celda 87, uso en la celda 108

---

## 1. Qué es Zernike, en una frase

Es un **cambio de base**: igual que Fourier, pero sobre un disco en vez de sobre una línea.

Fourier dice *«esta señal tiene tanta energía a 50 Hz, tanta a 100 Hz»*. Zernike dice *«esta
superficie tiene tanto desenfoque, tanto astigmatismo, tanta coma»*. En los dos casos el dato no
cambia: se reescribe en una base donde cada componente tiene nombre y significado físico.

### Los índices

Los polinomios se indexan con dos enteros, $Z_n^m(\rho,\theta)$:

- **$n$ — orden radial.** Cuánto detalle radial. $n$ bajo = formas suaves y globales;
  $n$ alto = estructura fina.
- **$m$ — frecuencia angular.** Cuántos lóbulos aparecen al girar. $m=0$ es rotacionalmente
  simétrico. Se cumple $|m|\le n$ y $n-|m|$ par.

La forma explícita es un polinomio radial multiplicado por una sinusoidal azimutal:

$$Z_n^m(\rho,\theta) = N_n^m \; R_n^{|m|}(\rho)\;\begin{cases}\cos(|m|\theta) & m\ge 0\\ \sin(|m|\theta) & m<0\end{cases}$$

$$R_n^{|m|}(\rho) = \sum_{k=0}^{(n-|m|)/2} \frac{(-1)^k\,(n-k)!}{k!\left(\frac{n+|m|}{2}-k\right)!\left(\frac{n-|m|}{2}-k\right)!}\;\rho^{\,n-2k}$$

### Los primeros modos son las aberraciones clásicas

| $n$ | modos | nombre óptico |
|---|---|---|
| 0 | $Z_0^0$ | pistón (constante) |
| 1 | $Z_1^{\pm1}$ | inclinación (*tilt*) |
| 2 | $Z_2^0$, $Z_2^{\pm2}$ | desenfoque, astigmatismo |
| 3 | $Z_3^{\pm1}$, $Z_3^{\pm3}$ | **coma**, trefoil |
| 4 | $Z_4^0$ | aberración esférica |

### Por qué esta base y no otra

Son **ortogonales sobre el disco unitario**, que es la forma natural de una pupila. Ortogonal
significa que los coeficientes son independientes: puedes leer *«hay 8.96 µm de astigmatismo»*
sin que ese número dependa de cuántos otros modos incluyeras en el ajuste. Un polinomio
cualquiera no tiene esa propiedad, y entonces los coeficientes se contaminan entre sí.

### La normalización importa

Se usa la normalización de Noll,

$$N_n^m = \sqrt{\frac{2(n+1)}{1+\delta_{m0}}}$$

con la que **cada polinomio tiene RMS unitario sobre el disco**. La consecuencia práctica es que
el coeficiente en µm *es directamente la contribución RMS de ese modo*, y por tanto se pueden
comparar entre sí y sumar en cuadratura:

$$\text{RMS total} = \sqrt{\sum_{i\neq\text{pistón}} c_i^2}$$

Comprobación sobre estos datos: la suma en cuadratura da **11.626 µm** y el RMS del mapa medido
es **11.534 µm**. Cierra al 0.8 %.

---

## 2. Para qué se usa aquí

El problema original era: tengo un mapa de elevación de decenas de micras y un defecto corneal
que debería valer unas pocas. ¿Dónde está el resto?

Se podría adivinar. Zernike lo contesta con un **presupuesto de energía cerrado**, y ese es todo
su papel en este trabajo: **es el instrumento de diagnóstico, no un paso del procesamiento.**

```
energia por orden radial [um rms]
   orden 1:  5.409     inclinacion  (descentramiento de captura)
   orden 2: 10.220     desenfoque + astigmatismo
   orden 3:  1.112     coma y trefoil
   orden 4:  0.399
   orden 5:  0.198
   orden 6:  0.118
   residuo no explicado: 0.313 um rms   <- ruido de medida
```

Con eso ya no hay que decidir nada a ojo: los términos que se comen el rango dinámico están
identificados y cuantificados.

---

## 3. El punto que lo une todo: Zernike **es** el ajuste

Esto es lo que suele costar ver. **Quitar el toro de mejor ajuste y quitar Zernike de orden
$n\le 2$ son literalmente la misma operación.**

Los conjuntos

$$\{1,\; x,\; y,\; r^2,\; x^2-y^2,\; xy\} \qquad\text{y}\qquad \{Z_0^0,\; Z_1^{\pm1},\; Z_2^0,\; Z_2^{\pm2}\}$$

generan **el mismo subespacio de funciones**. Restar la proyección sobre uno o sobre el otro da
idéntico resultado. Verificado numéricamente:

```
residuo quitando el polinomio {r2, x2-y2, xy, x, y, 1} : rms 1.223818 um
residuo quitando Zernike n <= 2                        : rms 1.223818 um
diferencia maxima punto a punto                        : 3.286e-14 um
```

Esos $3\times10^{-14}$ µm son precisión de punto flotante, no una coincidencia aproximada.

**Entonces, ¿por qué usar las dos formas?**

- En el **código** se resta el polinomio, porque son cuatro líneas de `lstsq` y no hay que
  construir polinomios radiales.
- En el **informe** se usan los coeficientes de Zernike, porque ahí cada número tiene nombre y
  unidades interpretables. El coeficiente $a = -0.0042$ del paraboloide no le dice nada a nadie;
  *«astigmatismo 0/90 = +8.96 µm»* sí.

---

## 4. Lo que **no** es: esto no es lo que hace Pentacam

Conviene dejarlo claro porque es una confusión natural. En Pentacam son **dos salidas
independientes** del mismo dato, ninguna derivada de la otra:

- **BFS / BFTE** — ajuste geométrico directo de una esfera o un elipsoide tórico por mínimos
  cuadrados sobre una zona fija (8.0 o 9.0 mm). No hay Zernike de por medio. Está descrito en
  `../23-11-KeratoconusMC/Metodologia_BFS_EBFS.md`.
- **Zernike** — vive en el módulo de aberrometría corneal, como producto aparte para cuantificar
  aberraciones de alto orden.

El flujo clínico **no** calcula uno a partir del otro.

Lo que se hizo aquí es distinto, y venía de un problema que un topógrafo comercial no tiene: en
Pentacam la zona de ajuste y la alineación están fijadas y calibradas por el instrumento. Aquí no
—había descentramiento entre capturas que nadie había medido— y no se sabía si el mapa estaba
mal calculado o solo mal escalado. Zernike fue el instrumento de medida para averiguarlo.

> **Formulación para la tesis.** *El ajuste de la superficie de referencia es un mínimo cuadrado
> polinómico de segundo orden; la descomposición de Zernike se emplea para caracterizar y
> reportar el contenido aberrométrico, y es equivalente al ajuste en sus primeros seis modos.*

---

## 5. Cómo se calcula: mínimos cuadrados lineales

El ajuste es lineal en los coeficientes, así que es un `lstsq` y nada más.

### Construcción del sistema

Para los $N$ nodos válidos se normaliza el radio al máximo de la apertura,
$\rho_i = r_i/r_{\max}$, y se monta la matriz de diseño evaluando cada modo en cada nodo:

$$\mathbf{A} = \begin{bmatrix}
Z_{n_1}^{m_1}(\rho_1,\theta_1) & Z_{n_2}^{m_2}(\rho_1,\theta_1) & \cdots \\
Z_{n_1}^{m_1}(\rho_2,\theta_2) & Z_{n_2}^{m_2}(\rho_2,\theta_2) & \cdots \\
\vdots & \vdots & \ddots
\end{bmatrix} \in \mathbb{R}^{253\times 28}$$

y se resuelve el sistema sobredeterminado

$$\mathbf{A}\,\mathbf{c} = \mathbf{z}, \qquad
\hat{\mathbf{c}} = \arg\min_{\mathbf{c}}\lVert\mathbf{A}\mathbf{c}-\mathbf{z}\rVert_2^2$$

`numpy.linalg.lstsq` lo hace por **SVD** (rutina LAPACK `gelsd`), calculando la pseudoinversa
$\mathbf{A}^+ = \mathbf{V}\boldsymbol\Sigma^+\mathbf{U}^{\mathsf T}$. Es la vía robusta: funciona
aunque $\mathbf{A}$ fuera deficiente de rango y entonces devuelve la solución de norma mínima.
El tratamiento matricial detallado —ecuaciones normales, QR, SVD, condicionamiento— está en
`../23-11-KeratoconusMC/Metodologia_BFS_EBFS.md`, sección 3.

### El código

```python
def zernike(n, m, rho, theta):
    """Polinomio de Zernike Z(n,m) con normalizacion de Noll (rms unitario)."""
    import math
    Rad = np.zeros_like(rho); am = abs(m)
    for k in range((n - am)//2 + 1):
        Rad += ((-1)**k * math.factorial(n - k) /
                (math.factorial(k) * math.factorial((n + am)//2 - k) *
                 math.factorial((n - am)//2 - k))) * rho**(n - 2*k)
    norm = np.sqrt((2*(n + 1)) / (1 + (m == 0)))
    return norm * Rad * (np.cos(am*theta) if m >= 0 else np.sin(am*theta))


def descomponer_zernike(Z, mascara, X_pts, Y_pts, n_max=6):
    """Ajusta Zernike hasta orden n_max. Devuelve (modos, coeficientes [um], residuo [um])."""
    m = np.asarray(mascara, dtype=bool)
    x = np.asarray(X_pts, float)[m]; x = x - x.mean()
    y = np.asarray(Y_pts, float)[m]; y = y - y.mean()
    z = np.asarray(Z, float)[m] * 1000.0            # um
    r = np.hypot(x, y); rho = r / r.max(); th = np.arctan2(y, x)
    modos = [(n, mm) for n in range(n_max + 1) for mm in range(-n, n + 1, 2)]
    A = np.column_stack([zernike(n, mm, rho, th) for n, mm in modos])
    coef, *_ = np.linalg.lstsq(A, z, rcond=None)
    return modos, coef, z - A @ coef
```

Tres detalles del código que conviene notar:

1. `x - x.mean()` centra la apertura. Sin centrar, la inclinación y el desenfoque se mezclan.
2. `rho = r / r.max()` normaliza al disco unitario, donde los polinomios están definidos.
3. `* 1000.0` convierte mm a µm, para que los coeficientes salgan en las unidades en que se
   reportan las aberraciones.

---

## 6. Las dos validaciones que hacen creíbles los números

Esta apertura **no es un disco perfecto**: es una rejilla 19×19 recortada en diamante, 253
nodos. Zernike es ortogonal sobre el disco *continuo*, así que sobre una apertura discreta y
recortada la ortogonalidad solo es aproximada y los coeficientes podrían estar acoplados. Hay
que medirlo, no suponerlo.

### 6.1 Ortogonalidad efectiva

Se normalizan las columnas de $\mathbf{A}$ y se calcula la matriz de Gram
$\mathbf{G} = \hat{\mathbf{A}}^{\mathsf T}\hat{\mathbf{A}}$. Si los modos fueran perfectamente
ortogonales, $\mathbf{G}$ sería la identidad.

```
modos = 28   nodos = 253
producto interno cruzado maximo = 0.0839   (peor par: Z(2,2) con Z(6,6))
cond(A) = 1.31
```

Condición 1.31 es prácticamente ortonormal —un sistema bien planteado tiene condición cercana a
1— y el acoplamiento del 8 % está entre dos modos de orden alto que aquí valen casi nada.

### 6.2 Estabilidad frente al truncamiento

Es la consecuencia práctica de lo anterior, y lo que de verdad importa: si los coeficientes
cambiaran según dónde cortas la serie, no significarían nada.

```
n_max        Z(2,0)      Z(2,2)     Z(3,-1)    residuo
  2          -0.388      +8.956           -      1.224
  4          -0.380      +8.954      -0.216      0.391
  6          -0.381      +8.956      -0.218      0.313
  8          -0.388      +8.953      -0.208      0.210
```

Los coeficientes no se mueven en la tercera cifra. **Eso es lo que convierte a Zernike en un
instrumento de medida y no en un simple ajuste.** Si al añadir el orden 6 el astigmatismo pasara
de 8.96 a 5 µm, no habría nada que reportar.

---

## 7. Resultados de este ojo

Apertura de 3.600 mm de radio (7.20 mm de diámetro), 253 nodos, mapa de 59.70 µm de pico a valle.

```
Z(2,+2)    +8.956 um   astigmatismo 0/90
Z(2,-2)    +4.910 um   astigmatismo oblicuo
Z(1,-1)    -3.893 um   inclinacion Y
Z(1,+1)    -3.755 um   inclinacion X
Z(3,+3)    +0.854 um   trefoil
Z(3,+1)    +0.524 um   coma horizontal
Z(3,-3)    +0.431 um   trefoil
Z(2,+0)    -0.381 um   desenfoque
Z(4,+0)    +0.286 um   aberracion esferica
Z(3,-1)    -0.218 um   COMA VERTICAL
residuo tras Zernike-6: 0.313 um rms
```

### Interpretación

**Dominan los dos astigmatismos** (8.96 y 4.91 µm), que combinados dan 10.2 µm rms de orden 2.
El ajuste geométrico independiente lo confirma: 2.28 D de astigmatismo corneal con el meridiano
curvo a 104°.

**El desenfoque es casi nulo** (−0.38 µm). No es un hallazgo sobre el ojo: es una consecuencia
de cómo se construyó la referencia. El paso de la rejilla sintética se ajustó al propio patrón
medido, lo que fuerza a cero la diferencia media de magnificación. Por eso este experimento no
da queratometría absoluta — ver `codefiles/rejilla_referencia.py`.

**La inclinación es pequeña** (5.4 µm rms frente a los 85.4 µm que arrastraba el experimento
23-11). También por construcción: la referencia se ancló al nodo central del patrón medido.

**No hay firma de queratocono.** El orden 3 vale 1.11 µm rms y la coma vertical −0.218 µm. Para
comparar, el queratocono de 23-11 daba 6.20 µm de orden 3 y −5.82 µm de coma vertical: un factor
27 en la coma. Y la dispersión entre métodos aquí es del mismo orden que la propia señal:

```
coma vertical Z(3,-1)          astigmatismo Z(2,2)
   Southwell  -0.218 um            +8.956 um
   Campbell   -0.122 um            +9.067 um
   Malacara   +0.124 um            +9.202 um
```

La coma **cambia de signo** entre métodos, así que no es medible con estos datos. El
astigmatismo, en cambio, coincide dentro del 2.7 %: ese sí es un resultado sólido.

Este ojo tiene astigmatismo regular apreciable y aberración de alto orden en el nivel del ruido.
Es el comportamiento esperable de un ojo sano.

---

## 8. Limitaciones que deben declararse

1. **El desenfoque no es medible** en este experimento, por la construcción de la referencia
   sintética. Cualquier afirmación sobre potencia corneal absoluta sería circular.
2. **El orden 3 y superiores están al borde del ruido.** Residuo no explicado 0.313 µm frente a
   1.112 µm de orden 3: SNR ≈ 3.6, y la dispersión entre métodos de reconstrucción es comparable
   a los propios coeficientes. Reportables como cota superior, no como medida.
3. **La ortogonalidad es aproximada** (8 % de acoplamiento máximo). Irrelevante para los órdenes
   1–3, que es lo que se reporta, pero crecería si se ampliara la serie más allá de $n=6$ sobre
   esta apertura.
4. **El radio de normalización es el máximo de la apertura** ($r_{\max} = 3.600$ mm). Los
   coeficientes de Zernike **no son comparables** con los de otro trabajo si el radio de
   normalización difiere: escalan con potencias de $r_{\max}$. Siempre hay que reportarlo.

---

## Apéndice — Modos ajustados

Con `n_max = 6` se ajustan 28 modos: los $(n,m)$ con $0\le n\le 6$, $|m|\le n$ y $n-|m|$ par.

| $n$ | modos | cuántos |
|---|---|---|
| 0 | $(0,0)$ | 1 |
| 1 | $(1,\pm1)$ | 2 |
| 2 | $(2,0), (2,\pm2)$ | 3 |
| 3 | $(3,\pm1), (3,\pm3)$ | 4 |
| 4 | $(4,0), (4,\pm2), (4,\pm4)$ | 5 |
| 5 | $(5,\pm1), (5,\pm3), (5,\pm5)$ | 6 |
| 6 | $(6,0), (6,\pm2), (6,\pm4), (6,\pm6)$ | 7 |
| | | **28** |

Con 253 nodos y 28 incógnitas el sistema está sobredeterminado por un factor 9, que es holgado.
