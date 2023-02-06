# ORALE CALCOLO NUMERICO

I numeri successivi alla **X[0-9]** sono le volte in cui in diversi orali
ha fatto quelle domande, quindi è più probabile che chieda.

Lista di domande:

1. [Polinomio interpolante (**X8**)](#polinomio-interpolante-x8)
2. [Condizionamento nelle sistemi lineari (**X5**)](#condizionamento-nelle-sistemi-lineari-x5)
3. [Condizionamento del problema in approssimazione (**X5**)](#condizionamento-del-problema-in-approssimazione-x5)
4. [Spline (def, teorema dimensionalità, proprietà delle cubiche) (**X4**)](#spline-def-teorema-dimensionalità-proprietà-delle-cubiche-x4)
5. [Risoluzione sistemi lineari (**X4**)](#risoluzione-sistemi-lineari-x4)
6. [Newton-Cotes (**X4**)](#newton-cotes-x4)
7. [Chebyshev (**X4**)](#chebyshev-x4)
8. [Matrici diagonali dominanti + dim (**X4**)](#matrici-diagonali-dominanti--dim-x4)
9. [Cancellazione numerica e somma algebrica (**X3**)](#cancellazione-numerica-e-somma-algebrica-x3)
10. [Dimostrazione $A=LDL^T$ (**X2**)](#dimostrazione-x2)
11. [Precisione di macchina (**X2**)](#precisione-di-macchina-x2)
12. [Condizionamento (di un problema) cap 1 (**X2**)](#condizionamento-di-un-problema-cap-1-x2)
13. [Fattorizzazione $QR$ (**X2**)](#fattorizzazione-x2)
14. [Matrice ortogonale (**X2**)](#matrice-ortogonale-x2)
15. [Dimostrazione $sdp \Rightarrow LU$ (**X2**)](#dimostrazione-x2-1)
16. [def splitting regolare di matrice (**X2**)](#def-splitting-regolare-di-matrice-x2)
17. [metodo di newton (dim convergenza, molteplicità, Aitken) (**X2**)](#metodo-di-newton-dim-convergenza-molteplicità-aitken-x2)
18. [condizionamento del problema nei sistemi lineari](#condizionamento-del-problema-nei-sistemi-lineari)
19. [calcolare la norma uno e la norma infinito di una matrice](#calcolare-la-norma-uno-e-la-norma-infinito-di-una-matrice)
20. [matrici triangolari con codice](#matrici-triangolari-con-codice)
21. [come si ottiene $Va=f$ (risp: prodotto scalare tra gli elementi di $V$ ed $a$)](#come-si-ottiene-risp-prodotto-scalare-tra-gli-elementi-di-ed)
22. [Differenze tra $LU$ e $LDL^T$](#differenze-tra-e)
23. [equazione retta tangente](#equazione-retta-tangente)
24. [def norma indotta su matrice](#def-norma-indotta-su-matrice)
25. [codice matrice ortogonale](#codice-matrice-ortogonale)

Domande fuori corso:

1. [Metodo iterativo Google pagerank (**X4**)](#metodo-iterativo-google-pagerank-x4)
2. [Raggio spettrale (**X4**)](#raggio-spettrale-x4)
3. [Metodo iterativo applicato ai sistemi lineari (**X2**)](#metodo-iterativo-applicato-ai-sistemi-lineari-x2)
4. [metodi iterativi per sistemi lineari](#metodi-iterativi-per-sistemi-lineari)
5. [M matrici e matrici monotone in generale (**X2**)](#m-matrici-e-matrici-monotone-in-generale-x2)

[Return to questions list](#orale-calcolo-numerico)

## Polinomio interpolante (**X8**)

> Approssimazione polinomiale (inizio cap 4)

Ascisse tra di loro distinte:
$$a \leq x_0 < x_1 < ... < x_n \leq b$$

### Definizione

Un **_polinomio interpolate_** di $f(x)$ sulle ascisse, è $p(x_i) = f_i$, dove
$f_i \equiv f(x_i)$, ed $i = 0,1,...,n$.

### Teorema 4.1 - Esistenza e unicità di $p(x)$

Date le ascisse, esiste ed è unico il polinomio $p(x) \in \prod_n$ che
soddisfa definizione di polinomio interpolante.

##### Dimostrazione

Un generico $p(x) \in \prod_n$ avrà la seguente forma:
$$p(x) = \sum_{k = 0}^n a_k x^k$$
in cui coefficienti $\{a_k\}$ sono da determinare, in modo da soddisfare
la definizione. Cosi facendo, si perviene al seguente sistema di equazioni
lineari $Va = f$ in cui:

$$
\begin{pmatrix}
  x_0^0 & x_0^1 & ... & x_0^n \\
  x_1^0 & x_1^1 & ... & x_1^n \\
  \vdots & \vdots &  & \vdots \\
  x_n^0 & x_n^1 & ... & x_n^n
\end{pmatrix},\quad a=\begin{pmatrix}
  a_0 \\
  a_1 \\
  \vdots \\
a_n\end{pmatrix},\quad f=\begin{pmatrix}
  f_0 \\
  f_1 \\
  \vdots \\
  f_n
\end{pmatrix}
$$

La matrice $V$ è una matrice di $Vandermonde$ (trasposta), che è univocamente
definita dalle ascisse $\{x_i\}$. Una delle proprietà di essa è
$$det(V) = \prod_{i > j}(x_i - x_j)$$
Dato che le ascisse sono tra loro distante, allora $V$ risulta essere
nonsingolare. Pertanto, esiste ed è unica la soluzione dei sistema lineare,
ovvero esiste ed è unico il polinomio soddisfacente la definizione di
**_polinomio interpolante_**

[Return to questions list](#orale-calcolo-numerico)

## Condizionamento nelle sistemi lineari (**X5**)

> Condizionamento delle matrici

### Condizionamento del problema della valutazione del polinomio interpolante

In questa sezione, studieremo in che modo perturbazioni sui dati del sistema
lineare (3.1), con $m = n$, si ripercuotono sulla sua soluzione. Più in
dettaglio, studieremo il sistema lineare perturbato

$$(A+\Delta A)(x+\Delta x)=b+\Delta b$$

in cui le perturbazioni sui dati, $\Delta A$ e $\Delta b$, determinano la
perturbazione $\Delta x$ sulla soluzione di sistema lineare. Per semplicità di
esposizione, senza tuttavia perdere in generalità. considereremo il caso in cui
i dati del problema dipendono da un parametro scalare di perturbazione
$\varepsilon \approx 0$

$$A(\varepsilon) = A + \varepsilon F, F \in \R^{n \times n} \Longrightarrow \Delta A = \varepsilon F,$$

$$b(\varepsilon) = b + \varepsilon f, f \in \R \Longrightarrow \Delta b = \varepsilon f,$$

Conseguentemente la sistema lineare perturbato sarà

$$A(\varepsilon)x(\varepsilon)=b(\varepsilon)$$

in qui $x(\varepsilon)$ è la corrispondente soluzione del sistema lineare
perturbato. Osserviamo che

$$A(0) = A, b(0) = b \Longrightarrow x(0) = x$$

Sviluppando nell'origine, si ottiene pertanto che, per $\varepsilon$
sufficientemente piccolo

$$x(\varepsilon) = x + \varepsilon \dot{x}(0) + O(\varepsilon^2) \approx x + \varepsilon \dot{x}(0)$$

ovvero,

$$\Delta x \equiv x(\varepsilon) - x \approx \varepsilon \dot{x}(0)$$

Dalla sistema lineare perturbato segue, inoltre (derivata del prodotto)
$(Ax)' = b' \rightarrow A'x + Ax' = b'$

$$\dot{A}(\varepsilon)x(\varepsilon) + A(\varepsilon)\dot{x}(\varepsilon) = \dot{b}(\varepsilon)$$

e, quindi

$$\dot{A}(0)x + A\dot{x}(0)= \dot{b}(0)$$

Considerando il caso in cui i dati del problema dipendono da un parametro
scalare di perturbazione, questa permette di ottenere

$$\dot{x}(0) = A^{-1}(f-Fx)$$

Sostituendo in sviluppo nell'origine con formula precedente si ottiene infine che:

$$\frac{\Vert \Delta x\Vert }{\Vert x\Vert } \approx \frac{\Vert A^{-1}(\varepsilon f - \varepsilon Fx)\Vert }{\Vert x\Vert } \equiv \frac{\Vert A^{-1}(\Delta b - \Delta Ax)\Vert }{\Vert x\Vert } \leq \frac{\Vert A^{-1}\Vert (\Vert \Delta b\Vert  + \Vert \Delta A\Vert  \Vert x\Vert )}{\Vert x\Vert } =$$

$$= \Vert A^{-1}\Vert \left(\frac{\Vert \Delta b\Vert }{\Vert x\Vert } + \Vert \Delta A\Vert \right) = \Vert A\Vert \Vert A^{-1}\Vert  \left( \frac{\Vert \Delta b\Vert }{\Vert A\Vert  \cdot \Vert x\Vert } + \frac{\Vert \Delta A\Vert }{\Vert A\Vert } \right) \leq$$

$$\leq \Vert A\Vert \Vert A^{-1}\Vert  \left( \frac{\Vert \Delta b\Vert }{\Vert b\Vert } + \frac{\Vert \Delta A\Vert }{\Vert A\Vert } \right)$$

Considerato che:

- $\frac{\Vert \Delta x\Vert }{\Vert x\Vert }$ può essere, in senso lato, assimilato ad una
  sorta di “errore relativo” sul risultato e, similmente.
- $\frac{\Vert \Delta A\Vert }{\Vert A\Vert }$ e $\frac{\Vert \Delta b\Vert }{\Vert b\Vert }$ possono essere
  assimilati a corrispondenti “errori relativi’ sui dati di ingresso.

si ottiene che la quantità

$$\kappa(A) \equiv \Vert A\Vert  \cdot \Vert A^{-1}\Vert $$

definisce il numero di condizionamento del problema.

[Return to questions list](#orale-calcolo-numerico)

## Condizionamento del problema in approssimazione (**X5**)

Consideriamo che le ascisse come parametri fissati, riguardando le $f_i$ come
gli unici dati di ingresso. In questo caso l'analisi di **_condizionamento_**
verrà condotta sugli errori assoluti. Sono dunque:

$$p(x)=\sum_{k=0}^n f_k L_{k n}(x)$$
$$\tilde{p}(x)=\sum_{k=0}^n \tilde{f_k} L_{k n}(x)$$

i polinomi interpolanti, esperessi nella forma di Lagrange, costruiti a
partire dai dati esatti $f_i$ e quelli perturbati $\tilde{f_i}$.
Si ottiene pertanto:

$$
\begin{align*}
  \vert p(x) - \tilde{p}(x)\vert
    & = {forma\ completa} = {si\ racoglie\ Lagrange} \\
    & \leq {\vert Lagrange\vert \ e\ \vert (f_k - \tilde{f_k})\vert } \\
    & \leq {si\ prende\ come\ costante\ il\ max_k(f_k - \tilde{f_k})} \\
    & = \{{\lambda\ =\ \sum\ Lagrange}\} \\
    & = {\lambda_n(x) * max_k(f_k - \tilde{f_k})}\\
\end{align*}
$$

in qui $\lambda_n(x)$ è detta funzione di Lebesgue.

Pertanto definiamo la norma ($\infty$) in $C^{(0)}$

$$\Vert f\Vert  = \underset{a \leq x \leq b}{max}\vert f(x)\vert $$

allora

$$\Vert p - \tilde{p}\Vert \leq \Vert \lambda_n\Vert  * \Vert f - \tilde{f}\Vert \equiv \Lambda_n * \Vert f - \tilde{f}\Vert$$

La costante di Lebesgue $\Lambda_n$, che misura la massima amplificazione
sul risultato dell'errore sui dati di ingresso, definisce, pertanto,
il **_numero di condizionamento_** del problema,

- $\Lambda_n \geq O(\log n) \rightarrow \infty$, per $n \rightarrow \infty$. Pertanto problema diventa
  progressivamente malcondizionato, al crescere di `n`.

[Return to questions list](#orale-calcolo-numerico)

## Spline (def, teorema dimensionalità, proprietà delle cubiche) (**X4**)

$C^{(k)}$ denota l'insieme delle funzioni $f: \R \rightarrow \R$ derivabili
$k$ volte, con derivata $k$-esima continua. Ove necessario, il dominio di $f$
può essere ristretto ad un particolare intervallo $[a, b] \subset \R$

$$\Delta = \{a = x_0 < x_1 < ... < x_n = b\}$$

### Definizione

La nuova funzione sarà una funzione _polinomiale_ a _tratti_. Più interpolante
esattamente, se:

1. $s_m(x) \in C^{m-1}$ sull'intervallo $[a,b]$, e inoltre,
2. $s_m\vert_{[x_{i-1}, x_i]}(x) \in \prod_m,\quad i = 1, ..., n$

Allora diremo che $s_m(x)$ è una **_spline di grado m_** sulla partizione
$\Delta$. Se, inoltre

$$s_m(x_i) = f_i,\quad i = 0, 1, ..., n$$

Allora diremo che la **_spline_** interpola la funzione $f(x)$ nei nodi di
tale partizione.

### Teorema (4.10)

Se $s_m(x)$ è una spline di grado $m$ sulla partizione $\Delta$, allora
$s_m'(x)$ è una spline di grado $m - 1$ sulla stessa partizione.

### Teorema (4.11?) dimensionalità

L'insieme delle funzioni `spline` di grado $m$ definite sulla partizione
($\Delta$) è uno spazio vettoriale di dimensione $m + n$.

Una conseguenza di questa teorema è che sono necessarie $m + n$ condizioni
(indipendenti) per individuare univocamente la spline interpolante una
funzione sulla partizione $\Delta$ assegnata.

### Proprieta delle spline cubiche

1. Spline **naturale**
   $$s_3''(a) = 0,\quad s_3''(b) = 0$$

2. Spline **completa**
   $$s_3'(a) = f'(a),\quad s_3'(b) = f'(b)$$

3. Spline **periodica**
   $$s_3'(a) = s_3'(b),\quad s_3''(a) = s_3''(b)$$

- Condizioni **not-a-knot**
  In questo caso, per evitare altre 3 condizioni precedenti, basta che vale
  una di questi due (sono equivalenti, osservando che, in virtù del teorema
  4.10, $s_3'''\vert_{x_{t-1},x_t}(x) \in \prod_0$):

  $$s_3'''\vert_{[x_0,x_1]}(x_1) = s_3'''\vert_{[x_1,x_2]}(x_1),\quad s_3'''\vert_{[x_{n-2},x_{n-1}]}(x_{n-1}) = s_3'''\vert_{[x_{n-1},x_{n}]}(x_{n-1})$$

  $$\frac{s_3''(x_1) - s_3''(x_0)}{x_1 - x_0} = \frac{s_3''(x_2) - s_3''(x_1)}{x_2 - x_1},\quad \frac{s_3''(x_{n-1}) - s_3''(x_{n-2})}{x_{n-1} - x_{n-2}} = \frac{s_3''(x_n) - s_3''(x_{n-1})}{x_n - x_{n-1}}$$

[Return to questions list](#orale-calcolo-numerico)

## Risoluzione sistemi lineari (**X4**)

> I sistemi lineari in generale, cap 3 fino a dim unicità di $A=LU$ (**X2**)

### Sistemi lineari

E noto che tali sistemi possono essere scritti nella forma $A x = b$.
Dove $A = (a_{ij}) \in \R^{m \times n}$ è la matrice dei coefficienti,
$b = (b_i) \in \R^m$ è il vettore dei termini noti e, infine,
$x = (x_i) \in \R^n$ è il vettore delle incognite.
Pertatno, la soluzione del sistema lineare esiste ed è unica
$$x = A^{-1}b$$
Tuttavia, questa espressione formale della soluzione non induce,
generalmente, un metodo di risoluzione efficiente.

Nella trattazione seguente sarà sempre assunto che $m > n$ e, inoltre,
che la matrice $A$ abbia rango massimo, ovvero $rank(A) = n$. Queste
assunzioni coprono una significativa parte dei problemi che derivano
dalle applicazioni.

### Casi semplici

1. Matrici diagonali
   - $x_i = \frac{b_i}{a_{ii}}$
   - $n$ flop
   - $n$ spazio
2. Matrici triangolari
   - $x_i = \frac{b_i - \sum_{j=1}^{i-1}(a_{ij} x_j)}{a_{nn}}$
   - $\sim n^2$ flop
   - $\frac{n^2}{2}$ spazio
3. Matrici ortogonali
   - $x = A^T b$, in matrici ortogonali $A^{-1} = A^T$
   - $2n^2$ flop
   - $n^2$ spazio

### Fattorizazione

Idea base ottenere una fattorizzazione della matrice di coeficienti
$A = F_1 * F_2 * ... * F_k$ con $F_i \in R^{n \times n}$
dove i fattori ( $F_i$ ) sono matrici nonsingolari. Soluzione può
essere calcolato, risolvendo seguenti sistemi lineari:

$$F_1 * x_1 = b$$
$$F_2 * x_2 = x_1$$
$$...$$
$$F_k * x_k = x_{k-1}$$
$$x \equiv x_k$$

Acluni metodi di fattorizzazione:

- $LU$
- $QR$
- $LDL^T$ per simmettriche definite positive

[Return to questions list](#orale-calcolo-numerico)

## Newton-Cotes (**X4**)

Si considera l'approssimazione di $f(x)$ fornita dal polinomio interpolante
su $n+1$ ascisse equidistanti.

$$p(x_i) = f(x_i),\quad i = 0, 1, ..., n$$
$$x_i = a + i * h,\quad h = \frac{b - a}{n}$$

Considerando la forma di Lagrange di polinomio, si ha:

$$I(f) \approx \int_a^b \sum_{k=0}^{n}(f_k L_{kn}(x))dx =\sum_{k=0}^n(f_k \int_a^b L_{kn}(x)dx\ = h \sum_{k=0}^n(f_k \int_a^b \prod_{j=0,j \neq k}^{n}\frac{t - j}{k - j}dt)$$

Nel ultimo passaggio utilizzata la trasformazione $x_t = a + th$.
Pertanto la formula è
$$I_n(f) \equiv \frac{b - a}{n} \sum_{k=0}^n(c_{kn} * f_k)$$
in cui
$$c_{kn} = \int_a^b\prod_{j=0,j \neq k}^{n}\frac{t - j}{k - j}dt,\quad k = 0, 1, ..., n$$
definisce l'approssimazione di $I(f)$ cercata. Essa difinisce la generica
_**formula di quadratura** di Newton-Cotes_.

Si distinguano casi:

- $n = 1$ allora è la _formula dei trapezi_
- $n = 2$ allora è la _formula di Simpson_

[Return to questions list](#orale-calcolo-numerico)

## Chebyshev (**X4**)

### Definizione

$$T_0(x) \equiv 1$$
$$T_1(x) = x$$
$$T_{k+1}(x) = 2xT_k(x) - T_{k - 1}(x),\quad k = 1, 2, ...$$

1. $T_k(x)$ è un polinomio di grado esatto $k$,
2. Il coefficiente principale di $T_k(x)$ è $2^{k - 1},\quad k = 1, 2, ...$
3. La famiglia di polinomi { $\hat{T}_k$ }, in cui
   $$\hat{T}_0(x) = T_0(x),\quad \hat{T}_k(x) = 2^{1 - k}T_k(x),\quad k = 1, 2, ...$$
4. Ponendo $x = \cos \theta,\quad \theta \in [0,\pi].$ per parametrizzare
   i punti dell'intervallo [-1,1] rispetto a $\theta$, e considerando che
   $\cos(k\theta+\theta) + \cos(k\theta-\theta) = 2\cos(k\theta)\cos(\theta)$,
   si ottiente:
   $$T_k(x) \equiv T_k(\cos(\theta)) = \cos(k\theta),\quad k = 0, 1, ...$$

### Gli zeri

Gli zeri di $T_k(x)$, tra loro tutti distinti, sono dati da:
$$x_{i}^{(k)} = \cos\left(\frac{(2i + 1)\pi}{2k}\right), \quad i = 0, 1, ..., k - 1$$

Inoltre, la costante di Legesgue è $\Lambda_n \approx \frac{2}{\pi}\log n$

[Return to questions list](#orale-calcolo-numerico)

## Matrici diagonali dominanti + dim (**X4**)

### Definizione

Data una matrice $A = (a\_{ij}) \in \R^{x \times n}$, si dice che essa è:

- diagonale dominate per righe se
  $$\vert a_{ii}\vert  > \sum_{j \neq i}\vert a_{ij}\vert ,\quad i = 1, ..., n$$
- diagonale dominate per colonne se
  $$\vert a_{ii}\vert  > \sum_{j \neq i}\vert a_{ji}\vert ,\quad i = 1, ..., n$$

### Lemma 3.4

Se una matrice $A$ è diagonale dominante per righe (rispettiva-mente, per colonne), allora tali sono tutte le sui sotto matrici principali.

##### Dimostrazione lemma 3.4

Dalla definizione di sottomatrici principali, si ha che tali sono quelli che
alla diagonale principale hanno elementi che sono anche elementi di diagonale
principale di martrice sorgente.

Dalla definizione di sottomatrici principale è ovvio che righe e collone di
matrice sorgente non cambiano tranne perdere elementi, ma questo significa
che somma di valori assoluti diminuisce o rimane uguale (in caso elementi
nulli), ma valore di elementi in diagonale rimane sempre stesso.

### Lemma 3.5

Una matrice $A$ è diagonale dominate per riche (rispettivamente, per colonne)
se e solo se $A^T$ è diagonale dominante per colonne(rispettivamente, per righe)

##### Dimostrazione lemma 3.5

Questo è ovvio, dato che matrice trasposta nient'altro che la matrice sorgente
con collone fatti di righe, e righe fatte di collone, questo significa che
solamente elementi in prindcipale diagonale rimangono uguali, ma da qui
è ovvio che se matrice principale è diagonale dominante per righe (colonne),
allora sua trasposta sarà diagonale dominate per colonne (righe).

[Return to questions list](#orale-calcolo-numerico)

## Metodo iterativo Google pagerank (**X4**)

Prolbema può essere riformulata in seguente sistema lineare:
$$A\hat{x} \equiv (I - pS)\hat{x} = \frac{1-p}{n}e \equiv b$$
Dalla dimensione di $A$ è impensabile applicare la fattorizzazione diretta.
La matrice $A$ ha una importate caratteristica, essere scritta in forma:
$$A = I - B,\quad B \geq 0,\quad \rho(B) < 1$$
Infatti, nel nostro caso, $S \geq 0$, $\rho(S) = 1$ e $p < 1$.

### Splitting regolari di matrici

$$M^{-1} \geq 0,\quad N \geq 0$$

### Lemma 6.1

Siano $A, B \in \R^{n \times n}$, $A \geq B \geq 0$. Allora
$A^i \geq B^i \geq 0$, $i \geq 0$.

### Lemma 6.2

Siano $A, B \in \R^{n \times n}$, $A \geq B \geq 0$. Allora
$\rho(A) \geq \rho(B)$.

### Criterio di arresto

$$r_k = Ax_k - b$$

### I metodi di Jacobi e Gauss-Seidel

$$A = D - L - U$$
in cui:

- $D$ è diagonale;
- $L$ è strettamente triangolare inferiore;
- $U$ è strettamente triangolare superiore.

[Return to questions list](#orale-calcolo-numerico)

## Raggio spettrale (**X4**)

$$ \hat{x} = (H + v \Delta^T)\hat{x} \equiv S\hat{x},$$
  dove
  $$v = \frac{1}{n}e,\quad e = (1, ..., 1)^T \in \R^n.$$

### Teorema 6.1

La matrice $S$ definita prima soddisfa le seguenti proprietà (e
disuguaglianze si intendono valere per ogni elemento):

1. $S \leq 0$;
2. $e^T S = e^T$;
3. $\lambda = 1$ è il **_raggio spettrale_** di $S$.

### Dimostrazione

La dimostrazione dei primi due punti è immediata. Riguardo al'ultimo punto,
osserviamo che dal secondo segue che $\lambda = 1$ è autovalore di $S^T$ e,
quindi, di $S$. Osservando che $p(S) \leq \\vert S\\vert$ per ogni norma indotta su
matrice, la tesi si completa in virtù del punto $1$, da cui si ottiene:
$1 = \\vert e^T S\\vert_{\infty} = \\vert S\\vert_1$.

[Return to questions list](#orale-calcolo-numerico)

## Cancellazione numerica e somma algebrica (**X3**)

**_Cancellazione numerica_** è la conseguenza più grave della rappresentazione
con precisione finita dei numeri reali all'interno di un calcolatore.

Tale fenomeno consiste nella perdita di cifre significative, dovuta a un'operazione di sottrazione tra due numeri "quasi uguali". Il termine "quasi uguali", indica che i due operandi hanno le prime $t$ cifre uguali con $t \in \N$, $t > 0$.

Con **_somma algebrica_** si intende l'operazione di addizione o sottrazione di
numeri complessi (quindi anche reali e a maggior ragione anche interi).

Il problema è quello di stuidare il condizionamento di
$$y = x_1 + x_2,\quad x_1, x_2 \in \R,\quad x_1 + x_2 \neq 0$$
Denotando con $\varepsilon_1$ e $\varepsilon_2$ gli errori relativi sui dati
iniziali, ed assumendo che nessun nuovo errore venga introdotto nel calcoli,
si ottiene:

$$y(1 + \varepsilon_y) = x_1(1 + \varepsilon_1) + x_2(1 + \varepsilon_2) = x_1 + x_2 + x_1\varepsilon_1 + x_2\varepsilon_2$$

Si ricava:

$$\vert \varepsilon_y\vert  \leq \frac{\vert x_1\vert  + \vert x_2\vert }{\vert x_1 + x_2\vert }\varepsilon_x \equiv \kappa \varepsilon_x,\quad \varepsilon_x = max\{\vert \varepsilon_1\vert , \vert \varepsilon_2\vert \}$$

[Return to questions list](#orale-calcolo-numerico)

## Dimostrazione $A=LDL^T$ (**X2**)

### Matrice simmetrica definita positiva

Una matrice $A \in \R^{n \times n}$ è sdp se è simmetrica (cioè, $A = A^T$) e,
per ogni $x \in \R^n$, $x \neq 0$, risulta $x^T A x > 0$

#### Lemma 3.7

Tutte le sottomatrici principali di una matrice sdp sono sdp.

#### Lemma 3.8

Una matrice sdp è nonsingolare.

### Teorema 3.4

Se $A$ è sdp, allora è fattorizzabile $LU$.

#### Dimostrazione sdp $LU$ fattorizzabile

Dal Lemma 3.7, tutte le sottomatrici principali sono
sdp e quindi, dal Lemma 3.8, segue che i corrispondenti minori
principali sono tutti non nulli.

### Teorema 3.5

Gli elementi diagonali di una matrice sdp sono positivi.

### Teorema 3.6

A è sdp se e solo se
$$A = LDL^T$$
Dove

- $L$ triangolare inferiore a diagonale unitaria.
- $D$ diagonale con elementi diagonali positivi.

#### Dimostrazione

$A$ è fattorizzabile $LU$. Inoltre il fattore $U$ può essere scritto
nella forma
$$U = D\hat{U}$$
con $D$ diagonale e $\hat{U}$ triangolare superiore a diagonale unitaria.
Essendo, inoltre $A = A^T$, segue pertanto che:
$$LD\hat{U} = A = A^T = (LD\hat{U})^T = \hat{U}^T D L^T$$
Per l'unicità della fattorizzazione $LU$, essendo $\hat{U}^T$ triangolare
inferiore a diagonale unitaria e $D L^T$ triangolare superiore, segue quindi
che $\hat{U}^T = L$. Pertanto la fattorizzazione di teorema è ben definita.
Rimane da dimostrare che gli elementi diagonali di $D$ sono positivi. In
virtù del Teorema 3.5, basta dimostrare che $D$ è sdp. Evidentemente, $D$
è simmetrica. Inoltre, comunque si fissi $x \neq 0$, esiste ed è unico il
vettore $y \neq 0$ tale che $L^T y = x$. Segue pertanto che
$$x^T D x = (L^T y)^T D (L^T y) = y^T LDL^T y = y^T A y > 0$$
essendo $A$ sdp

[Return to questions list](#orale-calcolo-numerico)

## Precisione di macchina (**X2**)

### Teorema 1.3

Il più piccolo ed il più grand (in valore assoluto), tra i numeri di
macchina diversi da 0, sono rispettivamente dati da:

$$
\begin{align*}
  &r_1 = b^{-v},\\
  &r_2 = (1 - b^{-m})b^{\varphi},\quad \varphi = b^s - v
\end{align*}
$$

Numeri di machina sono $L = [-r_2, -r_1]\cup\{0\}\cup[r_1, r_2]$

### Teorema 1.4

Se $x \in L$, $x \neq 0$, allora
$$fl(x) = x(1 + \varepsilon_x),\quad \vert \varepsilon_x\vert  \leq u$$
dove

$$u = \begin{cases} b^{1-m},\ in\ caso\ di\ troncamento,\\ \frac{1}{2} b^{1-m},\ in\ caso\ di\ arrotondamento \end{cases}$$

##### Dimostrazione

La **_precisione di macchina_** è definita da quantità $u$ in Teorema 1.4

[Return to questions list](#orale-calcolo-numerico)

## Condizionamento (di un problema) cap 1 (**X2**)

> Condizionamento del problema in generale

$$\vert \varepsilon_y\vert  \approx \left\vert  f'(x)\frac{x}{y}\right\vert  \vert \varepsilon_x\vert \ \equiv \kappa \vert \varepsilon_x\vert$$

Il fattore di amplificazione $\kappa$, che misura di quanto gli errori
iniziali possono amplificarsi sul risultato finale, è denominato
**_numero di condizione_** del problema.

In generale, si distinguono i seguenti casi significativi:

- $\kappa \approx 1$: gli errori sul risultato finale cono dello stesso
  ordine di quelli iniziali. In tal caso il problema si dice
  **ben condizionato**
- $\kappa \gg 1$: gli errori sul risultato finale possono essere assai più
  grandi degli errori iniziali. In questo caso, il problema si dice
  **malcondizionato**

Osservare che:

  <!-- 1. nel caso in cui si utilizzi una precisione di macchina $u$ e si abbia
     $\kappa \approx u^{-1}$, qualunque risultato sarà privo di significato -->

[Return to questions list](#orale-calcolo-numerico)

## Fattorizzazione $QR$ (**X2**)

$$Ax = b,\quad A \in \R^{m \times n},\quad m > n \equiv rank(A)$$

### Teorema 3.8

Data la matrice $A$, esistono:

- $Q \in \R^{m \times m}$, ortogonale
- $\hat{R} \in \R^{n \times n}$, triangolare superiore e non singolare

tali che:
$$A = QR \equiv Q\binom{\hat{R}}{O}$$

[Return to questions list](#orale-calcolo-numerico)

## Metodo iterativo applicato ai sistemi lineari (**X2**)

[Google pagerank](#metodo-iterativo-google-pagerank-x4)

[Raggio spettrale](#raggio-spettrale-x4)

[Return to questions list](#orale-calcolo-numerico)

## Matrice ortogonale (**X2**)

**_Matrice ortogonale_** è una matrice invertibile tale che la sua trasposta coincide con la sua inversa.

[Return to questions list](#orale-calcolo-numerico)

## M matrici e matrici monotone in generale (**X2**)

Le $M$-matrici sono particolari matrici $monotone$, in quanto, se $A$ è una
$M$-matrice allora
$$Ax \leq C\quad \Rightarrow\quad I \leq A^{-1}C.\ I \leq CA^{-1}$$
dove, al solito, le diseguaglianze si intendono elemento per elemento.

[Google pagerank](#metodo-iterativo-google-pagerank-x4)

[Raggio spettrale](#raggio-spettrale-x4)

[Return to questions list](#orale-calcolo-numerico)

## Dimostrazione $sdp \Rightarrow LU$ (**X2**)

### [Dimostrazione Teorema 3.4](#dimostrazione-sdp-fattorizzabile)

[Return to questions list](#orale-calcolo-numerico)

## def splitting regolare di matrice (**X2**)

[Return to questions list](#orale-calcolo-numerico)

## metodo di newton (dim convergenza, molteplicità, Aitken) (**X2**)

[Return to questions list](#orale-calcolo-numerico)

## condizionamento del problema nei sistemi lineari

> Condizionamento cap 2

[Return to questions list](#orale-calcolo-numerico)

## calcolare la norma uno e la norma infinito di una matrice

[Return to questions list](#orale-calcolo-numerico)

## matrici triangolari con codice

[Return to questions list](#orale-calcolo-numerico)

## come si ottiene $Va=f$ (risp: prodotto scalare tra gli elementi di $V$ ed $a$)

[Return to questions list](#orale-calcolo-numerico)

## Differenze tra $LU$ e $LDL^T$

### LU

Se la matrice A può essere scritta come il prodotto di due fattori
$$A = LU$$
con $L$ triangolare inferiore a diagonale unitaria, e $U$ triangolare
superiore, allora si dice che fattorizzabile LU.

#### Teorema 3.1 Unicità della fattorizzazione $LU$

Se la fattorizzazione esiste e A è nonsingolare, allora essa è unica

##### Dimostrazione

Infatti, se $A = L_1 U_1 = L_2 U_2$ fossero due fattorizzazioni $LU$ di $A$
allora seguirebbe che
$$0 \neq det(A) = det(L_2 U_2) = det(L_2)det(U_2) = det(U_2)$$
Pertanto $U_2$ è nonsingolare e, quindi,
$$L_1^{-1}L_2 = U_1 U_2^{-1} \equiv D$$
Tuttavia, essendo $L_1^{-1}L_2$ triangolare inferiore e $U_1 U_2^{-1}$
triangolare superiore, segue che $D$ è diagonale. Inoltre, essendo la
diagonale di $L_1^{-1}L_2$ unitaria, tale è anche quella di $D$, ovvero
$D = I$. Discende quindi immediatamente che $L_1 = L_2$ e $U_1 = U_2$
$.\square$

#### Teorema 3.1 Unicità della fattorizzazione $LU$

Se $A$ è nonsingolare, fattorizzazione esiste se e solo se tutti i minori
principali di $A$ sono non nulli.

Per definire le $L$ e $U$ dobbiamo applicare metodo di eliminazione
di Gauss ([1](https://www.youmath.it/lezioni/algebra-lineare/matrici-e-vettori/831-eliminazione-di-gauss.html),
[2](https://it.wikipedia.org/wiki/Metodo_di_eliminazione_di_Gauss)), con

$$g_i \equiv \frac{1}{a_{ii}^{(i)}}\ (0, ..., 0, a_{i+1,i}^{(i)}, ..., a_{ni}^{(i)})^T$$

$$
L = \begin{pmatrix}
  1      & 0      & \cdots    & 0      \\
  g_{21} & 1      & \ddots    & \vdots \\
  \vdots & \ddots & \ddots    & 0      \\
  g_{n1} & \cdots & g_{n,n-1} & 1
\end{pmatrix},\quad \ U = \begin{pmatrix}
  a_{11}^{(1)} & \cdots &                     & a_{1n}^{(1)}       \\
  0            & \ddots &                     & \vdots             \\
  \vdots       & \ddots & a_{n-1,n-1}^{(n)-1} & a_{n-1, n}^{(n-1)} \\
  0            & \cdots & 0                   & a_{nn}^{(n)}
\end{pmatrix}
$$

### $LDL^T$

Per ottenere $D$ basta prendere elementi diagonali di $U$ in $LU$

$$
L = \begin{pmatrix}
  1      & 0      & \cdots    & 0      \\
  g_{21} & 1      & \ddots    & \vdots \\
  \vdots & \ddots & \ddots    & 0      \\
  g_{n1} & \cdots & g_{n,n-1} & 1
\end{pmatrix},\quad \ D = \begin{pmatrix}
  a_{11}^{(1)} & 0      & \cdots              & 0            \\
  0            & \ddots & \ddots              & \vdots       \\
  \vdots       & \ddots & a_{n-1,n-1}^{(n)-1} & 0            \\
  0            & \cdots & 0                   & a_{nn}^{(n)}
\end{pmatrix}
$$

[Come funziona](https://yewtu.be/watch?v=8JdJoc3HMA8).

[Return to questions list](#orale-calcolo-numerico)

## metodi iterativi per sistemi lineari

[Return to questions list](#orale-calcolo-numerico)

## Equazione retta tangente

> [Testo originale](https://www.youmath.it/lezioni/analisi-matematica/derivate/225-calcolare-la-retta-tangente-ad-una-funzione-in-un-punto.html)

L'obbiettivo è determinare l'equazione della retta tangente al grafico di una
funzione in un punto.

Supponiamo di avere una funzione y=f(x) e che ci venga chiesto di calcolare
l'equazione della retta tangente al suo grafico nel punto di ascissa x=x.

La procedura che descriviamo qui è semplicemente la messa in pratica di quanto
spiegato nell'articolo sul significato geometrico della derivata di una
funzione in un punto.

1. Calcoliamo la derivata della funzione $y=f(x)$, _come funzione_ , e
   chiamiamola $y=f'(x)$.
   $$calcolo\ y=f'(x)$$
2. $$valuto\ y_0=f(x_0)$$
3. $$considero\ y=mx+q$$
4. $$valuto\ m=f'(x_0)$$
5. $$f(x_0)=f'(x_0) \cdot x_0+q$$
   $$q=f'(x_0)-f'(x_0) \cdot x_0$$
6. $$y = f'(x_0) \cdot x + f(x_0) - f'(x_0) \cdot x_0$$

[Return to questions list](#orale-calcolo-numerico)

## def norma indotta su matrice

[Return to questions list](#orale-calcolo-numerico)

## codice matrice ortogonale
