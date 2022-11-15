<!-- Add LaTeX support to Markdown
<script type="text/javascript"
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/latest/MathJax.js?config=TeX-AMS_CHTML">
</script>
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({
    tex2jax: {
      inlineMath: [['$','$'], ['\\(','\\)']],
      processEscapes: true},
      jax: ["input/TeX","input/MathML","input/AsciiMath","output/CommonHTML"],
      extensions: ["tex2jax.js","mml2jax.js","asciimath2jax.js","MathMenu.js","MathZoom.js","AssistiveMML.js", "[Contrib]/a11y/accessibility-menu.js"],
      TeX: {
      extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"],
      equationNumbers: {
      autoNumber: "AMS"
      }
    }
  });
</script>
 -->

## ORALE CALCOLO NUMERICO

I numeri successivi alla **X[0-9]** sono le volte in cui in diversi orali
ha fatto quelle domande, quindi è più probabile che chieda.

- <details><summary>Polinomio interpolante (<b>X8</b>)</summary>

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
    x_0^0  & x_0^1  & ... & x_0^n  \\
    x_1^0  & x_1^1  & ... & x_1^n  \\
    \vdots & \vdots &     & \vdots \\
    x_n^0  & x_n^1  & ... & x_n^n
  \end{pmatrix},\ \ \ \ \
  a =
  \begin{pmatrix}
    a_0 \\
    a_1 \\
    \vdots \\
    a_n
  \end{pmatrix},\ \ \ \ \
  f =
  \begin{pmatrix}
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

  </details>

- <details><summary>Condizionamento nell'interpolazione (<b>X5</b>)</summary>

  ### Condizionamento del problema della valutazione del polinomio interpolante

  Consideriamo che le ascisse come parametri fissati, riguardando le $f_i$ come
  gli unici dati di ingresso. In questo caso l'analisi di **_condizionamento_**
  verrà condotta sugli errori assoluti. Sono dunque:

  $$p(x)=\sum_{k=0}^n f_k L_{k n}(x)$$
  $$p'(x)=\sum_{k=0}^n f_k' L_{k n}(x)$$

  i polinomi interpolanti, esperessi nella forma di Lagrange, costruiti a
  partire dai dati esatti $f_i$ e quelli perturbati $f_i'$.
  Si ottiene pertanto:

  $$
  \begin{align*}
  \vert p(x) - p'(x)\vert  & = {forma\ completa} = {si\ racoglie\ Lagrange} \\
  & \leq {\vert Lagrange\vert \ e\ \vert (f_k - f_k')\vert } \\
  & \leq {si\ prende\ come\ costante\ il\ max_k(f_k - f_k')} \\
  & = \{{\lambda\ =\ \sum\ Lagrange}\} \\
  & = {\lambda_n(x) * max_k(f_k - f_k')}\\
  \end{align*}
  $$

  in qui $\lambda_n(x)$ è detta funzione di Lebesgue.

  $$\Vert p - p'\Vert  \leq \Vert \lambda_n\Vert  * \Vert f - f'\Vert  \equiv \Lambda_n * \Vert f - f'\Vert $$

  La costante di Lebesgue $\Lambda_n$, che misura la massima amplificazione
  sul risultato dell'errore sui dati di ingresso, definisce, pertanto,
  il **_numero di condizionamento_** del problema,

  - $\Lambda_n \geq O(\log n) \rightarrow \infty$, per $n \rightarrow \infty$. Pertanto problema diventa
  progressivamente malcondizionato, al crescere di `n`.
  </details>

- <details><summary>Spline (def, teorema dimensionalità, proprietà delle cubiche) (<b>X4</b>)</summary>

  $C^{(k)}$ denota l'insieme delle funzioni $f: \R \rightarrow \R$ derivabili
  $k$ volte, con derivata $k$-esima continua. Ove necessario, il dominio di $f$
  può essere ristretto ad un particolare intervallo $[a, b] \subset \R$

  $$\Delta = \{a = x_0 < x_1 < ... < x_n = b\}$$

  ### Definizione

  La nuova funzione sarà una funzione _polinomiale_ a _tratti_. Più interpolante
  esattamente, se:

  1. $s_m(x) \in C^{m-1}$ sull'intervallo $[a,b]$, e inoltre,
  2. $s_m\vert_{[x_{i-1}, x_i]}(x) \in \prod_m,\ \ \ \ i = 1, ..., n$

  Allora diremo che $s_m(x)$ è una **_spline di grado m_** sulla partizione
  $\Delta$. Se, inoltre

  $$s_m(x_i) = f_i,\ \ \ \ \ i = 0, 1, ..., n$$

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
     $$s_3''(a) = 0,\ \ \ \ s_3''(b) = 0$$

  2. Spline **completa**
     $$s_3'(a) = f'(a),\ \ \ \ s_3'(b) = f'(b)$$

  3. Spline **periodica**
     $$s_3'(a) = s_3'(b),\ \ \ \ s_3''(a) = s_3''(b)$$

  - Condizioni **not-a-knot**
    In questo caso, per evitare altre 3 condizioni precedenti, basta che vale
    una di questi due (sono equivalenti, osservando che, in virtù del teorema
    4.10, $s_3'''\vert_{x_{t-1},x_t}(x) \in \prod_0$):
    $$
      s_3'''\vert_{[x_0,x_1]}(x_1) = \
      s_3'''\vert_{[x_1,x_2]}(x_1),\ \ \ \ \
      s_3'''\vert_{[x_{n-2},x_{n-1}]}(x_{n-1}) = \
      s_3'''\vert_{[x_{n-1},x_{n}]}(x_{n-1})
    $$
    $$
      \frac{s_3''(x_1) - s_3''(x_0)}{x_1 - x_0} = \
      \frac{s_3''(x_2) - s_3''(x_1)}{x_2 - x_1},\ \ \ \ \
      \frac{s_3''(x_{n-1}) - s_3''(x_{n-2})}{x_{n-1} - x_{n-2}} = \
      \frac{s_3''(x_n) - s_3''(x_{n-1})}{x_n - x_{n-1}}
    $$

  </details>

- <details><summary>Risoluzione sistemi lineari (<b>X4</b>)</summary>

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
  dove i fattori ($F_i$) sono matrici nonsingolari. Soluzione può
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

   </details>

- <details><summary>Newton-Cotes (<b>X4</b>)</summary>

  Si considera l'approssimazione di $f(x)$ fornita dal polinomio interpolante
  su $n+1$ ascisse equidistanti.

  $$p(x_i) = f(x_i),\ \ \ \ i = 0, 1, ..., n$$
  $$x_i = a + i * h,\ \ \ \ h = \frac{b - a}{n}$$

  Considerando la forma di Lagrange di polinomio, si ha:

  $$
  I(f) \approx \int_a^b \sum_{k=0}^{n}(f_k L_{kn}(x))dx\
  =\sum_{k=0}^n(f_k \int_a^b L_{kn}(x)dx\\
  = h \sum_{k=0}^n(f_k \int_a^b \prod_{j=0,j!=k}^{n}\frac{t - j}{k - j}dt)\\
  $$

  Nel ultimo passaggio utilizzata la trasformazione $x_t = a + th$.
  Pertanto la formula è
  $$I_n(f) \equiv \frac{b - a}{n} \sum_{k=0}^n(c_{kn} * f_k)$$
  in cui
  $$c_{kn} = \int_a^b\prod_{j=0,j!=k}^{n}\frac{t - j}{k - j}dt,\ \ \ \ k = 0, 1, ..., n$$
  definisce l'approssimazione di $I(f)$ cercata. Essa difinisce la generica
  _**formula di quadratura** di Newton-Cotes_.

  Si distinguano casi:

  - $n = 1$ allora è la _formula dei trapezi_
  - $n = 2$ allora è la _formula di Simpson_
  </details>

- <details><summary>Chebyshev (<b>X4</b>)</summary>

  ### Definizione

  $$T_0(x) \equiv 1$$
  $$T_1(x) = x$$
  $$T_{k+1}(x) = 2xT_k(x) - T_{k - 1}(x),\ \ \ \ k = 1, 2, ...$$

  1. $T_k(x)$ è un polinomio di grado esatto $k$,
  2. Il coefficiente principale di $T_k(x)$ è $2^{k - 1},\ \ \ \ k = 1, 2, ...$
  3. La famiglia di polinomi {$\hat{T}_k$}, in cui
     $$\hat{T}_0(x) = T_0(x),\ \ \ \ \hat{T}_k(x) = 2^{1 - k}T_k(x),\ \ \ \ k = 1, 2, ...$$
  4. Ponendo $x = \cos \theta,\ \ \ \ \theta \in [0,\pi].$ per parametrizzare
     i punti dell'intervallo [-1,1] rispetto a $\theta$, e considerando che
     $\cos(k\theta+\theta) + \cos(k\theta-\theta) = 2\cos(k\theta)\cos(\theta)$,
     si ottiente:
     $$T_k(x) \equiv T_k(\cos(\theta)) = \cos(k\theta),\ \ \ \ k = 0, 1, ...$$

  ### Gli zeri

  Gli zeri di $T_k(x)$, tra loro tutti distinti, sono dati da:
  $$x_{i}^{(k)} = \cos\left(\frac{(2i + 1)\pi}{2k}\right), \ \ \ \ i = 0, 1, ..., k - 1$$

  Inoltre, la costante di Legesgue è $\Lambda_n \approx \frac{2}{\pi}\log n$

  </details>

- <details><summary>Matrici diagonali dominanti + dim (<b>X4</b>)</summary>

  ### Definizione

  Data una matrice $A = (a\_{ij}) \in \R^{x \times n}, si dice che essa è:

  - diagonale dominate per righe se
    $$\vert a_{ii}\vert  > \sum_{j \neq i}\vert a_{ij}\vert ,\ \ \ \ i = 1, ..., n$$
  - diagonale dominate per colonne se
    $$\vert a_{ii}\vert  > \sum_{j \neq i}\vert a_{ji}\vert ,\ \ \ \ i = 1, ..., n$$

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
  </details>

- <details><summary id="google-pagerank">Metodo iterativo Google pagerank (<b>X4</b>)</summary>

  Prolbema può essere riformulata in seguente sistema lineare:
  $$A\hat{x} \equiv (I - pS)\hat{x} = \frac{1-p}{n}e \equiv b$$
  Dalla dimensione di $A$ è impensabile applicare la fattorizzazione diretta.
  La matrice $A$ ha una importate caratteristica, essere scritta in forma:
  $$A = I - B,\ \ \ \ B \geq 0,\ \ \ \ \rho(B) < 1$$
  Infatti, nel nostro caso, $S \geq 0$, $\rho(S) = 1$ e $p < 1$.

  ### Splitting regolari di matrici

  $$M^{-1} \geq 0,\ \ \ \ N \geq 0$$

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

  </details>

- <details><summary id="raggio-spettrale">Raggio spettrale (<b>X4</b>)</summary>

  $$ \hat{x} = (H + v \Delta^T)\hat{x} \equiv S\hat{x},$$
  dove
  $$v = \frac{1}{n}e,\ \ \ \ e = (1, ..., 1)^T \in \R^n.$$

  ### Teorema 6.1

  La matrice $S$ definita prima soddisfa le seguenti proprietà (e
  disuguaglianze si intendono valere per ogni elemento):

  1. $S \leq 0$;
  2. $e^T S = e^T$;
  3. $\lambda = 1$ è il **_raggio spettrale_** di $S$.

  ### Dimostrazione

  La dimostrazione dei primi due punti è immediata. Riguardo al'ultimo punto,
  osserviamo che dal secondo segue che $\lambda = 1$ è autovalore di $S^T$ e,
  quindi, di $S$. Osservando che $p(S) \leq \\vert S\\vert $ per ogni norma indotta su
  matrice, la tesi si completa in virtù del punto $1$, da cui si ottiene:
  $1 = \\vert e^T S\\vert_{\infty} = \\vert S\\vert_1$.

  </details>

- <details><summary>Cancellazione numerica e somma algebrica (<b>X3</b>)</summary>

  **_Cancellazione numerica_** è la conseguenza più grave della rappresentazione
  con precisione finita dei numeri reali all'interno di un calcolatore.

  Tale fenomeno consiste nella perdita di cifre significative, dovuta a un'operazione di sottrazione tra due numeri "quasi uguali". Il termine "quasi uguali", indica che i due operandi hanno le prime $t$ cifre uguali con $t \in \N$, $t > 0$.

  Con **_somma algebrica_** si intende l'operazione di addizione o sottrazione di
  numeri complessi (quindi anche reali e a maggior ragione anche interi).

  Il problema è quello di stuidare il condizionamento di
  $$y = x_1 + x_2,\ \ \ \ x_1, x_2 \in \R,\ \ \ \ x_1 + x_2 \neq 0$$
  Denotando con $\varepsilon_1$ e $\varepsilon_2$ gli errori relativi sui dati
  iniziali, ed assumendo che nessun nuovo errore venga introdotto nel calcoli,
  si ottiene:

  $$
  y(1 + \varepsilon_y) = x_1(1 + \varepsilon_1) + x_2(1 + \varepsilon_2)\
  = x_1 + x_2 + x_1\varepsilon_1 + x_2\varepsilon_2
  $$

  Si ricava:

  $$
  \vert \varepsilon_y\vert  \leq \frac{\vert x_1\vert  + \vert x_2\vert }{\vert x_1 + x_2\vert }\varepsilon_x\
  \equiv \kappa \varepsilon_x,\ \ \ \ \varepsilon_x\
  = max\{\vert \varepsilon_1\vert , \vert \varepsilon_2\vert \}
  $$

  </details>

- <details><summary>Dimostrazione A=LDL<sup>T</sup> (<b>X2</b>)</summary>

  ### Matrice simmetrica definita positiva

  Una matrice $A \in \R^{n \times n}$ è sdp se è simmetrica (cioè, $A = A^T$) e,
  per ogni $x \in \R^n$, $x \neq 0$, risulta $x^T A x > 0$

  #### Lemma 3.7

  Tutte le sottomatrici principali di una matrice sdp sono sdp.

  #### Lemma 3.8

  Una matrice sdp è nonsingolare.

  ### Teorema 3.4

  Se $A$ è sdp, allora è fattorizzabile $LU$.

  #### <span id="sdp-LU-fattorizzabile">Dimostrazione</span>

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

  </details>

- <details><summary>Precisione di macchina (<b>X2</b>)</summary>

  ### Teorema 1.3

  Il più piccolo ed il più grand (in valore assoluto), tra i numeri di
  macchina diversi da 0, sono rispettivamente dati da:

  $$
  \begin{align*}
    &r_1 = b^{-v},\\ &r_2 = (1 - b^{-m})b^{\varphi},\ \ \ \ \varphi = b^s - v
  \end{align*}
  $$

  Numeri di machina sono $L = [-r_2, -r_1]\cup\{0\}\cup[r_1, r_2]$

  ### Teorema 1.4

  Se $x \in L$, $x \neq 0$, allora
  $$fl(x) = x(1 + \varepsilon_x),\ \ \ \ \vert \varepsilon_x\vert  \leq u$$
  dove

  $$
  u = \begin{cases}
      b^{1-m},\ in\ caso\ di\ troncamento,\\
      \frac{1}{2} b^{1-m},\ in\ caso\ di\ arrotondamento
  \end{cases}
  $$

  ##### Dimostrazione

  La **_precisione di macchina_** è definita da quantità $u$ in Teorema 1.4

  </details>

- <details><summary>Condizionamento (di un problema) cap 1 (<b>X2</b>)</summary>

  $$
  \vert \varepsilon_y\vert  \approx \left\vert  f'(x)\frac{x}{y}\right\vert  \vert \varepsilon_x\vert \
  \equiv \kappa \vert \varepsilon_x\vert
  $$

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

  </details>

- <details><summary>Fattorizzazione QR (<b>X2</b>)</summary>

  $$Ax = b,\ \ \ \ A \in \R^{m \times n},\ \ \ \ m > n \equiv rank(A)$$

  ### Teorema 3.8

  Data la matrice $A$, esistono:

  - $Q \in \R^{m \times m}$, ortogonale
  - $\hat{R} \in \R^{n \times n}$, triangolare superiore e non singolare

  tali che:
  $$A = QR \equiv Q\binom{\hat{R}}{O}$$

  </details>

- <details><summary>I sistemi lineari in generale, cap 3 fino a dim unicità du A=LU (<b>X2</b>)</summary>
  </details>
- <details><summary>Metodo iterativo applicato ai sistemi lineari (<b>X2</b>)</summary>

  <a href="#google-pagerank">Google pagerank</a>

  <a href="#raggio-spettrale">Raggio spettrale</a>

  </details>

- <details><summary>Matrice ortogonale (<b>X2</b>)</summary>

  **_Matrice ortogonale_** è una matrice invertibile tale che la sua trasposta coincide con la sua inversa.

  </details>

- <details><summary>M matrici e matrici monotone in generale (<b>X2</b>)</summary>

  Le $M$-matrici sono particolari matrici $monotone$, in quanto, se $A$ è una
  $M$-matrice allora
  $$Ax \leq C\ \ \ \ \Rightarrow\ \ \ \ I \leq A^{-1}C.\ I \leq CA^{-1}$$
  dove, al solito, le diseguaglianze si intendono elemento per elemento.

  <a href="#google-pagerank">Google pagerank</a>

  <a href="#raggio-spettrale">Raggio spettrale</a>

  </details>

- <details><summary>Dimostrazione sdp &rArr; LU (<b>X2</b>)</summary>

  ### <a href="#sdp-LU-fattorizzabile">Dimostrazione Teorema 3.4</a>

  </details>

- <details><summary>def splitting regolare di matrice (<b>X2</b>)</summary>
  </details>
- <details><summary>metodo di newton (dim convergenza, molteplicità, Aitken) (<b>X2</b>)</summa
  ry></details>
- <details><summary>condizionamento del problema nei sistemi lineari</summary>
  </details>
- <details><summary>calcolare la norma uno e la norma infinito di una matrice</summary>
  </details>
- <details><summary>matrici triangolari con codice</summary>
  </details>
- <details><summary>condizionamento delle matrici</summary>
  </details>
- <details><summary>come si ottiene Va=f (risp: prodotto scalare tra gli elementi di V ed a)</summary>
  </details>
- <details><summary>Differenze tra LU e LDL<sup>T</sup></summary>

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

  $$
  g_i \equiv \frac{1}{a_{ii}^{(i)}}\
  (0, ..., 0, a_{i+1,i}^{(i)}, ..., a_{ni}^{(i)})^T
  $$

  $$
  L = \begin{pmatrix}
    1      & 0      & \cdots    & 0      \\
    g_{21} & 1      & \ddots    & \vdots \\
    \vdots & \ddots & \ddots    & 0      \\
    g_{n1} & \cdots & g_{n,n-1} & 1
  \end{pmatrix},\ \ \ \ \
  U = \begin{pmatrix}
    a_{11}^{(1)} & \cdots &                     & a_{1n}^{(1)}       \\
    0            & \ddots &                     & \vdots             \\
    \vdots       & \ddots & a_{n-1,n-1}^{(n)-1} & a_{n-1, n}^{(n-1)} \\
    0            & \cdots & 0                   & a_{nn}^{(n)}
  \end{pmatrix}
  $$

  ### LDL<sup>T</sup>

  Per ottenere $D$ basta prendere elementi diagonali di $U$ in $LU$

  $$
  L = \begin{pmatrix}
    1      & 0      & \cdots    & 0      \\
    g_{21} & 1      & \ddots    & \vdots \\
    \vdots & \ddots & \ddots    & 0      \\
    g_{n1} & \cdots & g_{n,n-1} & 1
  \end{pmatrix},\ \ \ \ \
  D = \begin{pmatrix}
    a_{11}^{(1)} & 0      & \cdots              & 0            \\
    0            & \ddots & \ddots              & \vdots       \\
    \vdots       & \ddots & a_{n-1,n-1}^{(n)-1} & 0            \\
    0            & \cdots & 0                   & a_{nn}^{(n)}
  \end{pmatrix}
  $$

  [Come funziona](https://yewtu.be/watch?v=8JdJoc3HMA8).

  </details>

- <details><summary>metodi iterativi per sistemi lineari</summary>
  </details>
- <details><summary>condizionamento cap 2</summary>
  </details>
- <details><summary>equazione retta tangente</summary>
  </details>
- <details><summary>def norma indotta su matrice</summary>
  </details>
- <details><summary>approssimazione polinomiale (inizio cap 4)</summary>
  </details>
- <details><summary>condizionamento del problema in generale</summary>
  </details>
- <details><summary>codice matrice ortogonale</summary>
  </details>
