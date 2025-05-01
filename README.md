# Calcolo delle Molteplicità Algebrica e Geometrica degli Autovalori

## 📘 Descrizione
Progeto sviluppato come parte integrante del corso di Calcolo Numerico al Secondo Anno di Informatica presso l'Uniersità degli Studi di Padova
Questo progetto implementa in MATLAB un insieme di funzioni per studiare la **diagonalizzabilità** di una matrice quadrata reale o complessa, calcolando:

- La **molteplicità geometrica** di un autovalore.
- La **molteplicità algebrica** di un autovalore.
- La funzione caratteristica `f_A(z)` e la sua derivata, utilizzando un approccio numericamente stabile.
- L'effettiva **diagonalizzabilità** confrontando le due molteplicità.

## ⚙️ Funzionalità principali

### 1. `multgeo(A, λ, toll)`
Calcola la **molteplicità geometrica** dell'autovalore `λ` della matrice `A`, come numero di righe approssimativamente nulle nella fattorizzazione LU di `(A − λI)`.  
**Input:**
- `A`: matrice quadrata (reale o complessa)
- `λ`: autovalore da analizzare
- `toll`: tolleranza per determinare nullità

**Output:**
- `k`: dimensione dell’autospazio di `λ`

---

### 2. `myobjective(z, A)`
Calcola la funzione caratteristica `f_A(z)` e la quantità stabile `g_A(z) = -f_A(z)/f_A'(z)`.  
**Input:**
- `z`: punto complesso
- `A`: matrice quadrata

**Output:**
- `f`: valore di `det(A − zI)`
- `g`: valore di `g_A(z)`

---

### 3. `multalg(A, l0, toll, it, maxit)`
Stima la **molteplicità algebrica** di un autovalore usando un metodo di Newton modificato e una stima dell'ordine di convergenza.  
**Input:**
- `A`: matrice quadrata
- `l0`: punto iniziale
- `toll`: tolleranza di arresto
- `it`: numero minimo di iterazioni iniziali
- `maxit`: massimo numero di iterazioni

**Output:**
- `l`: autovalore calcolato
- `m`: molteplicità stimata
- `flag`: 1 se successo, 0 se errore

---

## 🧪 Sperimentazione

Per testare l'algoritmo, si suggerisce di costruire matrici nella forma:
    **A = Q' * J * Q**

dove `Q` è ortogonale (es. tramite fattorizzazione QR) e `J` è una matrice a blocchi di Jordan.

Sono raccomandati test con:
- Autovalori ben separati vs. clusterizzati
- Molteplicità basse (1–2) e alte (3–7)

---

## 🧠 Fondamenti teorici

- **Molteplicità geometrica**: dimensione del kernel di `(A − λI)`, cioè `n − rank(A − λI)`
- **Molteplicità algebrica**: ordine della radice `λ` della funzione caratteristica `f_A(x)`
- **Diagonalizzabilità**: `A` è diagonalizzabile se, e solo se, `multgeo(λ) = multalg(λ)` per ogni autovalore `λ`

---

## 🛠️ Requisiti

- MATLAB (con supporto per la fattorizzazione LU e operazioni complesse)
- Nessun toolbox esterno richiesto
