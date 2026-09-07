# 🚀 
# 🚀 Numerical Physics Engine : FEA Beam & 2D CFD Solver

Ce dépôt regroupe une suite d'outils de simulation numérique haute performance développés en **C++17**, couvrant à la fois la **mécanique des structures (1D Éléments Finis)** et la **mécanique des fluides (2D Volumes Finis / Navier-Stokes & LES)**.

---

## 🌊 1. Simulation Numérique de Fluides (2D CFD / LES)

Ce module (`drone_les_sim/`) implémente un solveur complet des équations de **Navier-Stokes incompressibles** en 2D avec modélisation de la turbulence (LES) et méthode de projection.

### 🛠️ Choix Numériques & Physique
* **Méthode des Volumes Finis (FVM)** sur grille colocalisée.
* **Algorithme de Projection (Chorin-Temam)** : découplage vitesse-pression par un pas de prédiction (advection-viscosité) et un pas de correction (gradient de pression).
* **Équation de Poisson pour la Pression** : résolue par méthode itérative **SOR** (*Successive Over-Relaxation*) avec conditions aux limites de Neumann ($\frac{\partial P}{\partial n} = 0$).
* **Turbulence (LES)** : sous-maille modélisée par viscosité turbulente.
* **Gestion Mémoire Optimisée** : utilisation de la classe `Field2D` encapsulant un `std::vector<double>` 1D contigu en mémoire (*row-major*) pour maximiser le cache CPU et éviter les copies grâce à `std::swap`.

---

## 🏗️ 2. Simulation de Flexion de Poutre (1D Éléments Finis)

Ce module étudie la déformation d'une poutre soumise à des charges externes en comparant la précision des approximations et l'efficacité des résolveurs matriciels.

### 🛠️ Modélisation
* **P1 (Lagrange Linéaire)** : approximation par segments droits, générant des matrices **tridiagonales**.
* **P2 (Lagrange Quadratique)** : approximation par paraboles, offrant une précision d'ordre supérieur et générant des matrices **pentadiagonales**.

### 💡 Structure de Données : Vecteurs Séparés
Plutôt que d'utiliser une matrice dense (`N x N`), le projet stocke uniquement les diagonales utiles sous forme de vecteurs indépendants (`d_centrale`, `d_inf1`, `d_sup1`, `d_inf2`, `d_sup2`).
* **Optimisation Mémoire** : pour $N = 10\,000$, la mémoire passe de ~800 Mo à ~0,4 Mo.
* **Complexité** : réduction des calculs de $\mathcal{O}(N^3)$ à une complexité linéaire **$\mathcal{O}(N)$**.

### 💻 Solveurs Matriciels Implémentés
* **Directs** : Algorithme de Thomas (LU) et Décomposition de Cholesky.
* **Itératifs** : Jacobi et Gauss-Seidel.

---

## 📁 Architecture du Projet

```text
.
├── 🏗️ Poutre_FEA/           # Code de simulation de la poutre (Éléments Finis 1D)
│   ├── Poutre.cpp / .hpp
│   ├── Methode.cpp / .hpp
│   └── Resolution.cpp / .hpp
│
├── 🌊 drone_les_sim/        # Solveur de dynamique des fluides (CFD / LES 2D)
│   ├── include/
│   │   ├── core/           # Field2D, Mesh
│   │   └── solvers/        # PoissonSolver, NavStokesSolver, ProjectionMethod
│   └── src/                # Implémentations .cpp et main.cpp
│
└── main.cpp                # Point d'entrée pour les tests FEA









Simulation de Flexion de Poutre : Analyse Numérique & Éléments Finis

Ce projet implémente une chaîne complète de simulation numérique pour l'étude de la déformation d'une poutre soumise à des charges externes. Il permet de comparer la précision des modèles (**P1 vs P2**) et l'efficacité des solveurs (**Directs vs Itératifs**).

## 🛠 Méthodes de Modélisation

Le projet repose sur la discrétisation de l'équation différentielle de la flexion par la méthode des **Éléments Finis** :
* **P1 (Lagrange Linéaire)** : Approximation par segments droits, générant des matrices **tridiagonales**.
* **P2 (Lagrange Quadratique)** : Approximation par paraboles, offrant une précision d'ordre supérieur et générant des matrices **pentadiagonales**.

## 🏗 Choix de Structure de Données : Vecteurs Séparés

Plutôt que d'utiliser une matrice dense (`n x n`), ce projet stocke uniquement les diagonales utiles sous forme de vecteurs indépendants :
* `std::vector<double> d_centrale` : Diagonale principale (taille $n$).
* `std::vector<double> d_inf1`, `d_sup1` : Diagonales de voisinage direct (taille $n-1$).
* `std::vector<double> d_inf2`, `d_sup2` : Diagonales de voisinage étendu pour P2 (taille $n-2$).

### Pourquoi ce choix ?
1.  **Optimisation Mémoire** : Pour $n=10,000$, une matrice dense utilise 100 millions de doubles (~800 Mo). Notre structure n'en utilise que ~50,000 (~0.4 Mo).
2.  **Performance Cache** : L'accès aux données est linéaire et contigu, maximisant l'efficacité du processeur.
3.  **Complexité** : Les algorithmes passent d'une complexité $O(n^3)$ à une complexité **linéaire $O(n)$**.

## 💻 Solveurs Implémentés

### Méthodes Directes
* **Algorithme de Thomas (LU)** : Une variante simplifiée de l'élimination de Gauss pour les systèmes rubanés.
* **Décomposition de Cholesky** : Version optimisée pour les matrices symétriques définies positives, idéale pour les problèmes de structures stables.

### Méthodes Itératives
* **Jacobi** : Méthode de point fixe utilisant le vecteur de l'itération précédente pour calculer la nouvelle approximation.
* **Gauss-Seidel** : Utilise les valeurs déjà mises à jour au cours de l'itération pour accélérer la convergence.

## 🚀 Utilisation

### Compilation
```bash
g++ -Wall -O3 -o simulation main.cpp Poutre.cpp Methode.cpp Resolution.cpp