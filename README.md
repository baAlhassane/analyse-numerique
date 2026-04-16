# 🚀 Simulation de Flexion de Poutre : Analyse Numérique & Éléments Finis

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