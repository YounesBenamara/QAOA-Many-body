# Analyse Variationnelle du Modèle de Heisenberg XXZ via QAOA

Ce projet explore les limites des ansatz variationnels (QAOA) pour simuler des systèmes quantiques fortement corrélés, spécifiquement la chaîne de spins Heisenberg XXZ.

## 1. Le Système Physique : Chaîne de Spins 1D (Modèle XXZ)

On considère une chaîne linéaire de $N$ sites (qubits), où chaque site $j$ porte un spin $1/2$. Le système est régi par l'Hamiltonien de Heisenberg XXZ.

**L'Hamiltonien $H$ :**
Pour des conditions aux limites ouvertes (OBC), l'opérateur s'écrit :

$$
H_{XXZ}(\Delta) = J \sum_{j=0}^{N-2} \left( \sigma^x_j \sigma^x_{j+1} + \sigma^y_j \sigma^y_{j+1} + \Delta \sigma^z_j \sigma^z_{j+1} \right)
$$

Où :
* $J > 0$ est la constante de couplage (Antiferromagnétique). On fixe $J=1$.
* $\sigma^x, \sigma^y, \sigma^z$ sont les matrices de Pauli.
* $\Delta$ est le **paramètre d'anisotropie** (paramètre de contrôle).

## 2. Décomposition des Termes et Compétition Physique

La difficulté réside dans la non-commutativité des termes :

* **Terme d'Ising (Diagonal) :** $\Delta \sum \sigma^z_j \sigma^z_{j+1}$
    * Favorise l'alignement classique (ordre de Néel $\uparrow \downarrow \uparrow \downarrow$).
    * Équivalent au problème MaxCut.
* **Terme d'Échange (Hors-Diagonal) :** $\sum (\sigma^x_j \sigma^x_{j+1} + \sigma^y_j \sigma^y_{j+1})$
    * Provoque le "saut" des excitations ($|01\rangle \leftrightarrow |10\rangle$).
    * Génère l'intrication quantique et la superposition.

**Problème :** Compétition entre l'ordre classique imposé par $\Delta$ et les fluctuations quantiques ($XX+YY$).

## 3. Objectif de l'Optimisation

L'objectif est d'approximer l'énergie du fondamental $E_0(\Delta)$.

**Problème variationnel :**
Trouver les paramètres $\vec{\theta}$ minimisant l'espérance de l'Hamiltonien :

$$
E_{QAOA} = \min_{\vec{\theta}} \langle \psi(\vec{\theta}) | H_{XXZ}(\Delta) | \psi(\vec{\theta}) \rangle
$$

Cette valeur est comparée à la "Vérité Terrain" $E_{exact}$ obtenue par diagonalisation exacte.

## 4. L'Ansatz QAOA (Le Circuit)

L'état d'essai $|\psi(\vec{\gamma}, \vec{\beta})\rangle$ est généré par l'application alternée de deux unitaires sur l'état initial $|+\rangle^{\otimes N}$ :

1.  **Unitaire de Coût ($U_C$) :** Basé uniquement sur la partie Ising.
    $$U_C(\gamma_k) = e^{-i \gamma_k \sum \sigma^z_j \sigma^z_{j+1}}$$
2.  **Unitaire de Mélange ($U_M$) :** Mélangeur transverse standard.
    $$U_M(\beta_k) = e^{-i \beta_k \sum \sigma^x_j}$$

Le circuit de profondeur $p$ est :

$$
|\psi(\vec{\gamma}, \vec{\beta})\rangle = \prod_{k=1}^p \left( U_M(\beta_k) U_C(\gamma_k) \right) |+\rangle^{\otimes N}
$$

## 5. Méthodologie

1.  **Input :** Choix de $\Delta$ (Diagramme de phase).
2.  **Modèle :** Construction de la matrice $H_{XXZ}$ (Opérateur de mesure).
3.  **Boucle VQE :** Optimisation des paramètres $\gamma, \beta$ via un optimiseur classique (COBYLA/SPSA).
4.  **Analyse :** Calcul du Ratio d'Approximation $r = E_{QAOA} / E_{exact}$.
