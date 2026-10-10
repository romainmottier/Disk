# Stratégie de parallélisation CPU et GPU de `unfitted_HHO`

Ce document fixe une stratégie pour paralléliser l'application `unfitted_HHO`
(méthodes HHO non conformes stabilisées par extension polynomiale ou par
agglomération), sur CPU (OpenMP, MPI) et sur GPU (CUDA, HIP, Kokkos,
OpenACC). Il part de mesures faites sur le code actuel. Il indique, pour chaque
étape, les endroits du code à modifier, le gain attendu et la manière de valider
le résultat.

## 1. État des lieux

### 1.1 Parallélisme actuel

Le code est séquentiel, à une exception près : la résolution du système global
par Pardiso (MKL), qui est multithreadée (voir le `README.md` de
l'application, réglage *solver*, pour l'édition de liens de la couche de threads
de MKL). Il n'y a
aucune directive OpenMP, ni dans l'application, ni dans le code cut de diskpp
(`cut_mesh.hpp`, `cut_geometry.hpp`, `cut_bases.hpp`).

Le build `Release` de diskpp compile pour un x86-64 générique (SSE2) : Eigen
n'utilise ni l'AVX2 ni l'AVX-512. Le type de build `ReleaseNative`
(`-O3 -march=native`) existe déjà dans le `CMakeLists.txt` racine.

### 1.2 Où part le temps

Temps par étape au niveau le plus fin, k = 3, extension polynomiale (`-p 0`),
test (56), réglages par défaut. Ces runs ont été faits avec 3 threads chacun et
8 runs simultanés sur 12 cœurs : les valeurs absolues sont pessimistes, mais la
répartition est représentative.

| cas | cellules | maillage + géométrie | assemblage (dont cellules coupées) | résolution | erreurs |
|---|---|---|---|---|---|
| tétraèdres, sphère | 48 000 | 1,7 s (1 %) | 88 s (24 %), dont 81 s | 232 s (63 %) | 46 s (12 %) |
| tétraèdres, cube arrondi (raffiné) | 48 812 | 48 s (10 %) | 115 s (24 %), dont 107 s | 219 s (46 %) | 96 s (20 %) |
| hexaèdres, sphère | 8 000 | 0,6 s (1 %) | 25 s (46 %), dont 25 s | 19 s (34 %) | 10 s (19 %) |
| prismes, sphère | 5 600 | 0,5 s (1 %) | 27 s (36 %), dont 21 s | 27 s (36 %) | 20 s (27 %) |
| quadrangles (2D), cercle | 6 400 | 0,06 s | 0,22 s | 0,24 s | 0,06 s |

Ordres de grandeur pour k = 3 au niveau 2 :

* hexaèdres : 824 cellules coupées, environ 250 000 inconnues de face ;
* tétraèdres : 3 804 cellules coupées, environ 1 million d'inconnues de face,
  13,9 Go de mémoire (essentiellement la factorisation) ;
* environ 3 000 points de quadrature par partie de cellule coupée.

Conclusions :

* la 2D coûte moins d'une seconde : la parallélisation vise la 3D ;
* en 3D, tout ce qui n'est pas la résolution est du **travail local par
  cellule** (géométrie de coupe, opérateurs locaux, erreurs), parfaitement
  parallèle. Plus de 90 % de l'assemblage se passe dans les cellules coupées,
  car les opérateurs des cellules non coupées sont calculés une fois par forme
  (`-u 1`) ;
* la résolution représente l'autre moitié du temps, et sa mémoire limite la
  taille des problèmes 3D (factorisation directe).

### 1.3 Machine de développement

* CPU AMD Ryzen AI 9 HX PRO 370 (Zen 5), 12 cœurs, 24 threads, AVX-512,
  86 Go de mémoire, un seul nœud NUMA ;
* GPU intégré Radeon 890M (gfx1150), qui partage la mémoire : ROCm n'est pas
  installé (seulement le runtime HSA) et gfx1150 n'est pas officiellement pris
  en charge. Pas de GPU NVIDIA : CUDA n'est pas utilisable ici ;
* installés : OpenMPI 4.1.6, MUMPS 5.6 (séquentiel et parallèle), METIS,
  Scotch. Absents : PETSc, hypre, Kokkos, MAGMA.

Le développement CPU (OpenMP, MPI) peut se faire sur cette machine ; les calculs
GPU demandent une machine équipée (NVIDIA A100/H100 ou AMD MI250/MI300).

## 2. Principes

1. **Séparer la géométrie de l'algèbre.**
   * La géométrie de coupe est irrégulière : recherche de racines, découpages
     récursifs des polygones et des polyèdres, repli par subdivision. Elle
     reste sur CPU.
   * L'algèbre locale est régulière : évaluation des bases aux points de
     quadrature, produits B₁ᵀ diag(w) B₂, petits systèmes denses (Cholesky,
     LU), condensation statique. Elle se prête au traitement par lots, donc au
     GPU.
2. **Exploiter le maillage structuré.** Les cellules non coupées d'une même
   forme ont le même opérateur local (clé `uncut_cell_key` dans
   `operators.hpp`). Sur GPU, appliquer leur contribution revient à appliquer
   une seule petite matrice dense à toutes ces cellules : un seul gros produit
   matriciel.
3. **Unité de travail : le patch.** Un patch est une cellule et les cellules mal
   coupées qu'elle stabilise (composante connexe de la relation S ~ N(S), voir
   `doc/static_condensation.pdf`), ou une cellule agglomérée. Les patchs sont
   indépendants : c'est la granularité naturelle pour OpenMP, pour la
   répartition MPI et pour les lots GPU.
4. **Un seul code pour CPU et GPU** quand c'est possible (Kokkos ou PETSc),
   pour ne pas maintenir deux versions.

## 3. Feuille de route

### Phase 0 : gains immédiats (environ 1 jour)

* Compiler en `ReleaseNative` (`-march=native`) : vectorisation AVX-512 des
  produits Eigen des opérateurs locaux et du calcul des erreurs.
* Construire le maillage coupé une fois par niveau. Aujourd'hui,
  `pre::msh_generation` est appelé pour chaque couple (k, niveau) dans
  `Elliptic_Interface` (`unfitted_HHO_run.hpp`) : la géométrie de coupe, la
  détection des cellules non coupées une fois et le raffinement local sont
  refaits pour chaque k, alors qu'ils ne dépendent pas de k (sauf
  l'interface polynomiale `-g -2`, dont le degré vaut k+2). Il suffit
  d'inverser les boucles ou de garder le maillage de chaque niveau.
* Mesurer chaque étape avec les minuteries existantes (`wall_timer`).

Gain attendu : un facteur 1,5 à 3 sur l'assemblage et les erreurs, à mesurer ;
jusqu'à un facteur 4 sur le temps de maillage des runs raffinés.

### Phase 1 : OpenMP (environ 1 semaine)

C'est l'étape au meilleur rapport gain/effort.

**Boucles à paralléliser** (`#pragma omp parallel for schedule(dynamic)`) :

| étape | fonction | remarque |
|---|---|---|
| géométrie des faces | `detect_cut_faces` (`cut_mesh.hpp`) | calculer d'abord les intersections des arêtes (`edge_intersections`) dans une passe parallèle sur les arêtes, puis les faces en lecture seule |
| géométrie des cellules | `detect_cut_cells` | indépendante par cellule |
| raffinement local | `cut_once_violations` | indépendante par face puis par cellule |
| appariement | `make_polynomial_extension`, `find_good_neighbor` | choix indépendants par cellule mal coupée, puis construction des listes de dépendants |
| opérateurs locaux et condensation | boucle sur les patchs de `Elliptic_Interface` | indépendante par patch |
| erreurs | `postprocessor::compute_errors`, `compute_errors_agglo` | réduction (`reduction(+:...)`) |

**États partagés à rendre sûrs entre threads :**

* caches statiques des règles de quadrature et des exposants :
  `gauss_legendre_01`, `reference_triangle_rule`, `reference_tetrahedron_rule`
  (`cut_geometry.hpp`), `cut_monomial_exponents` (`cut_bases.hpp`). Les remplir
  au démarrage pour tous les degrés utilisés, ou les protéger (verrou,
  `std::call_once` par degré) ;
* cache LRU mutable des quadratures coupées (`cut_mesh::quadrature_cache`) et
  caches des repères de bases (`cut_basis_cache` : `frames`, `uncut_frames`,
  `face_frames`, `uncut_face_frames`). Pendant l'assemblage, les remplacer par
  des tableaux indexés par (cellule, côté), calculés dans une passe parallèle
  préalable puis lus seulement. C'est plus simple et plus sûr qu'un cache
  concurrent ;
* cache des opérateurs des cellules non coupées (`uncut_cache` de
  `unfitted_hho_parameters`) : calculer d'abord une cellule de chaque forme,
  puis lire seulement.

**Assemblage global.** Les assembleurs (`unfitted_elliptic_interface_assembler`
dans `cuthho_assemblers.hpp`, `local_condensation` dans
`static_condensation.hpp`) accumulent des triplets dans un vecteur unique :

* version simple : un vecteur de triplets par thread, concaténés avant
  `setFromTriplets` (mémoire multipliée par le nombre de threads pour les
  triplets seulement) ;
* version pour les gros maillages : structure CSR des inconnues de face
  calculée une fois, puis remplissage parallèle sans atomiques grâce à un
  coloriage des patchs (deux patchs de même couleur ne partagent aucune face).

**Gain attendu** : un facteur 8 à 10 sur la géométrie, l'assemblage et les
erreurs avec 12 cœurs, soit un facteur 2,5 à 4 sur le temps total en 3D. La
résolution devient alors l'étape dominante.

### Phase 2 : le solveur (1 à 2 semaines)

Le système condensé (`-c 1`) ne porte que les inconnues de face et il est
symétrique défini positif.

* Sur un nœud : Pardiso sur tous les cœurs physiques, ou MUMPS (installé), qui
  sert aussi en MPI.
* En 3D plus fin, la factorisation directe sature la mémoire : 13,9 Go déjà pour
  les tétraèdres au niveau 2 avec k = 3, et environ 8 millions d'inconnues de
  face au niveau 3, avec un remplissage qui croît comme N^{4/3}. Il faut alors
  un gradient conjugué préconditionné par l'une des options suivantes :
  * multigrille algébrique : hypre BoomerAMG, ou GAMG de PETSc ;
  * p-multigrille (de k vers 0) suivie d'une multigrille algébrique sur l'ordre
    bas, adaptée aux méthodes d'ordre élevé ;
  * BDDC (`PCBDDC` de PETSc), naturel pour les méthodes hybrides (les
    inconnues de face forment l'interface entre sous-domaines) et robuste au
    contraste des coefficients avec la mise à l'échelle « deluxe ».
* Passer par **PETSc** : le même code sert ensuite pour MPI (phase 3) et pour
  le GPU (phase 4). La fonction `solve_linear_system` (`unfitted_HHO_run.hpp`)
  devient une interface vers un objet `KSP`.

### Phase 3 : MPI (2 à 4 semaines)

* **Partition** du maillage de fond structuré par blocs, ou par METIS ou Scotch
  (installés). Chaque patch appartient au rang de sa cellule stabilisatrice,
  pour ne jamais couper un patch.
* **Couche fantôme** de 2 cellules : l'appariement utilise les voisins par face
  et par sommet, et un patch peut déborder d'une cellule.
* **Raffinement local et équilibrage 2:1** (`refine_until_cut_once`) sur la
  description complète du maillage de fond, peu coûteuse, avant la partition.
* **Inconnues de face** numérotées globalement (propriétaire et fantômes) ;
  assemblage dans une matrice PETSc distribuée (`MatMPIAIJ`, ou `MatIS` pour
  BDDC) ; erreurs par `MPI_Allreduce`.
* **Hybride MPI et OpenMP** : un rang par complexe de cœurs, des threads à
  l'intérieur.

### Phase 4 : GPU (1 à 2 mois)

**Solveur.** Les backends GPU de PETSc (CUDA, HIP ou Kokkos), avec un
gradient conjugué et GAMG ou hypre sur GPU, ou bien cuDSS (direct) et AMGX sur
NVIDIA. Si la phase 2 est faite avec PETSc, ce point se réduit à des options
d'exécution.

**Opérateur sans assemblage (matrix-free).** Exploiter le maillage structuré :

* toutes les cellules non coupées d'une même forme partagent la même matrice
  condensée : l'appliquer à toutes ces cellules est un seul produit matriciel
  (matrice locale × matrice des inconnues de face de toutes les cellules) ;
* les cellules coupées gardent leurs propres matrices, appliquées par lots ;
* le préconditionneur est une p-multigrille ou une multigrille algébrique sur
  l'ordre bas.

Cette approche demande très peu de mémoire et atteint un débit très élevé sur
GPU.

**Algèbre locale des cellules coupées.** La géométrie reste sur CPU (OpenMP) ;
le CPU transmet par cellule les points et poids de quadrature et les repères des
bases. Le GPU traite les cellules par lots :

* évaluation des bases de Legendre aux points de quadrature (un thread par
  point) ;
* produits B₁ᵀ diag(w) B₂ (GEMM par lots, de tailles variables) ;
* petits systèmes denses : Cholesky ou LU de 35×35 à 300×300 pour k = 3
  (reconstruction du gradient, stabilisations, condensation du patch) ;
* renvoi des contributions condensées sur les inconnues de face.

On pipeline les deux : le CPU calcule la géométrie du lot suivant pendant que le
GPU calcule le lot courant.

Ordre de grandeur : pour k = 3, 824 cellules coupées à environ 3 000 points par
partie représentent de l'ordre de 10² GFLOP d'algèbre locale. Cela fait environ
1 s sur les 12 cœurs Zen 5 en AVX-512, et 10 à 50 ms sur un A100, un H100 ou un
MI300.

**Erreurs.** Mêmes noyaux (évaluation des bases aux points de quadrature),
suivis d'une réduction.

**Choix des outils.**

| outil | intérêt | limites |
|---|---|---|
| Kokkos et KokkosKernels | un seul source pour CPU (OpenMP), NVIDIA (CUDA), AMD (HIP), Intel (SYCL) ; noyaux denses par lots | dépendance supplémentaire |
| CUDA, avec cuBLAS ou MAGMA par lots | performance maximale sur NVIDIA | NVIDIA seulement, inutilisable sur la machine de développement |
| HIP, avec rocBLAS ou MAGMA | AMD (MI250, MI300) | gfx1150 non pris en charge officiellement |
| OpenMP `target` (déport) | directives, portable (GCC, Clang, compilateurs AMD et NVIDIA) | performances inégales selon les compilateurs |
| OpenACC | directives simples pour les boucles | vraiment praticable seulement avec NVHPC, donc NVIDIA ; prise en charge limitée dans GCC |

Recommandation : **Kokkos** pour les noyaux maison (algèbre locale, erreurs) et
**PETSc** pour les solveurs. Le développement se fait sur cette machine avec le
backend OpenMP de Kokkos ; les calculs GPU sur un cluster.

## 4. Validation

* Les sommes parallèles et l'ordre d'assemblage changent l'arrondi : les
  résultats ne sont plus identiques bit à bit. On compare à la version
  séquentielle avec une tolérance relative (par exemple 1e-10 sur les erreurs),
  comme pour les comparaisons de non-régression faites jusqu'ici.
* On relance la campagne de convergence sur maillages structurés (6 maillages,
  interfaces régulières, tests (56), (59), 5, et (57) à (58b) pour le cercle et
  la sphère, deux stabilisations). Elle vérifie les ordres de convergence par
  pente des moindres carrés sur les niveaux 1 et suivants.
* On mesure l'accélération de chaque étape (géométrie, assemblage, résolution,
  erreurs) avec les minuteries existantes, pour 1, 2, 4, 8 et 12 cœurs.

## 5. Ordre recommandé

1. Phase 0 et phase 1 : quelques jours de travail pour un facteur 3 à 4 sur le
   temps total en 3D sur la machine de développement.
2. Phases 2 et 3 avec PETSc : résolution des gros problèmes 3D et passage à
   plusieurs nœuds.
3. Phase 4 : GPU, d'abord le solveur (options PETSc), puis l'opérateur sans
   assemblage et l'algèbre locale des cellules coupées. C'est la phase la plus
   lourde, et elle demande une machine équipée.
