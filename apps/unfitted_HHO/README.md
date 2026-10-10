# unfitted_HHO

Unfitted HHO methods for the elliptic interface problem, ported from ProtoN.
This app reproduces the numerical results of Section 5 of *Unfitted hybrid
high-order methods stabilized by polynomial extension for elliptic interface
problems* (E. Burman, A. Ern, R. Mottier), and runs the same methods on the
other meshes of diskpp, in 2D and 3D (see *Meshes and dimension*).

## Prototypes

| `-p` | method |
|------|--------|
| `0`  | stabilization of the ill-cut cells by polynomial extension |
| `1`  | stabilization by cell agglomeration ([12] in the paper) |

## Options

| option | meaning | default |
|--------|---------|---------|
| `-M` | mesh: 0 quadrangles (paper), 1 triangles, 2 polygons (hexagons), 3 hexahedra, 4 tetrahedra, 5 polyhedra (hexagonal prisms) | 0 |
| `-k` | face degree k (cell degree k+1); degrees `-K`..`k` are run | 0 |
| `-l` | refinement level (quadrangles: mesh of (10·2^l)² cells, see below for the other meshes); levels `-L`..`l` are run | 0 |
| `-K`, `-L` | first degree / first level of the loops | 0 |
| `-i` | interface: 0 circle (sphere in 3D), 1 flower (55b) (3D: sphere with 4 petals), 2 square (2D), 3 flower of ProtoN (− c cos nθ, 2D), 4 ellipse (ellipsoid), 5 rounded square (cube), 6 annulus (torus); see *Interfaces* | 0 |
| `-R` | radius of the circle / flower / sphere | 1/3 |
| `-d` | square interface shifted by δ = 0.5·10^-d (coarsest mesh only: on finer meshes its sides lie on mesh lines) | 2 |
| `-t` | test case: 0 (56), 1 (57), 2 (58a), 3 (58b), 4 (59), 5 contrast for any interface: u_i = ∏ sin(π x_d) / κ_i (1, 2, 3 need the circle or the sphere) | 0 |
| `-m` | contrast κ₂ = 10^m κ₁ | 0 |
| `-g` | interface in the cut cells: -1 exact (root of the level set), 0 polyline of `-r` (paper), q ≥ 2 polynomial of degree q, -2 polynomial of degree k+2 (2D) | -1 |
| `-r` | polyline of 2^r segments in each cut cell (`-g 0`, and fallback of the curved interfaces) | 8 |
| `-a` | flagging parameter ϑ of the ill-cut cells (0: no flagged cell) | 0.3 |
| `-b` | cell bases: 0 scaled monomials of the cells, 1 scaled monomials of the subcells (paper), 2 orthonormal | 2 |
| `-P` | neighbor of an ill-cut cell (pairing, agglomeration): 0 as in ProtoN (paper), 1 largest part on the same side | 1 |
| `-e` | η_N of s^N: x ≥ 0 fixed value, x < 0 \|x\|·(k+1)² (`-e -1`: (k+1)² of ProtoN, used by `run_paper.sh`) | -20, i.e. 20·(k+1)² |
| `-c` | static condensation of the cell unknowns: 0 none, 1 local (patches), 2 algebraic | 1 |
| `-s` | solver: 0 SimplicialLDLT, 1 conjugate gradient (Jacobi), 2 Pardiso (Intel MKL) | 2 with MKL, else 0 |
| `-u` | 1: local operators of the uncut cells computed once per shape | 1 |
| `-x` | local refinement of the cells that the interface does not cut once, at most x times, if there are any (polytopal meshes, see *Local refinement*); 0: none | 3 |
| `-q` | 1: compute the condition number, also after a diagonal scaling (dense, coarse meshes only) | 0 |
| `-z` | 1: write the block sparsity pattern (`sparsity_k*_l*.txt`) | 0 |
| `-f` | 1: debug outputs (silo, okc arrows, `cells_l*.txt`) | 0 |

The errors are written to `solution_error_file.txt` (ProtoN format, plotted by
the generated `solution_error_file.py`). The `H1-norm error` and `L2-norm error`
lines extend the polynomial of the paired cell to the ill-cut subcells, as in
ProtoN and in the figures of the paper; the lines `... with own unknowns on
ill-cut subcells` use the own unknowns of these subcells (norm of Theorem 4.6).

## Reproducing the paper

```
scripts/run_paper.sh <path/to/unfitted_HHO> results
python3 scripts/plot_paper_results.py results
```

writes `fig5.png` ... `fig12.png` in `results`, the cell-agglomeration curves
being dashed. The script uses the parameters of ProtoN, with which the figures
of the paper were computed: η_N = (k+1)² (`-e -1`), scaled monomial bases
centered on the subcells (`-b 1`), the pairing of ProtoN (`-P 0`), the
polyline interface (`-g 0`), the full system solved by SimplicialLDLT
(`-c 0 -s 0 -u 0`), and the meshes of the paper, without local refinement
(`-x 0`).

## Default discretization

The defaults of the app improve the conditioning and the robustness of the
method of the paper at high degree, without changing its formulation:

* **cell bases** (`-b 2`): on each region (cell, subcell, ill-cut subcell
  merged with its stabilizing cell, agglomerated cell), tensor Legendre
  polynomials in the bounding box of the region in its principal axes of
  inertia, orthonormalized (Cholesky, twice) for the mean L2 product of the
  region; on the uncut cells they are the Legendre polynomials of the cell.
  The discrete solution is the same as with the monomials (up to round-off),
  but the worst condition number of Fig. 12 drops, e.g. for k = 3 from 6e11
  to 6e6 (from 8e6 to 3e5 after a diagonal scaling), and for k = 5 from 4e18
  to 8e9; the gradient of an ill-cut subcell is assembled directly from its
  stiffness matrix;
* **pairing** (`-P 1`): ProtoN chooses the neighbor of an ill-cut cell with the
  smallest part on the other side, but takes the whole measure of an uncut
  neighbor, so that uncut neighbors are never chosen when a cut one exists. An
  ill-cut strip along a face was thus paired with a thin cut neighbor instead
  of the uncut cell across the face. Preferring the largest part on the same
  side removes the remaining peaks of Fig. 12 (k = 3: worst condition number
  6e6 -> 6e5) and improves the errors at high degree (k = 5, finest meshes:
  factor 50 to 200), the errors for k <= 2 being unchanged;
* **η_N** = 20·(k+1)² (see below);
* **exact interface** (`-g -1`): in each cut cell, the interface is the
  zero level set itself, as the graph of a height above its chord p0-p1,
  computed by a root of the level set along the normal of the chord at each
  quadrature point (Newton iterations safeguarded by bisection). Each subcell
  is integrated by sectors from a point with respect to which it is
  star-shaped (triangles on the straight edges, curved sector on the
  interface): all the quadrature points lie in the subcell and all the
  weights are positive. A cell where no such point is found is split by the
  line through the middle of its interface, recursively. The geometry is
  exact up to round-off (area of the disk and length of the circle: errors
  1e-15 on every mesh, against 1e-8 with the polyline, r = 8): with the
  polyline, the errors of the tests with a contrast or jumps (57)-(59)
  stagnate for k ≥ 3 (e.g. L2 error 2e-12 for (57), k = 3, level 4), whereas
  they converge at the optimal rate with the curved interface (2e-13, and
  down to round-off for k = 4, 5), and the cut cells need about 20 times
  fewer quadrature points. `-g -2` replaces the exact height by its
  polynomial interpolant of degree k+2 at the Chebyshev-Lobatto nodes (the
  previous default, same errors up to round-off). Cells where the interface
  is not a graph above its chord, e.g. where it crosses a face twice, keep
  the polyline (warning);
* **static condensation** (`-c 1`): the cell unknowns are eliminated patch by
  patch before the assembly, a patch being a cell together with the ill-cut
  cells it stabilizes (connected components of S ~ N(S)) or an agglomerated
  cell; see `doc/static_condensation.pdf`. The condensed system is much better
  conditioned (Fig. 12 configurations, level 0: condition number 266 for k = 3
  instead of 1.1e5 for the full system) and 3 times smaller for k = 3;
* **solver** (`-s 2`): Pardiso (multithreaded) when Intel MKL is
  available. The Intel threading layer of MKL linked by diskpp needs the
  OpenMP runtime of Intel, but `libiomp5` is the `libomp` of LLVM on Debian
  and Ubuntu: the multithreaded factorization failed. With GCC, the
  `CMakeLists.txt` of the app links the GNU threading layer of MKL and
  `libgomp` first, which take precedence (the root fix would be in
  `cmake/FindMKL.cmake`, which links the layered libraries and `mkl_rt`
  together). Should the factorization fail anyway, it is retried with one
  thread, then SimplicialLDLT is used;
* **cache of the uncut cells** (`-u 1`): the local matrices of the uncut
  cells without dependent cells only depend on the shape of the cell (up to a
  translation), on the numbering of its faces and on κ; they are computed
  once per shape (key: positions of the vertices of the cell and of its
  faces relative to its first vertex). On the structured meshes (Cartesian
  grids, triangles and tetrahedra of the Kuhn subdivision, prisms), a few
  shapes cover all the uncut cells. Since the same local matrix is used
  everywhere, its round-off errors add up coherently: the round-off floor of
  the errors can be 2 to 3 times higher (around 1e-13 for the L2 error,
  k = 3, level 4); `-u 0` computes them in every cell;
* **local operators**: the integrals are computed by matrix products of the
  values of the bases at all the quadrature points (B1ᵀ diag(w) B2), the bases
  being evaluated by blocks; the gradient reconstruction uses the mass matrix
  of the scalar basis for each component;
* **local refinement** (`-x 3`): if some cells of the meshes of the run are
  not cut once by the interface, these cells are refined (see *Local
  refinement*); otherwise the meshes are not changed.

With these defaults, the convergence is optimal up to k = 4 on the finest
meshes tested (flower, k = 4: H1 errors 3.9e-6, 1.1e-7, 3.2e-9, 1.0e-10 on the
levels 0 to 3); for k = 5, the errors stagnate around 1e-10 (round-off).

Computing time for k = 3 on the level 4 (circle, poly. extension, 24 cores):

| mesh, settings | assembly (cut cells) | solve | errors | total |
|----------------|----------------------|-------|--------|-------|
| 160² quadrangles, paper (`run_paper.sh`) | 2.7 s (1.8 s) | 4.3 s | 1.0 s | 8.0 s |
| 160² quadrangles, defaults with the polyline (`-g 0`) | 2.0 s (1.8 s) | 0.5 s | 0.9 s | 3.5 s |
| 160² quadrangles, defaults | 0.35 s (0.12 s) | 0.45 s | 0.12 s | 1.05 s |
| 51200 triangles (`-M 1`), defaults | 0.96 s (0.19 s) | 0.57 s | 0.62 s | 2.4 s |
| 16576 hexagons (`-M 2`), defaults | 0.72 s (0.12 s) | 0.67 s | 0.59 s | 2.1 s |

(With the previous version of the code, the defaults took 1.9 s: the solve of
Pardiso with one thread, 1.2 s, and the assembly of the cut cells by rank-one
updates, 0.27 s.)

A quadrature of the polynomial integrals by integration on the boundary of
the subcells (Antonietti et al., Chin, Lasserre and Sukumar) would only speed
up the polynomial integrals of the cut cells, which take 0.12 s out of 1.05 s
in 2D (the right-hand side and the errors, with non-polynomial integrands,
need a volume quadrature anyway): it is not worth rewriting the local
operators in terms of polynomial moments. In 3D, the cut cells dominate the
assembly (about 3000 quadrature points per part for k = 3).

The method assumes that the mesh resolves the interface: each cell is cut
once, i.e. the interface crosses each face at most once. The detection of the
cut faces only checks the signs of the level set at the nodes; when the
refined interface of a cell leaves it (a face crossed twice, e.g. the flower
on level 1, where the interface dips 1.6 % of h below a face in 8 cells), the
cell is refined (`-x`, default). With `-x 0`, a warning is printed and the
points of the interface outside the cell are projected onto its faces, as in
ProtoN: the parts of the cell are then self-touching polygons whose
triangulations have negative weights, and the errors of the test cases with
jumps are unreliable in these cells (flower, (59), level 1, k = 1: negative
sum of the squared errors, NaN).

## Meshes and dimension

The cut geometry, the bases, the operators and the pairing are written for any
mesh of diskpp, in 2D and 3D (`-M`):

| `-M` | mesh | level l |
|------|------|---------|
| 0 | quadrangles (Cartesian grid, paper) | (10·2^l)² squares |
| 1 | triangles | (10·2^l)² squares split into 2 right triangles |
| 2 | polygons | hexagonal mesh of the level l+2 of the FVCA5 benchmark |
| 3 | hexahedra (Cartesian grid) | (5·2^l)³ cubes |
| 4 | tetrahedra | (5·2^l)³ cubes split into the 6 tetrahedra of the Kuhn subdivision |
| 5 | polyhedra | hexagonal prisms: level l+1 of FVCA5 extruded in 5·2^l layers |

**Interfaces** (`-i`), centered at (0.5, 0.5[, 0.5]); the new ones are
rotated (angle 0.4 in 2D, rotations of angles 0.5, 0.3, 0.7 about z, y, x in
3D) so as not to be aligned with the meshes:

| `-i` | 2D | 3D |
|------|----|----|
| 0 | circle of radius `-R` | sphere of radius `-R` |
| 1 | flower (55b): ρ² − R² + 0.03 cos 8θ | sphere with 4 petals: ρ² − R² + 0.03 Re(w⁴), w = (x̃ + i ỹ)/ρ |
| 2 | square (Fig. 12, quadrangles, level 0 only) | – |
| 3 | flower of ProtoN | – |
| 4 | ellipse, semi-axes 0.36, 0.2 | ellipsoid, semi-axes 0.36, 0.26, 0.18 |
| 5 | rounded square ∑ y_d⁴ = 0.3⁴ | rounded cube ∑ y_d⁴ = 0.3⁴ |
| 6 | annulus 0.18 < ρ < 0.42 (two circles) | torus R = 0.28, r = 0.14 |

The test cases (56)-(59) are the same functions in 3D (of ρ, and of
∏ sin(π x_d), e^x cos y); (57), (58a), (58b) need the circle or the sphere,
(56), (59) and the contrast test 5 hold for any interface (their jumps are
computed with the normal of the level set). On the meshes finer than the
level 0, the vertical sides of the square lie on mesh lines: the parts of
the cells along them have a zero area and the system is singular (the paper
only uses the square on the level 0).

**Structure.** `disk::cut_mesh<Mesh>` (`libdiskpp/.../mesh/cut_mesh.hpp`) is a
mesh of diskpp of any type with the description of its cut by the level set:
location of the nodes, faces and cells, geometry of the cut faces and cells,
cached quadratures, neighbors and pairing. The geometry itself
(`cut_geometry.hpp`) does not depend on the mesh:

* *polygons cut once* (cells of a 2D mesh, faces of a 3D mesh in the
  coordinates of their plane): exact interface as above, parts integrated by
  sectors from a point with respect to which they are star-shaped, recursive
  splitting by the line through the middle of the interface (thin slivers
  along an edge), polyline of 2^r segments and ear clipping as fallback;
* *convex polyhedra cut once* (cells of a 3D mesh): the interface curves of
  the cut faces form a loop; the interface is the graph, in the direction ν
  of the mean normal, of the level set above the region D of a plane
  enclosed by the projection of the loop. D is integrated by sectors from its
  center (its boundary being the projected curves), and the points are
  lifted onto the zero level set (weights |∇φ|/|∇φ·ν|, normals ∇φ/|∇φ|).
  Each part is integrated by cones from a point with respect to which it is
  star-shaped, the bases of the cones being its flat faces (whole faces, or
  parts of the cut faces) and the interface: positive weights, points in the
  part. A polyhedron for which D or a part is not star-shaped (thin strips of
  the interface in thin tetrahedra) is split by the plane through the middle
  of its interface orthogonal to its largest extent, recursively (3 levels);
  the last resort is a subdivision into tetrahedra with planar interfaces
  (never needed in the tests below). The splitting is a quadrature device:
  the cell keeps one part on each side. A face or a cell cut more than once
  (several loops) stops the computation: the mesh must resolve the
  interface.

The orthonormal bases are built in the principal axes of each region in any
dimension (with a canonical sign of the axes, so that the translates of a
cell get the same basis), the face bases in 3D are orthonormalized in the
plane of the face. The neighbors Δ(T) are the cells sharing a face, then a
vertex; on very coarse meshes, a cut diagonal neighbor that is not ill-cut on
the same side is accepted. An ill-cut side is flagged by its measure
(|T^i| < ϑ|T|; ProtoN's criterion on the cut faces is kept for the
quadrangles). In the agglomeration, groups targeting each other are merged,
and the interface of an agglomerated cell need not be connected (3D).

**Validation.** Volume of the ball, area of the sphere and flux of the normal
computed with the cut quadratures (degree 6): relative errors below 1e-13 on
all the 3D meshes from 10 cells per direction on (1e-9 on the coarsest
ones), all the quadrature points on the right side, positive weights in the
cut cells; no cell needs the subdivision (hexahedra: no split; tetrahedra:
19 %, 1.2 % and 1.6 % of the cut cells are split on the levels 0, 1, 2;
prisms: 10 % on the level 0, less than 1 % on the finer levels). Convergence of the
polynomial extension (default settings), test (56), orders of the L2 / H1
errors between the last two levels:

| mesh | k = 0 | k = 1 | k = 2 | k = 3 |
|------|-------|-------|-------|-------|
| triangles, l = 2→3 | 2.01 / 1.00 | 3.00 / 1.99 | 4.01 / 3.00 | 4.76 / 3.99 |
| hexagons, l = 2→3 | 1.96 / 1.06 | 2.96 / 1.98 | 4.25 / 3.04 | 5.29 / 4.03 |
| hexahedra, l = 1→2 | 2.05 / 1.12 | 3.15 / 2.09 | 4.15 / 3.12 | 4.97 / 4.02 |
| tetrahedra, l = 1→2 | 2.04 / 1.08 | 3.01 / 1.98 | 4.11 / 3.01 | 5.35 / 4.03 |
| prisms, l = 1→2 | 1.91 / 1.14 | 2.95 / 2.01 | 4.15 / 3.05 | 5.64 / 4.17 |

The contrast κ₂ = 100 κ₁ (57), the jumps (58a), (58b), (59) and the cell
agglomeration (`-p1`) converge at the same rates on these meshes (e.g.
tetrahedra, (59), k = 2: L2 errors 2.0e-4, 1.3e-5, 7.9e-7). Cost in 3D, k = 3,
hexahedra, level 2 (8000 cells, 824 cut cells): assembly 13 s, solve
(Pardiso, 250 000 face unknowns) 7.6 s, errors 5.7 s.

**Campaign on the structured meshes** (default settings, local refinement
`-x 3` triggered where needed): the 6 meshes, all the smooth interfaces
(`-i` 0, 1, 3, 4, 5, 6 in 2D, 0, 1, 4, 5, 6 in 3D), the tests (56), (59)
and 5 (contrast κ₂ = 100 κ₁) for every interface and (57), (58a), (58b) for
the circle and the sphere, both stabilizations (polynomial extension up to
k = 3, agglomeration up to k = 2), levels 0 to 3 (2D) or 2 (3D): 219 runs, of
which 109 used refined meshes. The least-squares orders over the levels ≥ 1
are within 0.25 of the optimal ones (k+2 for L2, k+1 for H1) in 217 runs. The
two others: the torus on the prisms with a contrast, k = 3 (4.45 / 3.41 from
level 1 to 2, 4.96 / 3.86 on average from level 0, pre-asymptotic on these
coarse prisms), and the square on the hexagons (not supported, see below).

**Limitations.** The cells must be convex for the splitting (other cells
fall back to the subdivision when the direct construction fails), and each
cell and face must be cut once: the local refinement (`-x`) ensures it for
smooth interfaces; the interfaces must be smooth (except the square on the
level 0) and must not touch the boundary of the domain. The background
meshes are meant to be structured (`-M` 0 to 5): unstructured meshes were
not pursued. In debug
builds, `measure()` of diskpp asserts on the hexahedra of `cartesian_mesh`
(`assert(pts.size() == 4)` instead of 8, not modified).

## Local refinement

The meshes do not resolve every interface (petals and torus on the coarsest
meshes, edges nearly tangent to the interface): some cells are not cut once.
`disk::cut_once_violations` (`cut_mesh.hpp`) finds the cells

* with an edge that the level set crosses more often than the signs at its
  ends tell (16 points per edge), in particular twice with no change of sign,
* with a face, or which are, cut more than once,
* where the interface leaves a face or the cell between two crossings of its
  boundary, or pokes through a face or the cell without crossing its edges,
* (3D) whose parts could only be integrated by subdivision.

A curved interface that cannot be built while the cell is cut once (corner
of the square) is not a violation: the polyline is then a simple polygon.

These cells are split as on the next level of the mesh, and so on, at most
`-x` times (`common/refinement.hpp`): triangles into 4 triangles, the other
polygons into quadrangles joining the middles of their edges to their
barycenter (quadrangles into 4, hexagons into 6), tetrahedra into 8 (the
inner octahedron split along its shortest diagonal), prisms and hexahedra by
splitting their base polygon and their height. The refinement is balanced:
two cells sharing a point differ by at most one level. Without this
balance, a large ill-cut cell could be paired with a cell refined 3 times,
whose polynomial, extended to the ill-cut cell, is extrapolated 8 times
beyond its cell (rounded cube, hexahedra, level 2, k = 3: H1 error with the
extension 4.2e-3 instead of 1.4e-5). The coarser neighbors of a refined cell
keep their faces, split at the new points, i.e. the refinement creates
hanging nodes. For HHO methods, a cell with a hanging
node is just a polygon (polyhedron) with aligned vertices, whose sub-faces
carry their own face unknowns: the refined mesh is a polytopal mesh
(`disk::generic_mesh`), with no conformity condition. Some geometric
computations of diskpp degenerate on aligned vertices (normal of a face from
its first three vertices, Kahan formula for a triangle of zero area,
convexity test of the polygons); for such elements of the generic meshes,
the overloads of `cut_mesh` compute the quadratures, measures, centroids and
normals from fans of triangles (tetrahedra) around the average of the
vertices (Newell normal), the other elements keeping the computations of
diskpp. The cells of the 3D generic meshes (refined meshes, prisms) are
integrated on tetrahedra with the rules of Xiao and Gimbutas (`arbq`, positive
weights, as on the tetrahedral meshes) instead of the Grundmann-Möller rules
of diskpp, whose weights are negative from the degree 2 on: they are exact
for the polynomials, but the computed squared errors of non-polynomial
functions could be negative (L2 error of the ellipsoid, tetrahedra, k = 3:
NaN, or underestimated 3 times).

By default (`-x 3`), the meshes of all the levels of the run are first
checked: if no cell violates the assumption, the run uses the meshes of
`-M` unchanged (same results as before); otherwise all the levels use the
locally refined meshes, run as polytopal meshes, and h is the largest
diameter of the cells (that of the background mesh). `-x 0` keeps the
meshes (projection of the interface onto the faces, see above).

## Choice of η_N

s^N = η_N κ_i h_T⁻² ‖v_S⁺ − v_T‖²_{T^i} was studied for k = 0..3: errors of the
solution (56) (circle and flower) and of the solution (57) with κ₂ = 10⁴κ₁
(circle) on the levels 0..4, condition numbers for the interfaces of Fig. 12
(circles and squares on the level 0, circles on the levels 1 and 2):

* the scaling h_T⁻² is the right one: the best η_N does not drift with h;
* η_N must grow with k, and η_N = 20·(k+1)² (20, 80, 180, 320) is the default.
  With the bases and pairing of the paper (`-b 1 -P 0`), it keeps the errors
  of Theorem 4.6 within 14 % (coarsest mesh) and 4 % (finer meshes) of those of
  the best η_N, and the condition number of the worst cut configuration within
  a factor 2 of the best one; with η_N = (k+1)², these factors reach 3.7 for
  the errors and 15 for the condition number on the coarsest mesh. With the
  default bases and pairing, the conditioning hardly depends on η_N up to
  40·(k+1)², and 20·(k+1)² keeps the errors of Theorem 4.6 within 13 % / 9 %
  of the best ones among c·(k+1)², c in [1, 200]: larger values bring the own
  unknowns of the ill-cut subcells closer to the extension of the paired
  polynomial, which is less accurate;
* the errors computed with the extension of the paired polynomial (ProtoN
  norms) reach their plateau only for η_N ≈ 50·(k+1)^{5/2}, at the expense of
  the conditioning: they are within a factor 2.5 of their best value with the
  default law, and up to 34 with η_N = (k+1)²;
* the peaks of the condition number in Fig. 12 (e.g. R ≈ 0.396) do not depend
  on η_N: they come from the scaled monomial bases on thin ill-cut subcells and
  from the pairing of these subcells with thin neighbors (see the default
  discretization above);
* with the scaled monomial bases, round-off errors dominate for k ≥ 4 from
  level 2 on.

## Differences with ProtoN

The local matrices and right-hand sides coincide with ProtoN's. The remaining
differences are:

* quadratures: diskpp uses rules exact for the requested degree (ProtoN's
  Dunavant table is shifted by one degree); this only affects the integration
  of the non-polynomial data, at the level of 1e-5 relative;
* on the coarsest mesh, the pairing of cells in symmetric configurations
  (exact ties) is resolved deterministically, whereas ProtoN resolves it by
  round-off;
* the orientation of the interface quadrature no longer depends on a test on
  the barycenter of the subcell, which failed for non-convex subcells (flower
  interface, levels 0 and 1);
* the quadrature on extremely thin subcells (square interface, small δ) uses
  signed cross-product areas instead of the Kahan formula, which returned NaN;
* with cell agglomeration, the Dirichlet jump of the right-hand side uses the
  same bases and penalty as the matrix;
* when the interface crosses a face twice (flower, level 1), the cells are
  refined (`-x`); with `-x 0`, the excursion out of the cell is projected onto
  the face (see above), whereas ProtoN integrated over the resulting
  self-intersecting polygons, which changes the errors of the flower on
  level 1 by at most 1e-3 relative for k <= 2 and 2 % for k = 3.
