#!/usr/bin/env python3
"""
gar6more.py — solutions de reference semi-analytiques de Gar6more2D (J. Diaz, A. Ezziani, INRIA, licence CeCILL,
https://gitlab.inria.fr/jdiaz/gar6more2d) pour les cas tests avec une ondelette de Ricker en condition initiale.

Gar6more2D (methode de Cagniard-de Hoop) calcule, dans deux demi-espaces (fluide y > 0, solide y < 0, interface y = 0),
la solution de  P_tt - c^2 Delta P = delta(x - x_s) f(t)  (et la vitesse V, dV/dt = -grad P / rho dans le fluide, la
vitesse du solide), source ponctuelle en x_s = (0, h), f(t) = derivee d'ordre `deriv` de exp(-pi^2 f0^2 (t - t_d)^2)
(au signe (-1)^(deriv+1) pres), nulle hors de (0, 2 t_d).

Equivalence avec les cas tests HHO (condition initiale, pas de source) : par decomposition en ondes planes, une
donnee initiale gaussienne w(0) = A exp(-a r^2), a = pi^2 / lambda^2, w_t(0) = 0 donne, pour t >> 1/f0, la meme solution
que la source ponctuelle de signature B d/dt exp(-pi^2 f0^2 t^2), B = c A lambda / sqrt(pi), f0 = c / lambda (erreur
d'ordre exp(-pi^2 f0^2 t^2)) ; on calcule donc la reponse retardee de t_d et on decale le temps de t_d.
  - flux initial (problem.ricker.field = "acoustic_flux") m0 = theta exp(-a r^2) (x - x_c), potentiel
    u0 = -rho theta / (2a) exp(-a r^2), p = du/dt : deriv = 2 (Ricker), p = S P_g6, S = c rho theta lambda^3 / (2 pi^5/2)
  - pression initiale ("acoustic_pressure") p0 = A exp(-a r^2) : deriv = 1, p = S P_g6, S = c A lambda / sqrt(pi)
  - vitesses (flux m du fluide, vitesse v du solide) : -S V_g6
Les figures de [M2AN25] sont tracees dans les unites de Gar6more2D (amplitude 1) : valeurs HHO divisees par S.

Compilation (find_or_build) : sources telechargees (git clone) dans papers_Mottier/gar6more2d/src, compilees dans
papers_Mottier/gar6more2d/build (cmake, gfortran, BLAS / LAPACK).
"""
import math, os, shutil, subprocess
from pathlib import Path

URL = "https://gitlab.inria.fr/jdiaz/gar6more2d.git"


def find_or_build(exe=None, here=None):
    if exe:
        return str(Path(exe).resolve())
    here = Path(here or Path(__file__).resolve().parent)
    root = here / "gar6more2d"
    built = root / "build" / "Gar6more2D.out"
    if built.exists():
        return str(built)
    found = shutil.which("Gar6more2D.out")
    if found:
        return found
    src = root / "src"
    if not (src / "CMakeLists.txt").exists():
        print(f"  Gar6more2D : telechargement de {URL} dans {src}", flush=True)
        root.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", "--depth", "1", URL, str(src)], check=True)
    b = root / "build"
    b.mkdir(parents=True, exist_ok=True)
    print(f"  Gar6more2D : compilation dans {b}", flush=True)
    with open(b / "build.log", "w") as log:
        subprocess.run(["cmake", str(src), "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
                        "-DCMAKE_Fortran_STANDARD_LIBRARIES=-llapack -lblas"], cwd=b, check=True, stdout=log,
                       stderr=subprocess.STDOUT)
        subprocess.run(["make", "-j", "8"], cwd=b, check=True, stdout=log, stderr=subprocess.STDOUT)
    if not built.exists():
        raise RuntimeError(f"compilation de Gar6more2D echouee (voir {b / 'build.log'})")
    return str(built)


def _read(path):
    rows = [[float(v) for v in l.split()] for l in path.read_text().splitlines() if l.strip()]
    return [r[0] for r in rows], [r[1] for r in rows]


class Reference:
    """Reference Gar6more2D d'un cas test : fluide (rho, c) au-dessus de y = 0, solide (rho, cp, cs) en dessous, ondelette
    centree en center (dans le fluide), frequence f0 (lambda = c_fluide / f0), champ initial field ("acoustic_flux" ou
    "acoustic_pressure") d'amplitude amplitude ; capteurs [(nom, x, y)] ; temps 0..T, pas dt ; nint : intervalles de la
    convolution numerique de Gar6more2D ; delay : retard du signal source, en periodes 1/f0. Le resultat ne devrait pas
    en dependre, mais la contribution de l'onde conique de Gar6more2D est faussee a partir de son arrivee moins delay :
    le delai fixe ou tombe ce defaut (1.5 reproduit les references publiees de [M2AN25])"""
    def __init__(self, name, fluid, solid, center, f0, field, amplitude, sensors, T, dt, nint=4000, delay=2.0):
        self.name, self.fluid, self.solid, self.center, self.f0 = name, fluid, solid, center, f0
        self.field, self.amplitude, self.sensors, self.T, self.dt, self.nint = field, amplitude, sensors, T, dt, nint
        self.delay = delay

    def scale(self):
        rho, c = self.fluid
        lam = c / self.f0
        if self.field == "acoustic_flux":
            return 2, c * rho * self.amplitude * lam ** 3 / (2.0 * math.pi ** 2.5)
        if self.field == "acoustic_pressure":
            return 1, c * self.amplitude * lam / math.sqrt(math.pi)
        raise ValueError(self.field)

    def describe(self):
        d, S = self.scale()
        return (f"fluide rho, c = {self.fluid}, solide rho, cp, cs = {self.solid}, source {self.center}, f0 = {self.f0}, "
                f"derivee {d}, S = {S:.6g}, retard {self.delay:g} / f0, capteurs {self.sensors}, T = {self.T}, dt = {self.dt}")

    def compute(self, exe, root, quick=False):
        out = Path(root) / self.name
        out.mkdir(parents=True, exist_ok=True)
        deriv, S = self.scale()
        rho_f, c_f = self.fluid
        rho_s, cp, cs = self.solid
        mu_s = rho_s * cs * cs
        lam_s = rho_s * cp * cp - 2.0 * mu_s
        td = self.delay / self.f0               # source negligeable en t = 0 (exp(-pi^2 delay^2)) et nulle apres 2 td
        xc, yc = self.center
        nint = min(self.nint, 1000) if quick else self.nint
        summary = {}
        for name, x, y in self.sensors:
            d = out / name
            d.mkdir(exist_ok=True)
            # capteur d'interface : evalue cote solide a 1e-3 L de l'interface (L = |x_c| + |y_c|) ; plus pres, l'integrale
            # de Gar6more2D devient instable (oscillations de l'ordre de la solution a 1e-9 L)
            Y = y if y != 0.0 else -1e-3 * max(1.0, abs(xc) + abs(yc))
            lines = ["4", "0", str(deriv), repr(self.f0), "1", repr(td), repr(yc), repr(Y), repr(x - xc), repr(x - xc), "1",
                     repr(td), repr(td + self.T), repr(self.dt), str(nint), f"{rho_f * c_f * c_f!r} {rho_f!r}", "1",
                     f"{mu_s!r} {lam_s!r} {rho_s!r}"]
            (d / "Gar6more2D.dat").write_text("\n".join(lines) + "\n")
            with open(d / "log.txt", "w") as log:
                subprocess.run([exe], cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
            t, P = _read(d / "P.dat")
            _, Ux = _read(d / "Ux.dat")
            _, Uy = _read(d / "Uy.dat")
            t = [ti - td for ti in t]
            fluid = y > 0.0
            cols = {"t": t}
            if fluid:
                cols.update({"p": [S * v for v in P], "mx": [-S * v for v in Ux], "my": [-S * v for v in Uy]})
                cols.update({"g6_p": P, "g6_vx": Ux, "g6_vy": Uy})
            else:
                cols.update({"vx": [-S * v for v in Ux], "vy": [-S * v for v in Uy], "g6_vx": Ux, "g6_vy": Uy})
            with open(out / f"{name}.csv", "w") as f:
                head = f"# Gar6more2D, capteur {name} ({x}, {y})" + (f" evalue en y = {Y:g}" if Y != y else "")
                f.write(head + f" ; valeurs HHO (p = S P, v = -S V, S = {S:.12g}) et brutes g6_* de Gar6more2D "
                        "(amplitude 1)\n" + ",".join(cols) + "\n")
                for i in range(len(t)):
                    f.write(",".join(f"{cols[c][i]:.12e}" for c in cols) + "\n")
            summary[name] = {"x": x, "y": y, "fluid": fluid}
        (out / "reference.txt").write_text(self.describe() + "\n")
        return summary


def load(root, name, sensor):
    """Colonnes d'une trace de reference calculee : {"t": [...], "p": [...], ...}"""
    path = Path(root) / name / f"{sensor}.csv"
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
    names = lines[0].split(",")
    cols = {n: [] for n in names}
    for l in lines[1:]:
        for n, v in zip(names, l.split(",")):
            cols[n].append(float(v))
    return cols
