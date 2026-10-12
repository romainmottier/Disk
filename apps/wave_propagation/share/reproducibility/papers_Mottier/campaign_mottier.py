#!/usr/bin/env python3
"""
campaign_mottier.py — machinerie commune des campagnes de reproduction des cas tests de deux articles sur le couplage
elasto-acoustique par methodes HHO, refaits avec l'application wave_propagation actuelle (prototype d'ordre un pilote
par fichier Lua) :

  [M2AN25]   R. Mottier, A. Ern, R. Khot, L. Guillot, Hybrid high-order methods for elasto-acoustic wave propagation
             in the time domain, ESAIM Math. Model. Numer. Anal. 59 (2025) 2685-2715 (dossier elastoacoustic_m2an2025)
  [CMAME25]  R. Mottier, A. Ern, L. Guillot, Elasto-acoustic wave propagation in geophysical media using hybrid
             high-order methods on general meshes, Comput. Methods Appl. Mech. Engrg. (2025) (dossier
             geophysics_cmame2025)

Un dossier par article : son script run_campaign.py (cas, references Gar6more, figures), un fichier d'entree Lua par
figure ou tableau, valeurs_article.json (valeurs relevees sur les figures publiees, valeurs exactes des tableaux) et
meshes/ (maillages generes par ../meshes.py, comprimes). Les sources ne contiennent que cela ; CMake recopie le dossier
dans le build, <build>/apps/wave_propagation/reproducibility/papers_Mottier, ou la campagne est lancee et ou elle ecrit
ses resultats, dans le dossier de l'article :
  figure_<n>/, tableau_<n>/    un dossier par figure ou tableau de l'article : la figure reproduite (PNG, panneaux
                               disposes comme dans l'article), la comparaison avec les valeurs relevees sur la figure
                               publiee, les variantes, et un dossier par cas :
    <cas>/<cas>.lua            entree complete du cas (fichier de la figure suivi des parametres du cas), executee
                               telle quelle : le cas se relance seul depuis son dossier
    <cas>/log.txt, <cas>.txt, <cas>.dat, <cas>.json   sortie, tableau, donnees, resultats relus par plot_mottier.py
    <cas>/<cas>_energy.dat, <cas>_<capteur>.csv, <cas>_cells.dat, <cas>_snapshot_<i>.dat   historiques, capteurs,
                               instantanes quand le cas les demande
  references/<ref>/            solutions de reference de Gar6more2D (../gar6more.py)
  meshes/                      maillages decompresses
  rapport.md                   pour chaque figure / tableau : parametres, ecarts, comparaison avec l'article, statut
"""
import argparse, bz2, json, os, re, shutil, subprocess, sys, threading, time
from pathlib import Path

HERE = Path(__file__).resolve().parent         # papers_Mottier (sources ou copie du build)
sys.path.insert(0, str(HERE))

PAPERS = {
    "M2AN25": "R. Mottier, A. Ern, R. Khot, L. Guillot, *Hybrid high-order methods for elasto-acoustic wave propagation "
              "in the time domain*, ESAIM Math. Model. Numer. Anal. 59 (2025) 2685-2715",
    "CMAME25": "R. Mottier, A. Ern, L. Guillot, *Elasto-acoustic wave propagation in geophysical media using hybrid "
               "high-order methods on general meshes*, Comput. Methods Appl. Mech. Engrg. (2025)",
}


def item_dir(item):
    """Dossier d'une figure ou d'un tableau de l'article : "Fig. 6" -> figure_6, "Table 2" -> tableau_2"""
    m = re.match(r"(Fig\.|Table)\s*(\d+)", item)
    return ("figure_" if m.group(1) == "Fig." else "tableau_") + m.group(2)


class Case:
    """Un calcul : entree de la figure (lua, dans le dossier de l'article) suivie des instructions du cas.
    item, panel : figure ou tableau de l'article et panneau ; meshes : maillages du dossier meshes/ de l'article
    (decompresses dans <sortie>/meshes, mesh.files = {...}) ; cost : duree estimee en minutes (campagne complete) ;
    threads : fils OpenMP / MKL du cas (1 par defaut ; CPU des etudes d'efficacite mesures sur un fil) ;
    quick : instructions ajoutees en mode rapide (quick_meshes : maillages du mode rapide, a la place de meshes) ;
    meta : donnees libres des figures (serie, legende, style...)"""
    def __init__(self, name, paper, item, panel, lua, statements=(), meshes=None, cost=1.0, threads=1, quick=(),
                 meta=None, variant="", mem=2.0, quick_meshes=None):
        self.name, self.paper, self.item, self.panel, self.lua = name, paper, item, panel, lua
        self.statements, self.meshes, self.cost, self.threads = list(statements), meshes, cost, threads
        self.mem = mem                                  # memoire de pointe estimee (Go) : limite --mem de la campagne
        self.quick, self.meta, self.variant = list(quick), dict(meta or {}), variant
        self.quick_meshes = quick_meshes
        self.folder = item_dir(item)


def case_dir(case, out):
    return out / case.folder / case.name


def in_build():
    """Copie du build : CMake y ecrit paths.json"""
    return (HERE / "paths.json").exists()


def find_executable(exe):
    if exe:
        return str(Path(exe).resolve())
    for cand in [p / "wave_propagation" for p in list(HERE.parents)[:3]] + [shutil.which("wave_propagation")]:
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return str(Path(cand).resolve())
    return None


def lua_string(s):
    return "'" + str(s).replace("\\", "\\\\").replace("'", "\\'") + "'"


_MESH_LOCK = threading.Lock()

# maillages construits a partir d'un autre (au lieu d'etre stockes) : nom -> (maillage de base, etiquettes raffinees)
DERIVED_MESHES = {"bassin_non_conforme.txt": ("bassin_quadrangles.txt", (2, 3))}


def prepare_meshes(case, out, paper_dir, write=True, quick=False):
    """Maillages du cas (dossier meshes/ de l'article, comprimes .bz2) decompresses dans <sortie>/meshes ; instruction
    Lua mesh.files = {...}, chemins relatifs au dossier du cas"""
    files = []
    for name in (case.quick_meshes if quick and case.quick_meshes else case.meshes):
        dst = out / "meshes" / name
        with _MESH_LOCK:
            if write and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                src = paper_dir / "meshes" / name
                tmp = dst.with_suffix(dst.suffix + ".tmp")
                if name in DERIVED_MESHES and not Path(str(src) + ".bz2").exists() and not src.exists():
                    import meshes
                    base, tags = DERIVED_MESHES[name]
                    pts, cells, tg = meshes.read_polygons(str(paper_dir / "meshes" / base) + ".bz2")
                    pts, cells, tg = meshes.refine(pts, cells, tg, set(tags))
                    meshes.write_polygons(tmp, pts, cells, tg, f"{base} raffine une fois dans les etiquettes {tags} "
                                                                "(noeuds pendants), construit par campaign_mottier.py")
                elif Path(str(src) + ".bz2").exists():
                    with bz2.open(str(src) + ".bz2", "rb") as fi, open(tmp, "wb") as fo:
                        shutil.copyfileobj(fi, fo)
                else:
                    shutil.copy(src, tmp)
                os.replace(tmp, dst)
        files.append(lua_string(os.path.relpath(dst, case_dir(case, out))))
    return "mesh.files = {" + ", ".join(files) + "}"


def case_statements(case, args, paper_dir):
    stmts = list(case.statements)
    if case.meshes is not None:
        stmts.append(prepare_meshes(case, args.out, paper_dir, write=not args.dry_run, quick=args.quick))
    if args.quick:
        stmts += case.quick
    return stmts + [f"output.file = {lua_string(case.name + '.txt')}", f"output.table = {lua_string(case.name + '.dat')}"]


def case_input(case, stmts, paper_dir):
    return (f"-- {case.name} ([{case.paper}] {case.item} {case.panel}) : entree complete ecrite par run_campaign.py, le\n"
            f"-- fichier {paper_dir.name}/{case.lua} suivi des parametres propres au cas. Relancer depuis ce dossier :\n"
            f"--   <build>/apps/wave_propagation/wave_propagation {case.name}.lua\n\n"
            + (paper_dir / case.lua).read_text().rstrip() + "\n\n-- parametres propres au cas\n" + "\n".join(stmts) + "\n")


# ------------------------------------------------------------------------------------------------------------------
# lecture des sorties

def read_table(path):
    if not path.exists():
        return None
    cols, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            cols = line[1:].split()
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    if cols is None:
        return None
    return {c: [r[i] if i < len(r) else None for r in rows] for i, c in enumerate(cols)}


def read_energy(path):
    if not path.exists():
        return None
    cols = ["t", "E", "E_acoustic", "E_elastic"]
    data, var = {c: [] for c in cols}, None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            m = re.search(r"max \|E - E0\| / E0 = (\S+)", line)
            if m:
                var = float(m.group(1))
        elif line.strip():
            for c, v in zip(cols, line.split()):
                data[c].append(float(v))
    data = {c: v for c, v in data.items() if v}
    data["max_variation"] = var
    return data


def read_csv(path):
    lines = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
    names = lines[0].split(",")
    cols = {n: [] for n in names}
    for l in lines[1:]:
        for n, v in zip(names, l.split(",")):
            cols[n].append(float(v))
    return cols


def read_sensors(d, name):
    out = {}
    for f in sorted(d.glob(f"{name}_*.csv")):
        out[f.stem[len(name) + 1:]] = read_csv(f)
    return out or None


def write_results(case, args, paper_dir, exe, cmd, status, seconds, order, stmts=()):
    d = case_dir(case, args.out)
    snaps = sorted(str(p.name) for p in d.glob(f"{case.name}_snapshot_*.dat"))
    res = {"case": case.name, "order": order, "paper": case.paper, "paper_ref": PAPERS[case.paper], "item": case.item,
           "folder": case.folder, "panel": case.panel, "variant": case.variant, "meta": case.meta,
           "status": status, "seconds": seconds, "quick": args.quick, "exe": exe, "command": cmd,
           "threads": case.threads, "input": f"{case.name}.lua", "figure_input": f"{paper_dir.name}/{case.lua}",
           "statements": list(stmts), "data": read_table(d / f"{case.name}.dat"),
           "energy": read_energy(d / f"{case.name}_energy.dat"), "sensors": read_sensors(d, case.name),
           "snapshots": snaps or None}
    (d / f"{case.name}.json").write_text(json.dumps(res, indent=1))


def estimated_minutes(cases, jobs):
    """Duree de la campagne si les cas, les plus longs d'abord, occupent leurs fils (cost minutes chacun)"""
    load = [0.0] * max(1, jobs)
    for c in sorted(cases, key=lambda c: -c.cost):
        i = load.index(min(load))
        load[i] += c.cost
    return max(load)


def case_needs(case, args):
    """(fils, Go) pris par un cas pendant son calcul"""
    return min(case.threads, args.jobs), 1.0 if args.quick else case.mem


def run_case(case, args, paper_dir, exe, order):
    d = case_dir(case, args.out)
    d.mkdir(parents=True, exist_ok=True)
    threads, mem = case_needs(case, args)
    print(f"  [debut {time.strftime('%H:%M')}] {case.name}" + ("" if args.quick else f" (~{case.cost:.0f} min)")
          + (f", {threads} fils" if threads > 1 else ""), flush=True)
    cmd, stmts, t0 = [], [], time.time()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    try:
        stmts = case_statements(case, args, paper_dir)
        (d / f"{case.name}.lua").write_text(case_input(case, stmts, paper_dir))
        cmd = [exe, f"{case.name}.lua"]
        with open(d / "log.txt", "w") as log:
            log.write("# " + " ".join(cmd) + "\n"); log.flush()
            try:
                rc = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env,
                                    timeout=args.timeout).returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
        status = "ok" if rc == 0 else f"echec ({rc})"
    except Exception as e:
        with open(d / "log.txt", "a") as log:
            log.write(f"[ERREUR] {e!r}\n")
        status = f"echec ({type(e).__name__})"
    sec = time.time() - t0
    write_results(case, args, paper_dir, exe, cmd, status, sec, order, stmts)
    return case.name, status, sec


class Slots:
    """Fils (--jobs) et memoire (--mem, Go) disponibles : un cas a n fils et m Go en prend autant. Les cas partent dans
    l'ordre donne (les plus longs d'abord) ; quand le suivant ne tient pas dans ce qui reste, le premier qui tient part
    a sa place (pas de fils inactifs derriere un gros cas). mem_file : fichier contenant la memoire utilisable, relu
    toutes les 30 s (budget ajustable pendant la campagne, quand d'autres calculs liberent ou prennent de la memoire)"""
    def __init__(self, n, mem, mem_file=None):
        self.free, self.mem, self.total_mem, self.mem_file = n, mem, mem, mem_file
        self.cv = threading.Condition()

    def _reread(self):
        try:
            total = float(Path(self.mem_file).read_text().split()[0])
        except Exception:
            return
        if total != self.total_mem:
            print(f"  [{time.strftime('%H:%M')}] memoire utilisable : {self.total_mem:g} -> {total:g} Go", flush=True)
            self.mem += total - self.total_mem
            self.total_mem = total

    def run(self, todo, needs, work):
        """work(c) pour chaque cas c de todo, needs(c) = (fils, Go) ; renvoie les resultats a mesure (generateur)"""
        pending, done, running = list(todo), [], [0]

        def target(c, k, m):
            try:
                r = work(c)
            finally:
                with self.cv:
                    self.free += k; self.mem += m; running[0] -= 1
                    done.append(r if "r" in locals() else (c.name, "echec (interne)", 0.0))
                    self.cv.notify_all()

        while pending or running[0] or done:
            with self.cv:
                while True:
                    if self.mem_file:
                        self._reread()
                    start = next((c for c in pending if needs(c)[0] <= self.free and needs(c)[1] <= self.mem), None)
                    if start is not None or done or (not running[0] and not pending):
                        break
                    if not running[0] and not self.mem_file:    # rien ne tourne et rien ne tient : on force le premier
                        start = pending[0]; break
                    self.cv.wait(timeout=30 if self.mem_file else None)
                results, done[:] = list(done), []
                if start is not None:
                    pending.remove(start)
                    k, m = needs(start)
                    self.free -= k; self.mem -= m; running[0] += 1
                    threading.Thread(target=target, args=(start, k, m), daemon=True).start()
            yield from results


def _default_mem():
    try:
        with open("/proc/meminfo") as f:
            kb = int(f.readline().split()[1])
        return 0.8 * kb / 2 ** 20
    except Exception:
        return 16.0


def main(all_cases, paper_dir, doc, references=(), figures_module=None):
    """Lance les cas et les references d'un article (script run_campaign.py de son dossier), puis ses figures et son
    rapport (plot_mottier.py)"""
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default=None, help="executable wave_propagation (defaut: celui du build)")
    ap.add_argument("--out", default=None, help="dossier des resultats (defaut, dans le build : ce dossier, ou rapide/ "
                                                 "avec --quick)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="fils utilises en parallele (defaut: la moitie des coeurs) ; un cas en prend `threads`")
    ap.add_argument("--quick", action="store_true", help="maillages et durees reduits (verification rapide)")
    ap.add_argument("--only", default=None, help="regex sur le nom des cas (et des references)")
    ap.add_argument("--list", action="store_true", help="liste des cas et leur duree")
    ap.add_argument("--dry-run", action="store_true", help="entree de chaque cas, sans rien executer ni ecrire")
    ap.add_argument("--no-plot", action="store_true", help="sans figures ni rapport")
    ap.add_argument("--no-run", action="store_true", help="figures et rapport seulement (resultats existants)")
    ap.add_argument("--skip-done", action="store_true", help="ne relance pas les cas dont le resultat (<cas>.json) est ok "
                                                             "(reprise d'une campagne interrompue)")
    ap.add_argument("--gar6more", default=None, help="executable Gar6more2D.out (defaut: compile dans le dossier "
                                                    "gar6more2d de papers_Mottier, telecharge si besoin)")
    ap.add_argument("--timeout", type=float, default=48 * 3600, help="duree maximale d'un cas en secondes (48 h)")
    ap.add_argument("--mem", type=float, default=_default_mem(), help="memoire utilisable par les cas simultanes, en Go "
                                                                      "(defaut: 80 %% de la memoire de la machine) ; "
                                                                      "les cas qui demandent plus ne sont pas lances")
    ap.add_argument("--mem-file", default=None, help="fichier contenant la memoire utilisable (Go), relu toutes les 30 s "
                                                     "(cree avec --mem s'il n'existe pas) : budget ajustable en cours de "
                                                     "campagne ; les cas trop gros attendent au lieu d'etre ecartes")
    ap.add_argument("--format", default="png", help="formats des figures (png, pdf, svg)")
    args = ap.parse_args()

    cases = [c for c in all_cases if not args.only or re.search(args.only, c.name)]
    refs = [r for r in references if not args.only or re.search(args.only, r.name)]
    if args.list:
        for c in cases:
            print(f"  {c.name:42s} {c.folder:11s} {c.panel:22s} {c.cost:7.1f} min {c.threads:2d} fil(s) {c.mem:5.1f} Go  {c.lua} "
                  + " ".join(c.statements) + (f"  [variante : {c.variant}]" if c.variant else ""))
        for r in refs:
            print(f"  {r.name:42s} reference Gar6more2D")
        tot = sum(c.cost * c.threads for c in cases)
        print(f"{len(cases)} cas, {tot / 60.0:.1f} h-fil de calcul (campagne complete)" if cases else "0 cas")
        return 0
    if args.out:
        args.out = Path(args.out).resolve()
    elif in_build():
        args.out = paper_dir / "rapide" if args.quick else paper_dir
    else:
        print("[ERREUR] les resultats s'ecrivent dans le build, pas dans les sources : lancer la copie de ce script dans "
              f"<build>/apps/wave_propagation/reproducibility/papers_Mottier/{paper_dir.name}, ou donner --out")
        return 1
    exe = find_executable(args.exe)
    if exe is None and not args.no_run:
        print("[ERREUR] executable wave_propagation introuvable (--exe)"); return 1
    if args.dry_run:
        for c in cases:
            print(f"[dry-run] {c.name} : {exe} {c.name}.lua, dans {case_dir(c, args.out)}\n"
                  f"    {c.name}.lua = {paper_dir.name}/{c.lua} + " + "; ".join(case_statements(c, args, paper_dir)))
        for r in refs:
            print(f"[dry-run] reference {r.name} : Gar6more2D, {r.describe()}")
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    n_fail = 0
    if not args.no_run:
        if refs:
            import gar6more
            g6 = gar6more.find_or_build(args.gar6more, HERE)
            for r in refs:
                try:
                    r.compute(g6, args.out / "references", quick=args.quick)
                    print(f"  reference {r.name:36s} ok", flush=True)
                except Exception as e:
                    n_fail += 1
                    print(f"  reference {r.name:36s} echec : {e!r}", flush=True)
        order = {c.name: i for i, c in enumerate(all_cases)}
        todo = sorted(cases, key=lambda c: -c.cost * c.threads)
        if args.skip_done:
            def is_done(c):
                try:
                    return json.loads((case_dir(c, args.out) / f"{c.name}.json").read_text()).get("status") == "ok"
                except Exception:
                    return False
            done_cases = [c for c in todo if is_done(c)]
            todo = [c for c in todo if c not in done_cases]
            print(f"{len(done_cases)} cas deja calcules (--skip-done), non relances", flush=True)
        if args.mem_file:                               # budget ajustable : valeur du fichier (cree avec --mem)
            mf = Path(args.mem_file)
            if mf.exists():
                args.mem = float(mf.read_text().split()[0])
            else:
                mf.write_text(f"{args.mem:g}\n")
        too_big = [] if args.mem_file else [c for c in todo if case_needs(c, args)[1] > args.mem]
        if too_big:                                     # jamais lances : relancer avec plus de memoire (--mem)
            todo = [c for c in todo if c not in too_big]
            print(f"{len(too_big)} cas non lances, memoire estimee au-dessus de --mem {args.mem:g} Go : "
                  + ", ".join(f"{c.name} ({c.mem:g} Go)" for c in too_big), flush=True)
        print(f"{len(todo)} cas, {args.jobs} fils, resultats dans {args.out}" + (" (mode rapide)" if args.quick else ""),
              flush=True)
        if not args.quick:
            print(f"duree estimee ~{estimated_minutes(todo, args.jobs):.0f} min (les cas les plus longs d'abord)", flush=True)
        t0, n_ok, n_done = time.time(), 0, 0
        slots = Slots(args.jobs, args.mem, args.mem_file)
        for name, status, sec in slots.run(todo, lambda c: case_needs(c, args),
                                           lambda c: run_case(c, args, paper_dir, exe, order[c.name])):
            n_ok += status == "ok"; n_fail += status.startswith("echec"); n_done += 1
            print(f"  [{n_done}/{len(todo)}] {name:40s} {status:20s} {sec:9.1f} s", flush=True)
        print(f"\n{n_ok} cas reussis, {n_fail} echec(s), {(time.time() - t0) / 60.0:.1f} min.", flush=True)

    if args.no_plot:
        print("Figures et rapport : python3 " + str(HERE / "plot_mottier.py") + f" {paper_dir.name} {args.out}")
        return 1 if n_fail else 0
    print("\nGeneration des figures et du rapport ...", flush=True)
    import plot_mottier
    rc = plot_mottier.make_all(paper_dir, args.out, figures_module, fmts=tuple(args.format.split(",")))
    return 1 if (n_fail or rc) else 0
