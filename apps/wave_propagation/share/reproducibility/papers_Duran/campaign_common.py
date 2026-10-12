#!/usr/bin/env python3
"""
campaign_common.py — machinerie commune des campagnes de verification des resultats de convergence de deux
articles anterieurs obtenus avec diskpp, refaites avec l'application wave_propagation actuelle :

  [CAMC22] E. Burman, O. Duran, A. Ern, Hybrid high-order methods for the acoustic wave equation in the
           time domain, Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, section 5.1
  [JSC21]  E. Burman, O. Duran, A. Ern, M. Steins, Convergence analysis of hybrid high-order methods for
           the wave equation, J. Sci. Comput. 87 (2021) 91, hal-02922720, section 6.1

Un dossier par article (acoustic_camc2022, elastic_jsc2021) : son script run_campaign.py (les cas : ce qui
varie d'une courbe a l'autre, ordres attendus) et un fichier d'entree Lua par figure (parametres de
l'article). Les sources ne contiennent que ces scripts et ces entrees ; CMake les recopie dans le dossier de
build, <build>/apps/wave_propagation/reproducibility/papers_Duran, ou la campagne est lancee et ou elle ecrit
ses resultats, dans le dossier de l'article :
  figure_<n>/                  un dossier par figure de l'article : la figure reproduite (<ART>_figure_<n>.png,
                               panneaux disposes comme dans l'article), chaque panneau seul, la comparaison avec
                               les valeurs relevees sur la figure publiee (..._comparaison.png), les variantes
                               (variantes/), les ordres observes (ordres_observes/), et un dossier par cas :
    <cas>/<cas>.lua            entree complete du cas (fichier de la figure suivi des parametres du cas),
                               executee telle quelle : le cas se relance seul depuis son dossier
    <cas>/log.txt, <cas>.txt, <cas>.dat, <cas>.json   sortie, tableau, donnees, resultats lus par
                               plot_convergence.py (et <cas>_energy.dat pour les historiques d'energie)
  rapport.md                   pour chaque figure : parametres, ecarts a la legende, rapports nous/article,
                               ordres attendus et observes, statut
  <article>_recapitulatif.png  toutes les figures de l'article
Les valeurs relevees sur les figures publiees sont dans valeurs_article.json (dossier de l'article).
"""
import argparse, bz2, concurrent.futures, json, os, re, shutil, subprocess, sys, threading, time
from pathlib import Path

HERE = Path(__file__).resolve().parent         # papers_Duran (sources ou copie du build)

PAPERS = {
    "CAMC22": "E. Burman, O. Duran, A. Ern, *Hybrid high-order methods for the acoustic wave equation in the time "
              "domain*, Commun. Appl. Math. Comput. 4 (2022) 597-633, hal-02922702, section 5.1",
    "JSC21":  "E. Burman, O. Duran, A. Ern, M. Steins, *Convergence analysis of hybrid high-order methods for the wave "
              "equation*, J. Sci. Comput. 87 (2021) 91, hal-02922720, section 6.1",
}

QUICK_SPACE = ["study.last_level = math.min(study.last_level, 4)"]
QUICK_TIME = ["study.time_levels = math.min(study.time_levels, 6)", "mesh.level = math.min(mesh.level, 5)"]


def col(name, norm, expected, window=None, tol=0.3, mode="min"):
    """Colonne d'erreur verifiee : ordre attendu, "min" (ordre optimal annonce : retenu >= attendu - tol) ou "equal"
    (ordre sous-optimal annonce : |retenu - attendu| <= tol), fenetre explicite de segments (python : (debut, fin)),
    sinon choix automatique de plot_convergence.py (plancher final exclu, deux derniers OCV)"""
    return {"name": name, "norm": norm, "expected": expected, "window": None if window is None else list(window),
            "tol": tol, "mode": mode}


def figure_number(panel):
    """Numero de la figure de l'article d'un panneau, "Fig. 3 (gauche)" -> 3"""
    return int(re.search(r"Fig\. (\d+)", panel).group(1))


class Case:
    # lua : fichier d'entree de la figure, dans le dossier de l'article ; cost : duree du cas dans la campagne
    # complete, en minutes sur un thread (mesuree sur une station de 24 coeurs) : les plus longs sont lances d'abord.
    # Figure reproduite : variant = "" pour la figure publiee, sinon nom de la variante (parametres de la legende quand
    # ils different de ceux de la figure publiee, ...) ; legend : titre de la legende du panneau (comme dans
    # l'article), labels : {colonne (ou "energy"): texte de la legende} ; style : {"color": 0..2 (bleu, rouge, orange),
    # "markers": marqueurs, "thick": trait epais clair (courbes d'energie tau = O(1/h))} ; ref : comparaison avec les
    # valeurs relevees sur la figure publiee, {"panel": cle de valeurs_article.json, "series": {colonne: serie}} ;
    # points : nombre de lignes attendues du tableau (verification : aucune valeur manquante ou non finie) ;
    # polymesher : maillages PolyMesher du dossier meshes/ de l'article (decompresses dans figure_<n>/meshes)
    def __init__(self, name, paper, panel, box, lua, statements, columns, series=("k", 0), x="h", quick=(),
                 fvca5=None, energy=False, energy_check=None, cost=1.0, note="", variant="", legend="", labels=None,
                 style=None, ref=None, points=None, polymesher=None):
        self.name, self.paper, self.panel, self.box, self.lua = name, paper, panel, box, lua
        self.statements, self.columns, self.series, self.x = list(statements), columns, series, x
        self.quick, self.fvca5, self.energy, self.energy_check = list(quick), fvca5, energy, energy_check
        self.cost, self.note = cost, note
        self.variant, self.legend, self.labels, self.style = variant, legend, dict(labels or {}), dict(style or {})
        self.ref, self.points, self.polymesher = ref, points, polymesher
        self.figure = figure_number(panel)


def case_dir(case, out):
    return out / f"figure_{case.figure}" / case.name


# ------------------------------------------------------------------------------------------------------------------

def in_build():
    """Copie du build : CMake y ecrit paths.json (emplacement des maillages de diskpp)"""
    return (HERE / "paths.json").exists()


def find_executable(exe):
    if exe:
        return str(Path(exe).resolve())
    # dans le build : <build>/apps/wave_propagation/reproducibility/papers_Duran/../../wave_propagation
    for cand in [p / "wave_propagation" for p in list(HERE.parents)[:3]] + [shutil.which("wave_propagation")]:
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return str(Path(cand).resolve())
    return None


def find_meshes(meshes):
    cands = [meshes] if meshes else []
    if in_build():
        cands.append(json.loads((HERE / "paths.json").read_text()).get("meshes"))
    cands += [p / "refactor_old_diskpp_code" / "meshes" for p in HERE.parents]     # depuis les sources de diskpp
    for c in cands:
        if c and (Path(c) / "2D_hex").is_dir():
            return Path(c).resolve()
    return None


def lua_string(s):
    return "'" + str(s).replace("\\", "\\\\").replace("'", "\\'") + "'"


_MESH_LOCK = threading.Lock()                   # les cas tournent dans des threads : une decompression par maillage


def prepare_fvca5(case, out, meshes, write=True):
    """Decompresse les maillages du cas dans <sortie>/figure_<n>/meshes (sauf si write = False), renvoie
    l'instruction Lua mesh.files = {...}, chemins relatifs au dossier du cas (ou le programme est lance)"""
    files = []
    for rel in case.fvca5:
        dst = out / f"figure_{case.figure}" / "meshes" / rel
        with _MESH_LOCK:
            if write and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                src = meshes / rel
                tmp = dst.with_suffix(dst.suffix + ".tmp")
                if Path(str(src) + ".bz2").exists():
                    with bz2.open(str(src) + ".bz2", "rb") as fi, open(tmp, "wb") as fo:
                        shutil.copyfileobj(fi, fo)
                else:
                    shutil.copy(src, tmp)
                os.replace(tmp, dst)
        files.append(lua_string(os.path.relpath(dst, case_dir(case, out))))
    return "mesh.files = {" + ", ".join(files) + "}"


def prepare_polymesher(case, out, paper_dir, write=True):
    """Maillages PolyMesher de l'article (dossier meshes/ de l'article, comprimes) decompresses dans
    <sortie>/figure_<n>/meshes (sauf si write = False) ; instruction Lua mesh.files = {...}, chemins relatifs au dossier
    du cas"""
    files = []
    for name in case.polymesher:
        dst = out / f"figure_{case.figure}" / "meshes" / name
        with _MESH_LOCK:
            if write and not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                src = paper_dir / "meshes" / name
                tmp = dst.with_suffix(dst.suffix + ".tmp")
                if Path(str(src) + ".bz2").exists():
                    with bz2.open(str(src) + ".bz2", "rb") as fi, open(tmp, "wb") as fo:
                        shutil.copyfileobj(fi, fo)
                else:
                    shutil.copy(src, tmp)
                os.replace(tmp, dst)
        files.append(lua_string(os.path.relpath(dst, case_dir(case, out))))
    return "mesh.files = {" + ", ".join(files) + "}"


def case_statements(case, args, meshes, paper_dir):
    """Parametres propres au cas, ajoutes au fichier de la figure : k, s, ..., maillages, mode rapide, sorties"""
    stmts = list(case.statements)
    if case.fvca5 is not None:
        stmts.append(prepare_fvca5(case, args.out, meshes, write=not args.dry_run))
    if case.polymesher is not None:
        stmts.append(prepare_polymesher(case, args.out, paper_dir, write=not args.dry_run))
    if args.quick:
        stmts += case.quick
    return stmts + [f"output.file = {lua_string(case.name + '.txt')}", f"output.table = {lua_string(case.name + '.dat')}"]


def case_input(case, stmts, paper_dir):
    """Entree Lua complete du cas : le fichier de la figure suivi des parametres propres au cas"""
    return (f"-- {case.name} ([{case.paper}] {case.panel}) : entree complete ecrite par run_campaign.py, le fichier\n"
            f"-- de la figure {paper_dir.name}/{case.lua} suivi des parametres propres au cas. Relancer depuis ce dossier :\n"
            f"--   <build>/apps/wave_propagation/wave_propagation {case.name}.lua\n\n"
            + (paper_dir / case.lua).read_text().rstrip() + "\n\n-- parametres propres au cas\n" + "\n".join(stmts) + "\n")


def read_table(path):
    if not path.exists():
        return None
    cols, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            cols = line[1:].split()
        elif line.strip():
            rows.append([float(v) for v in line.split()])
    if not cols:
        return None
    return {c: [r[i] if i < len(r) else None for r in rows] for i, c in enumerate(cols)}


def read_energy(path):
    if not path.exists():
        return None
    t, E, var = [], [], None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            m = re.search(r"max \|E - E0\| / E0 = (\S+)", line)
            if m:
                var = float(m.group(1))
        elif line.strip():
            a, b = line.split()[:2]
            t.append(float(a)); E.append(float(b))
    return {"t": t, "E": E, "max_variation": var}


def write_results(case, args, paper_dir, exe, cmd, status, seconds, order, stmts=()):
    d = case_dir(case, args.out)
    res = {"case": case.name, "order": order, "paper": case.paper, "paper_ref": PAPERS[case.paper], "panel": case.panel,
           "figure": case.figure, "box": case.box, "series": case.series[0], "series_value": case.series[1],
           "x": case.x, "columns": case.columns, "energy_check": case.energy_check, "note": case.note,
           "variant": case.variant, "legend": case.legend, "labels": case.labels, "style": case.style,
           "ref": case.ref, "points": None if args.quick else case.points, "paper_dir": paper_dir.name,
           "status": status, "seconds": seconds, "quick": args.quick, "exe": exe, "command": cmd,
           "input": f"{case.name}.lua", "figure_input": f"{paper_dir.name}/{case.lua}", "statements": list(stmts),
           "data": read_table(d / f"{case.name}.dat"),
           "energy": read_energy(d / f"{case.name}_energy.dat") if case.energy else None}
    (d / f"{case.name}.json").write_text(json.dumps(res, indent=1))


def estimated_minutes(cases, jobs):
    """Duree de la campagne complete si chaque cas, les plus longs d'abord, dure le temps mesure (cost)"""
    load = [0.0] * max(1, jobs)
    for c in sorted(cases, key=lambda c: -c.cost):
        load[load.index(min(load))] += c.cost
    return max(load)


def run_case(case, args, paper_dir, exe, meshes, order):
    d = case_dir(case, args.out)
    d.mkdir(parents=True, exist_ok=True)
    print(f"  [debut {time.strftime('%H:%M')}] {case.name}" + ("" if args.quick else f" (~{case.cost:.0f} min)"), flush=True)
    if case.fvca5 is not None and meshes is None:
        write_results(case, args, paper_dir, exe, [], "skipped", 0.0, order)
        return case.name, "ignore (maillages)", 0.0
    cmd, stmts, t0 = [], [], time.time()
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:                                        # une erreur dans un cas ne doit pas interrompre la campagne
        stmts = case_statements(case, args, meshes, paper_dir)
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


def main(all_cases, paper_dir, doc):
    """Lance les cas d'un article (script run_campaign.py de son dossier), puis ses figures et son rapport"""
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default=None, help="executable wave_propagation (defaut: celui du build, trouve dans les "
                                                 "dossiers parents, puis le PATH)")
    ap.add_argument("--out", default=None, help="dossier des resultats (defaut, dans le build : ce dossier, ou rapide/ "
                                                 "avec --quick)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help="cas lances en parallele, un thread chacun (defaut: la moitie des coeurs)")
    ap.add_argument("--quick", action="store_true", help="maillages <= niveau 4, 6 pas de temps (tolerance + 0.2)")
    ap.add_argument("--only", default=None, help="regex sur le nom des cas")
    ap.add_argument("--list", action="store_true", help="liste des cas et leur duree")
    ap.add_argument("--dry-run", action="store_true", help="entree de chaque cas, sans rien executer ni ecrire")
    ap.add_argument("--no-plot", action="store_true", help="sans figures ni rapport")
    ap.add_argument("--meshes", default=None, help="maillages de diskpp (refactor_old_diskpp_code/meshes)")
    ap.add_argument("--timeout", type=float, default=6 * 3600, help="duree maximale d'un cas en secondes (6 h)")
    ap.add_argument("--format", default="png", help="formats des figures, separes par des virgules : png, pdf, svg "
                                                    "(defaut: png)")
    args = ap.parse_args()

    cases = [c for c in all_cases if not args.only or re.search(args.only, c.name)]
    if args.list:
        for c in cases:
            print(f"  {c.name:38s} figure_{c.figure:<3d} {c.panel:24s} {c.cost:5.1f} min  {c.lua} {' '.join(c.statements)}"
                  + (f"  [variante : {c.variant}]" if c.variant else ""))
        tot = sum(c.cost for c in cases)
        print(f"{len(cases)} cas, " + (f"{tot / 60.0:.1f} h" if tot >= 60 else f"{tot:.1f} min")
              + f" de calcul sur un thread (campagne complete), cas le plus long {max(c.cost for c in cases):.1f} min"
              if cases else "0 cas")
        return 0
    if args.out:
        args.out = Path(args.out).resolve()
    elif in_build():
        args.out = paper_dir / "rapide" if args.quick else paper_dir
    else:
        print(f"[ERREUR] les resultats s'ecrivent dans le build, pas dans les sources : lancer la copie de ce script "
              f"dans <build>/apps/wave_propagation/reproducibility/papers_Duran/{paper_dir.name}, ou donner --out")
        return 1
    exe = find_executable(args.exe)
    if exe is None:
        print("[ERREUR] executable wave_propagation introuvable (--exe)"); return 1
    meshes = find_meshes(args.meshes)
    if meshes is None and any(c.fvca5 for c in cases):
        print("[WARN] maillages de diskpp introuvables (--meshes) : les cas sur maillages FVCA5 sont ignores")
    if args.dry_run:
        for c in cases:
            if c.fvca5 is not None and meshes is None:
                continue
            print(f"[dry-run] {c.name} : {exe} {c.name}.lua, dans {case_dir(c, args.out)}\n"
                  f"    {c.name}.lua = {paper_dir.name}/{c.lua} + " + "; ".join(case_statements(c, args, meshes, paper_dir)))
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    order = {c.name: i for i, c in enumerate(all_cases)}
    todo = sorted(cases, key=lambda c: -c.cost)
    print(f"{len(todo)} cas, {args.jobs} en parallele, resultats dans {args.out}"
          + (" (mode rapide)" if args.quick else ""), flush=True)
    if not args.quick:
        print(f"duree estimee ~{estimated_minutes(todo, args.jobs):.0f} min (les cas les plus longs d'abord) ; figures "
              f"dans les dossiers figure_<n> et rapport {args.out / 'rapport.md'} a la fin", flush=True)
    t0, n_ok, n_fail, n_done = time.time(), 0, 0, 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(run_case, c, args, paper_dir, exe, meshes, order[c.name]) for c in todo]
        for f in concurrent.futures.as_completed(futs):
            name, status, sec = f.result()
            n_ok += status == "ok"; n_fail += status.startswith("echec"); n_done += 1
            print(f"  [{n_done}/{len(todo)}] {name:34s} {status:20s} {sec:8.1f} s", flush=True)
    print(f"\n{n_ok} cas reussis, {n_fail} echec(s), {(time.time() - t0) / 60.0:.1f} min.", flush=True)

    plot = [sys.executable, str(HERE / "plot_convergence.py"), str(args.out), "--format", args.format]
    if args.no_plot:
        print("Figures et rapport : " + " ".join(plot))
        return 1 if n_fail else 0
    print("\nGeneration des figures et du rapport ...", flush=True)
    rc = subprocess.run(plot).returncode
    return 1 if (n_fail or rc) else 0
