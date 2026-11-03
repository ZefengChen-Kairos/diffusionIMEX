#!/usr/bin/env python3
"""One-shot convergence study of the travelling MHD vortex (paper 2, Sec. 4.1 / Table 1).

What it does, in one command:
  1. generates one source file per mesh from Final_emb_omp_patched.cpp (test = testVortex, nxCells = nyCells = cells),
  2. compiles them (in parallel),
  3. runs every (rho0, N) case (in parallel, OMP_NUM_THREADS per job), t_f = 1, CFL 0.2 (the value used for the
     table; --cfl changes it), gamma = 5/3,
     B~ = 1, v~ = sqrt(2 pi), u0 = (1, 1), p0 = 1,
     domain rule of Table 1: h = 10/N; rho0 >= 1e-1 -> [-7.5, 7.5]^2 with 1.5 N cells, rho0 <= 1e-2 -> [-5, 5]^2 with N cells,
  4. parses the "VORTEX N=... L2: rho=... u=... p=... Bx=... Az=..." line of every run (L2 errors of the cell averages
     against the exact translated vortex, computed inside the code) and
  5. prints / writes the Table-1-style EOC table (Markdown, CSV, LaTeX) and a results.json (used to resume).

Usage (typical, 32 cores):
    python3 run_vortex_eoc.py --jobs 4 --threads 8
Options:
    --src PATH          source file (default ../patch/Final_emb_omp_patched.cpp)
    --outdir DIR        work/output directory (default ./eoc_run)
    --rho0 LIST         comma list, default 1,0.1,0.01,0.001,0.0001,0.00001
    --N LIST            comma list of the nominal N (h = 10/N), default 32,64,128,256
    --domain paper|L5   'paper' = Table-1 rule above (default); 'L5' = [-5,5]^2 with N cells for every rho0
    --cfl X             default 0.2
    --orig              run the ORIGINAL scheme (plap=0 hllemflat=0 hllemcap=0) instead of fixes A + C
    --extra "k=v ..."   extra run-time options passed to every run (e.g. "pswitch=1")
    --max-N N           skip cases whose cell count exceeds N (smoke tests)
    --dry-run           show the plan only
    --no-build          reuse the binaries in outdir/bin
Re-running the same command skips cases already present in results.json (delete it to start over).
"""
import argparse, json, math, os, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = ['rho', 'u', 'p', 'Bx', 'Az']


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--src', default=os.path.join(HERE, '..', 'patch', 'Final_emb_omp_patched.cpp'))
    ap.add_argument('--outdir', default=os.path.join(os.getcwd(), 'eoc_run'))
    ap.add_argument('--rho0', default='1,0.1,0.01,0.001,0.0001,0.00001')
    ap.add_argument('--N', default='32,64,128,256')
    ap.add_argument('--domain', choices=['paper', 'L5'], default='paper')
    ap.add_argument('--cfl', type=float, default=0.2)
    ap.add_argument('--jobs', type=int, default=2, help='concurrent runs')
    ap.add_argument('--threads', type=int, default=4, help='OMP threads per run')
    ap.add_argument('--build-jobs', type=int, default=4)
    ap.add_argument('--orig', action='store_true')
    ap.add_argument('--extra', default='')
    ap.add_argument('--max-N', type=int, default=0)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--no-build', action='store_true')
    ap.add_argument('--cxx', default='g++')
    ap.add_argument('--cxxflags', default='-O3 -march=native -fopenmp -std=c++17 -I/usr/include/eigen3')
    return ap.parse_args()


def case_setup(rho0, N, rule):
    """cells per direction and half-width L for a nominal N (h = 10/N)."""
    if rule == 'paper' and rho0 >= 0.1:
        return int(round(1.5 * N)), 7.5
    return N, 5.0


def make_source(src, cells, dst):
    s = open(src).read()
    n = 0
    for old, new in [('Test test = testblastwave;', 'Test test = testVortex;'),
                     ('const int nxCells = 1000;', 'const int nxCells = %d;' % cells),
                     ('const int nyCells = 1000;', 'const int nyCells = %d;' % cells)]:
        if old not in s:
            # tolerate a source whose defaults differ: match the declaration pattern instead
            pat = {'Test test = testblastwave;': r'^Test test = \w+;', 'const int nxCells = 1000;': r'^const int nxCells = \d+;',
                   'const int nyCells = 1000;': r'^const int nyCells = \d+;'}[old]
            s, k = re.subn(pat, new, s, count=1, flags=re.M)
        else:
            s = s.replace(old, new, 1); k = 1
        n += k
    if n != 3:
        sys.exit('could not set test/mesh in %s' % src)
    open(dst, 'w').write(s)


def build(cxx, flags, srcfile, exe):
    cmd = [cxx] + flags.split() + [srcfile, '-o', exe]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-3000:]); sys.exit('compilation failed: %s' % srcfile)
    return exe


def run_case(exe, workdir, rho0, L, cfl, threads, extra, orig):
    os.makedirs(workdir, exist_ok=True)
    opts = ['cfl=%g' % cfl, 'rho0=%.17g' % rho0, 'vL=%g' % L, 'novtk=1', 'tag=eoc']
    if orig:
        opts += ['plap=0', 'hllemflat=0', 'hllemcap=0']
    if extra:
        opts += extra.split()
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    t0 = time.time()
    with open(os.path.join(workdir, 'run.log'), 'w') as log:
        r = subprocess.run([exe] + opts, cwd=workdir, stdout=log, stderr=subprocess.STDOUT, env=env)
    txt = open(os.path.join(workdir, 'run.log')).read()
    m = re.search(r'VORTEX N=(\d+) rho0=(\S+) t=(\S+) L2: rho=(\S+) u=(\S+) p=(\S+) Bx=(\S+) Az=(\S+)', txt)
    steps = len(re.findall(r'^step ', txt, flags=re.M))
    if r.returncode != 0 or not m:
        return {'ok': False, 'steps': steps, 'wall': time.time() - t0, 'log': os.path.join(workdir, 'run.log')}
    return {'ok': True, 'steps': steps, 'wall': time.time() - t0, 't': float(m.group(3)),
            'L2': dict(zip(FIELDS, [float(m.group(4 + k)) for k in range(5)]))}


def eoc_table(results, rho0s, Ns):
    """results[(rho0, N)] -> dict with 'L2'; returns (markdown, csv, latex)."""
    md = ['| rho0 | N | ' + ' | '.join('L2(%s) | EOC' % f for f in FIELDS) + ' |',
          '|---|---|' + '---|---|' * len(FIELDS)]
    csv = ['rho0,N,' + ','.join('L2_%s,EOC_%s' % (f, f) for f in FIELDS)]
    tex = ['\\begin{tabular}{rr' + 'cc' * len(FIELDS) + '}', '\\hline',
           '$\\rho_0$ & $N$ & ' + ' & '.join('$L^2(%s)$ & EOC' % {'rho': '\\rho', 'u': 'u', 'p': 'p', 'Bx': 'B_x', 'Az': 'A_z'}[f] for f in FIELDS) + ' \\\\', '\\hline']
    for rho0 in rho0s:
        prev = None
        for N in Ns:
            r = results.get((rho0, N))
            if not r or not r.get('ok'):
                continue
            e = r['L2']; cells = [('%g' % rho0) if prev is None else '', str(N)]
            md_c, csv_c, tex_c = [], [], []
            for f in FIELDS:
                eoc = math.log(prev[f] / e[f]) / math.log(2) if prev else float('nan')
                md_c += ['%.2e' % e[f], '-' if prev is None else '%.2f' % eoc]
                csv_c += ['%.6e' % e[f], '' if prev is None else '%.3f' % eoc]
                tex_c += ['%.2e' % e[f], '--' if prev is None else '%.2f' % eoc]
            md.append('| ' + ' | '.join(cells + md_c) + ' |')
            csv.append(','.join(['%g' % rho0, str(N)] + csv_c))
            tex.append(' & '.join(['$10^{%d}$' % round(math.log10(rho0)) if prev is None else '', str(N)] + tex_c) + ' \\\\')
            prev = e
        tex.append('\\hline')
    tex.append('\\end{tabular}')
    return '\n'.join(md), '\n'.join(csv), '\n'.join(tex)


def main():
    a = parse_args()
    rho0s = [float(x) for x in a.rho0.split(',')]
    Ns = [int(x) for x in a.N.split(',')]
    src = os.path.abspath(a.src); out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True); os.makedirs(os.path.join(out, 'bin'), exist_ok=True); os.makedirs(os.path.join(out, 'src'), exist_ok=True)
    resfile = os.path.join(out, 'results.json')
    results = {}
    if os.path.exists(resfile):
        for k, v in json.load(open(resfile)).items():
            r0, n = k.split('|'); results[(float(r0), int(n))] = v
    # plan
    cases = []
    for rho0 in rho0s:
        for N in Ns:
            cells, L = case_setup(rho0, N, a.domain)
            if a.max_N and cells > a.max_N:
                continue
            cases.append((rho0, N, cells, L))
    meshes = sorted(set(c[2] for c in cases))
    print('source  : %s' % src)
    print('scheme  : %s' % ('ORIGINAL (plap=0 hllemflat=0 hllemcap=0)' if a.orig else 'fixes A + C (patched defaults)') + (' + ' + a.extra if a.extra else ''))
    print('cases   : %d  (rho0 x N), meshes to build: %s' % (len(cases), meshes))
    for rho0, N, cells, L in cases:
        st = 'done' if results.get((rho0, N), {}).get('ok') else 'todo'
        print('   rho0=%-8g N=%-4d cells=%-4d L=%-4g  %s' % (rho0, N, cells, L, st))
    if a.dry_run:
        return
    # build
    exes = {}
    for cells in meshes:
        exes[cells] = os.path.join(out, 'bin', 'vortex_%d' % cells)
    todo_build = [c for c in meshes if not (a.no_build and os.path.exists(exes[c]))]
    if todo_build:
        print('building %d binaries with %d parallel compilations ...' % (len(todo_build), a.build_jobs)); t0 = time.time()
        with ThreadPoolExecutor(max_workers=a.build_jobs) as ex:
            futs = []
            for cells in todo_build:
                sf = os.path.join(out, 'src', 'vortex_%d.cpp' % cells); make_source(src, cells, sf)
                futs.append(ex.submit(build, a.cxx, a.cxxflags, sf, exes[cells]))
            for f in as_completed(futs):
                f.result()
        print('   built in %.0f s' % (time.time() - t0))
    # run
    todo = [c for c in cases if not results.get((c[0], c[1]), {}).get('ok')]
    todo.sort(key=lambda c: -c[2])   # largest meshes first
    print('running %d cases, %d at a time, %d threads each ...' % (len(todo), a.jobs, a.threads))

    def one(c):
        rho0, N, cells, L = c
        wd = os.path.join(out, 'run_rho%g_N%d' % (rho0, N))
        return c, run_case(exes[cells], wd, rho0, L, a.cfl, a.threads, a.extra, a.orig)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in as_completed([ex.submit(one, c) for c in todo]):
            (rho0, N, cells, L), r = f.result()
            results[(rho0, N)] = r
            json.dump({'%.17g|%d' % (k[0], k[1]): v for k, v in results.items()}, open(resfile, 'w'), indent=1)
            if r['ok']:
                print('   rho0=%-8g N=%-4d cells=%-4d steps=%-5d %6.0f s  L2: ' % (rho0, N, cells, r['steps'], r['wall']) + ' '.join('%s=%.3e' % (f, r['L2'][f]) for f in FIELDS))
            else:
                print('   rho0=%-8g N=%-4d FAILED after %d steps (%.0f s), see %s' % (rho0, N, r['steps'], r['wall'], r['log']))
    md, csv, tex = eoc_table(results, rho0s, Ns)
    open(os.path.join(out, 'eoc_table.md'), 'w').write(md + '\n'); open(os.path.join(out, 'eoc_table.csv'), 'w').write(csv + '\n'); open(os.path.join(out, 'eoc_table.tex'), 'w').write(tex + '\n')
    print('\nL2 errors of the cell averages at t = 1 (CFL %g, %s):\n' % (a.cfl, 'original scheme' if a.orig else 'fixes A + C'))
    print(md)
    print('\nwritten: %s/eoc_table.{md,csv,tex}, results.json' % out)


if __name__ == '__main__':
    main()
