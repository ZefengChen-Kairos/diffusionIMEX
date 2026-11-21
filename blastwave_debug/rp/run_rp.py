#!/usr/bin/env python3
"""One-shot runner for the 1D MHD Riemann problems (paper 2, Sec. 4.2 / Table 2).

In one command it builds, runs and plots the Riemann tests of the code.  The mapping between
the code's test names and the paper's labels is

    test1 = RP1 (Brio-Wu)   test2 = RP2 (Ryu-Jones)   test3 = RP3   test4 = RP4
    test5 = RP5             test6 = RP6               test7 = RP7

Every run is 1D: nxCells = --cells, nyCells = 1, domain [-0.5, 0.5], transmissive boundaries,
the final time built into the code for each test (RP5 0.04, RP6 0.03, RP7 0.25).  There is no
constrained transport in 1D (B_x is constant), so the magnetic field is advanced in flux form.

What it does:
  1. generates one source file per test from Final_emb_omp_patched.cpp (test, nxCells, nyCells set),
  2. compiles them in parallel,
  3. runs each (test, scheme) pair, where a scheme is a preset of run-time options:
        fix   : the patched defaults, i.e. fixes A + C                         (default)
        orig  : the original scheme, plap=0 hllemflat=0 hllemcap=0
        eps0  : eps_M = eps_B = 0, the fixed fully-implicit splitting ("standard IMEX" baseline;
                note the paper closes that baseline with the Rusanov flux (50), while this run
                keeps the HLLEM flux at eps = 0, so it is a close but not identical comparison)
     several schemes may be given and are run for every test,
  4. merges the code's per-variable field files into one ready-to-plot table per (test, scheme),
        <outdir>/data/rp<k>_<scheme>.dat
     with the columns
        x  rho  vx  vy  vz  p  Bx  By  Bz  E  eps_M  eps_B
     sorted in x, one row per cell, a '#' header line, whitespace separated (gnuplot / numpy /
     pgfplots read it directly).  No figures are produced.

Usage:
    python3 run_rp.py                                   # RP5, RP6, RP7 at 2000 cells, fixes A + C
    python3 run_rp.py --schemes fix,orig                # compare against the unfixed scheme
    python3 run_rp.py --tests 1,2,3,4 --schemes fix,eps0
    python3 run_rp.py --cells 500                       # quick check
Options:
    --src PATH      source file; by default Final_emb_omp_patched.cpp is looked for next to this
                    script, in the current directory, and in ../patch/
    --outdir DIR    work/output directory (default ./rp_run)
    --tests LIST    comma list of test numbers, default 5,6,7
    --schemes LIST  comma list of fix,orig,eps0 (default fix)
    --cells N       cells in x (default 2000)
    --cfl X         default 0.8
    --jobs / --threads / --build-jobs
    --extra "k=v"   extra run-time options appended to every run
    --dry-run, --no-build
Re-running skips (test, scheme) pairs already finished (delete <outdir>/results.json to redo).
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMES = {
    'fix':  ([], 'fixes A + C (patched defaults)'),
    'orig': (['plap=0', 'hllemflat=0', 'hllemcap=0'], 'original scheme'),
    'eps0': (['epsMmax=0', 'epsBmax=0'], 'eps_M = eps_B = 0 (fixed implicit splitting)'),
}
VARS = [('rho', r'$\rho$'), ('p', r'$p$'), ('By', r'$B_y$')]


def find_source(name='Final_emb_omp_patched.cpp'):
    tried = [os.path.join(HERE, name), os.path.join(os.getcwd(), name), os.path.join(HERE, '..', 'patch', name)]
    for p in tried:
        if os.path.exists(p):
            return os.path.abspath(p)
    sys.exit('%s not found. Looked in:\n  ' % name + '\n  '.join(os.path.abspath(p) for p in tried) +
             '\nPut it in one of these, or pass --src /full/path/to/%s' % name)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--src', default='')
    ap.add_argument('--outdir', default=os.path.join(os.getcwd(), 'rp_run'))
    ap.add_argument('--tests', default='5,6,7')
    ap.add_argument('--schemes', default='fix')
    ap.add_argument('--cells', type=int, default=2000)
    ap.add_argument('--cfl', type=float, default=0.8)
    ap.add_argument('--jobs', type=int, default=2)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--build-jobs', type=int, default=4)
    ap.add_argument('--extra', default='')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--no-build', action='store_true')
    ap.add_argument('--cxx', default='g++')
    ap.add_argument('--cxxflags', default='-O3 -march=native -fopenmp -std=c++17 -I/usr/include/eigen3')
    return ap.parse_args()


def make_source(src, test, cells, dst):
    s = open(src).read()
    subs = [(r'^Test test = \w+;', 'Test test = test%d;' % test),
            (r'^const int nxCells = \d+;', 'const int nxCells = %d;' % cells),
            (r'^const int nyCells = \d+;', 'const int nyCells = 1;')]
    for pat, new in subs:
        s, n = re.subn(pat, new, s, count=1, flags=re.M)
        if n != 1:
            sys.exit('could not apply "%s" to %s' % (new, src))
    open(dst, 'w').write(s)


def build(cxx, flags, srcfile, exe):
    r = subprocess.run([cxx] + flags.split() + [srcfile, '-o', exe], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-3000:]); sys.exit('compilation failed: %s' % srcfile)
    return exe


def run_case(exe, workdir, scheme, cfl, threads, extra):
    os.makedirs(workdir, exist_ok=True)
    opts = ['cfl=%g' % cfl, 'tag=%s' % scheme] + SCHEMES[scheme][0] + (extra.split() if extra else [])
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    t0 = time.time()
    with open(os.path.join(workdir, 'run.log'), 'w') as log:
        r = subprocess.run([exe] + opts, cwd=workdir, stdout=log, stderr=subprocess.STDOUT, env=env)
    txt = open(os.path.join(workdir, 'run.log')).read()
    steps = len(re.findall(r'^step ', txt, flags=re.M))
    # the code writes its field files into plot/<testname>/
    plotdirs = []
    pd = os.path.join(workdir, 'plot')
    if os.path.isdir(pd):
        plotdirs = [os.path.join(pd, d) for d in sorted(os.listdir(pd))]
    got = plotdirs[0] if plotdirs and os.path.exists(os.path.join(plotdirs[0], 'RK_rho.dat')) else None
    ok = (r.returncode == 0 and got is not None)
    warn = len(re.findall(r'WARNING|NEGP', txt))
    return {'ok': ok, 'steps': steps, 'wall': time.time() - t0, 'dir': got,
            'warnings': warn, 'log': os.path.join(workdir, 'run.log')}


def _column(d, fname):
    """RK_*.dat holds  x  y  value ; in 1D there is one y per x."""
    import numpy as np
    p = os.path.join(d, fname)
    if not os.path.exists(p):
        return None, None
    a = np.loadtxt(p)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    o = np.argsort(a[:, 0])
    return a[o, 0], a[o, 2]


# column name -> the file the code writes it to.  NOTE RK_p_line.dat does NOT hold the pressure:
# the original code writes u_pri[3] = v_z into it (the file name is a leftover and was kept).
COLUMNS = [('rho', 'RK_rho.dat'), ('vx', 'RK_vx.dat'), ('vy', 'RK_vy.dat'), ('vz', 'RK_p_line.dat'),
           ('p', 'RK_p.dat'), ('Bx', 'RK_Bx.dat'), ('By', 'RK_By.dat'), ('Bz', 'RK_Bz.dat'),
           ('E', 'RK_E.dat'), ('eps_M', 'RK_eps_M.dat'), ('eps_B', 'RK_eps_B.dat')]


def merge(d, dst, test, scheme, cells, cfl):
    """Merge the per-variable files of one run into a single table; returns the path."""
    import numpy as np
    x = None; cols = []; names = []
    for name, fname in COLUMNS:
        xi, v = _column(d, fname)
        if v is None:
            continue
        if x is None:
            x = xi
        elif len(xi) != len(x) or not np.allclose(xi, x):
            continue
        cols.append(v); names.append(name)
    if x is None:
        return None
    data = np.column_stack([x] + cols)
    hdr = ('RP%d, scheme %s (%s), %d cells on [-0.5, 0.5], CFL %g\n' % (test, scheme, SCHEMES[scheme][1], cells, cfl)
           + 'v_z is taken from RK_p_line.dat, which the code mislabels (it holds u_pri[3], not p)\n'
           + ' '.join(['x'] + names))
    np.savetxt(dst, data, fmt='%.10e', header=hdr)
    return dst


def main():
    a = parse_args()
    tests = [int(x) for x in a.tests.split(',')]
    schemes = [x.strip() for x in a.schemes.split(',')]
    for s in schemes:
        if s not in SCHEMES:
            sys.exit('unknown scheme "%s"; choose from %s' % (s, ', '.join(SCHEMES)))
    src = os.path.abspath(a.src) if a.src else find_source()
    if not os.path.exists(src):
        sys.exit('source file not found: %s' % src)
    out = os.path.abspath(a.outdir)
    for d in (out, os.path.join(out, 'bin'), os.path.join(out, 'src')):
        os.makedirs(d, exist_ok=True)
    resfile = os.path.join(out, 'results.json')
    results = json.load(open(resfile)) if os.path.exists(resfile) else {}

    print('source  : %s' % src)
    print('tests   : %s   (test<k> = RP<k>)' % ', '.join('RP%d' % t for t in tests))
    print('schemes : %s' % ', '.join('%s (%s)' % (s, SCHEMES[s][1]) for s in schemes))
    print('mesh    : %d x 1 cells on [-0.5, 0.5], CFL %g' % (a.cells, a.cfl))
    cases = [(t, s) for t in tests for s in schemes]
    for t, s in cases:
        print('   RP%-2d %-5s %s' % (t, s, 'done' if results.get('%d|%s' % (t, s), {}).get('ok') else 'todo'))
    if a.dry_run:
        return

    exes = {t: os.path.join(out, 'bin', 'rp%d_%d' % (t, a.cells)) for t in tests}
    todo_build = [t for t in tests if not (a.no_build and os.path.exists(exes[t]))]
    if todo_build:
        print('building %d binaries ...' % len(todo_build)); t0 = time.time()
        with ThreadPoolExecutor(max_workers=a.build_jobs) as ex:
            futs = []
            for t in todo_build:
                sf = os.path.join(out, 'src', 'rp%d_%d.cpp' % (t, a.cells)); make_source(src, t, a.cells, sf)
                futs.append(ex.submit(build, a.cxx, a.cxxflags, sf, exes[t]))
            for f in as_completed(futs):
                f.result()
        print('   built in %.0f s' % (time.time() - t0))

    todo = [c for c in cases if not results.get('%d|%s' % c, {}).get('ok')]
    print('running %d cases, %d at a time, %d threads each ...' % (len(todo), a.jobs, a.threads))

    def one(c):
        t, s = c
        return c, run_case(exes[t], os.path.join(out, 'rp%d_%s' % (t, s)), s, a.cfl, a.threads, a.extra)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for f in as_completed([ex.submit(one, c) for c in todo]):
            (t, s), r = f.result()
            results['%d|%s' % (t, s)] = r
            json.dump(results, open(resfile, 'w'), indent=1)
            if r['ok']:
                print('   RP%-2d %-5s %5d steps  %6.1f s%s' % (t, s, r['steps'], r['wall'],
                      '   (%d warnings in the log)' % r['warnings'] if r['warnings'] else ''))
            else:
                print('   RP%-2d %-5s FAILED after %d steps, see %s' % (t, s, r['steps'], r['log']))

    data = os.path.join(out, 'data'); os.makedirs(data, exist_ok=True)
    made = []
    for t in tests:
        for sch in schemes:
            r = results.get('%d|%s' % (t, sch))
            if not r or not r.get('ok'):
                continue
            dst = os.path.join(data, 'rp%d_%s.dat' % (t, sch))
            if merge(r['dir'], dst, t, sch, a.cells, a.cfl):
                made.append(dst)
    print('\nplot data (columns: x rho vx vy vz p Bx By Bz E eps_M eps_B):')
    for m in made:
        print('   %s' % m)
    print('\nraw per-variable files stay in <outdir>/rp<k>_<scheme>/plot/<name>/RK_*.dat;')
    print('per-step log (incl. divB) in the same directory, steplog_*.csv')


if __name__ == '__main__':
    main()
