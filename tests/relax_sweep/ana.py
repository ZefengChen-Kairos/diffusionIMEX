import numpy as np, re, json, os, sys
R='/home/user/diffusionIMEX/run/expo/'; PC='/home/user/diffusionIMEX/run/pc/'
EX=json.load(open(os.path.join(os.path.dirname(__file__),'exact_vals.json')))
C={'p':2,'u':3,'a2':6,'Y2':7,'T1':12,'T2':13}
_cache={}
def L(path):
    if path not in _cache: _cache[path]=np.loadtxt(path)
    return _cache[path]
def f(n): return n if n.startswith('/') else R+n+'.dat'
def log(n):
    p=(n[:-4] if n.endswith('.dat') else f(n)[:-4])+'.log'; s=open(p).readline()
    g=lambda pat: (re.search(pat,s).group(1) if re.search(pat,s) else '?')
    return dict(steps=g(r'steps (\d+)'),rej=g(r'rejects (\d+)'),stages=g(r'stages (\d+)'),nit=g(r'newton it (\d+)'),nps=g(r'\(([\d.]+)/stage\)'),
                res=g(r'maxres (\S+)'),fails=g(r'closure fails (\d+)'),dM2=g(r'dM2/M (\S+)'),wall=g(r'wall (\S+)'))
def logs(n): d=log(n); return '%s/%s/%s(%s) fails %s res %s' % (d['steps'],d['rej'],d['nit'],d['nps'],d['fails'],d['res'])
def ctr(n):
    d=L(f(n)); i=np.argsort(np.abs(d[:,0]-0.5))[:2]; m=d[i].mean(0)
    return dict(p=m[2],u=m[3],a2=m[6],Y2=m[7],T1=m[12],T2=m[13],dTmax=np.max(np.abs(d[:,13]-d[:,12])))
def dif(a,b,rel=True):
    A=L(f(a)); B=L(f(b)); o={}
    for k,c in C.items():
        d=np.max(np.abs(A[:,c]-B[:,c])); o[k]=d/np.max(np.abs(B[:,c])) if rel else d
    o['T']=max(o['T1'],o['T2']); return o
def g(x): return '%.7g'%x
def row(*a): print(' | '.join(str(x) for x in a))
mode=sys.argv[1] if len(sys.argv)>1 else 'all'
if mode in('E1','all'):
    print('## E1 centre values'); row('run','p','T1','T2','alpha2','Y2','max|T2-T1|','steps/rej/newton(per stage)')
    for n in ['E1_pTf_f0','E1_pTf_f1','E1_pTgf_f0','E1_pTgf_f1']:
        c=ctr(n); row(n,g(c['p']),'%.7f'%c['T1'],'%.7f'%c['T2'],g(c['a2']),g(c['Y2']),g(c['dTmax']),logs(n))
    print('## E1 BE vs EXP max|diff|/max|ref| (ref = BE)'); row('pair','p','u','alpha2','Y2','T')
    for a,b in [('E1_pTf_f1','E1_pTf_f0'),('E1_pTgf_f1','E1_pTgf_f0')]:
        o=dif(a,b); row(a+' vs f0',*[g(o[k]) for k in ('p','u','a2','Y2','T')])
    print('## E1 f0 vs run/pc bitwise?')
    for a,b in [('E1_pTf_f0','pel821_pTf'),('E1_pTgf_f0','pel821_pTgf')]:
        print(a, open(f(a),'rb').read()==open(PC+b+'.dat','rb').read())
if mode in('E2','all'):
    print('## E2 max|diff| (absolute) vs run/pc separate closures'); row('run','ref','p','u','alpha2','Y2','T','log')
    for k,ref in [('pT','pel821_pT'),('pTg','pel821_pTg')]:
        for fo in (1,0):
            for up in (1,0):
                n='E2_%s_f%d_up%d'%(k,fo,up); o=dif(n,PC+ref+'.dat',False); row(n,ref,*[g(o[q]) for q in ('p','u','a2','Y2','T')],logs(n))
        print('ref log',ref,logs(PC+ref+'.dat'))
    for c in (1,2):
        for r in (2,3):
            ref='ps%d_r%d'%(c,r)
            for up in (1,0):
                n='E2_%s_up%d'%(ref,up); o=dif(n,PC+ref+'.dat',False); row(n,ref,*[g(o[q]) for q in ('p','u','a2','Y2','T')],logs(n))
            print('ref log',ref,logs(PC+ref+'.dat'))
    print('## E2 rel differences max|diff|/max|ref|')
    for k,ref in [('pT','pel821_pT'),('pTg','pel821_pTg')]:
        for fo in (1,0):
            for up in (1,0):
                n='E2_%s_f%d_up%d'%(k,fo,up); o=dif(n,PC+ref+'.dat'); row(n,*[g(o[q]) for q in ('p','u','a2','Y2','T')])
    for c in (1,2):
        for r in (2,3):
            ref='ps%d_r%d'%(c,r)
            for up in (1,0):
                n='E2_%s_up%d'%(ref,up); o=dif(n,PC+ref+'.dat'); row(n,*[g(o[q]) for q in ('p','u','a2','Y2','T')])
    print('## E2 centre values (exact: pel pT p %.6g T %.7g ; pTg p %.6g T %.7g Y %.6g a2 %.6g)'%(EX['pel']['pT']['p'],EX['pel']['pT']['T'],EX['pel']['pTg']['p'],EX['pel']['pTg']['T'],EX['pel']['pTg']['Y'],EX['pel']['pTg']['a2']))
    row('run','p','T1','T2','alpha2','Y2','p-p_ex','T1-T_ex')
    for k,ref,ex in [('pT','pel821_pT','pT'),('pTg','pel821_pTg','pTg')]:
        E=EX['pel'][ex]
        for n in [PC+ref+'.dat']+['E2_%s_f%d_up%d'%(k,fo,up) for fo in (1,0) for up in (1,0)]:
            c=ctr(n); row(os.path.basename(n),g(c['p']),'%.7f'%c['T1'],'%.7f'%c['T2'],g(c['a2']),g(c['Y2']),g(c['p']-E['p']),g(c['T1']-E['T']))
    for c_ in (1,2):
        E=EX['ps%d'%c_]
        print('exact ps%d: pT %s ; pTg p %.6g T %.7g Y %.6g'%(c_, ('p %.6g T %.7g'%(E['pT']['p'],E['pT']['T'])) if E['pT']['centre'] else 'vacuum',E['pTg']['p'],E['pTg']['T'],E['pTg']['Y']))
        for r in (2,3):
            ref='ps%d_r%d'%(c_,r)
            for n in [PC+ref+'.dat','E2_%s_up1'%ref,'E2_%s_up0'%ref]:
                c=ctr(n); row(os.path.basename(n),g(c['p']),'%.7f'%c['T1'],'%.7f'%c['T2'],g(c['a2']),g(c['Y2']),g(c['p']-E['pTg']['p']),g(c['T1']-E['pTg']['T']))
def thname(th,fo,up=0):
    if th=='2000': return 'E1_pTf_f%d'%fo
    if th=='inf': return 'E2_pT_f%d_up%d'%(fo,up)
    return 'E3_th%s_f%d'%(th,fo)
def nuname(nu,fo,up=0):
    if nu=='1e-4': return 'E1_pTgf_f%d'%fo
    if nu=='inf': return 'E2_pTg_f%d_up%d'%(fo,up)
    return 'E3_nu%s_f%d'%(nu,fo)
THS=['0','20','200','2000','2e4','2e5','2e6','1e8','1e10','inf']; NUS=['1e-7','1e-6','1e-5','1e-4','1e-3','1e-2','1','10','100','1e4','inf']
if mode in('E3','all'):
    E=EX['pel']['pT']; c0=ctr(PC+'pel821_p.dat')
    print('## E3 theta sweep; exact 4-eq p %.7g T %.7g; exact 5-eq p %.7g; pel821_p numeric centre p %.7g T1 %.7f T2 %.7f'%(E['p'],E['T'],EX['pel']['p5']['p'],c0['p'],c0['T1'],c0['T2']))
    row('theta','form','p','T1','T2','p-p_ex4','T1-T_ex','T2-T_ex','log')
    for fo in (0,1):
        for th in THS:
            for up in ((0,1) if th=='inf' else (0,)):
                n=thname(th,fo,up); c=ctr(n); row(th+('(proj%d)'%up if th=='inf' else ''),'BE' if fo==0 else 'EXP',g(c['p']),'%.7f'%c['T1'],'%.7f'%c['T2'],g(c['p']-E['p']),g(c['T1']-E['T']),g(c['T2']-E['T']),logs(n))
    print('theta=0 vs pel821_p rel diff:',{k:g(v) for k,v in dif(thname('0',1),PC+'pel821_p.dat').items()}, 'BE:',{k:g(v) for k,v in dif(thname('0',0),PC+'pel821_p.dat').items()})
    E=EX['pel']['pTg']
    print('## E3 nu sweep; exact pTg p %.7g T %.7g Y %.6g'%(E['p'],E['T'],E['Y']))
    row('nu','form','p','T1','T2','Y2','p-p_ex','T1-T_ex','Y2-Y_ex','log')
    for fo in (0,1):
        for nu in NUS:
            for up in ((0,1) if nu=='inf' else (0,)):
                n=nuname(nu,fo,up); c=ctr(n); row(nu+('(proj%d)'%up if nu=='inf' else ''),'BE' if fo==0 else 'EXP',g(c['p']),'%.7f'%c['T1'],'%.7f'%c['T2'],g(c['Y2']),g(c['p']-E['p']),g(c['T1']-E['T']),g(c['Y2']-E['Y']),logs(n))
    print('## E3 BE-EXP domain rel diff per rate')
    for th in THS[1:-1]: o=dif(thname(th,1),thname(th,0)); row('theta',th,*[g(o[k]) for k in ('p','a2','T')])
    for nu in NUS[:-1]: o=dif(nuname(nu,1),nuname(nu,0)); row('nu',nu,*[g(o[k]) for k in ('p','a2','Y2','T')])
    print('## monotonicity of centre p vs rate (sign of successive differences), and |p - p(inf,proj0)| sequence')
    for fo in (0,1):
        P=[ctr(thname(t,fo))['p'] for t in THS]; print('theta',fo,np.sign(np.diff(P)).astype(int).tolist(),[g(abs(x-P[-1])) for x in P])
        T2=[ctr(thname(t,fo))['T2'] for t in THS]; print(' theta T2',fo,np.sign(np.diff(T2)).astype(int).tolist())
        T1=[ctr(thname(t,fo))['T1'] for t in THS]; print(' theta T1',fo,np.sign(np.diff(T1)).astype(int).tolist())
        P=[ctr(nuname(t,fo))['p'] for t in NUS]; print('nu',fo,np.sign(np.diff(P)).astype(int).tolist(),[g(abs(x-P[-1])) for x in P])
        Y=[ctr(nuname(t,fo))['Y2'] for t in NUS]; print(' nu Y',fo,np.sign(np.diff(Y)).astype(int).tolist())
        T=[ctr(nuname(t,fo))['T1'] for t in NUS]; print(' nu T1',fo,np.sign(np.diff(T)).astype(int).tolist())
if mode in('E4','all'):
    for k,base in (('pTf','E1_pTf_f%d'),('pTgf','E1_pTgf_f%d')):
        print('## E4',k); row('cfl','form','p','T1','T2','alpha2','Y2','log')
        nm=lambda c,fo: base%fo if c=='0.2' else 'E4_%s_c%s_f%d'%(k,c,fo)
        for c in ('0.2','0.1','0.05'):
            for fo in (0,1):
                n=nm(c,fo); x=ctr(n); row(c,'BE' if fo==0 else 'EXP',g(x['p']),'%.7f'%x['T1'],'%.7f'%x['T2'],g(x['a2']),g(x['Y2']),logs(n))
        row('cfl','|BE-EXP| centre p','T1','T2','a2','Y2','domain rel p','a2','Y2','T','| BE-ref(EXP .05) dom rel p','a2','Y2','T','| EXP-ref p','a2','Y2','T')
        ref=nm('0.05',1); cr=ctr(ref)
        for c in ('0.2','0.1','0.05'):
            a=ctr(nm(c,0)); b=ctr(nm(c,1)); o=dif(nm(c,1),nm(c,0)); r0=dif(nm(c,0),ref); r1=dif(nm(c,1),ref)
            row(c,*[g(abs(a[q]-b[q])) for q in ('p','T1','T2','a2','Y2')],*[g(o[q]) for q in ('p','a2','Y2','T')],*[g(r0[q]) for q in ('p','a2','Y2','T')],*[g(r1[q]) for q in ('p','a2','Y2','T')])
        print(' centre diff to ref (EXP 0.05): ')
        for c in ('0.2','0.1','0.05'):
            for fo in (0,1):
                a=ctr(nm(c,fo)); row(c,fo,*[g(a[q]-cr[q]) for q in ('p','T1','T2','a2','Y2')])
