import json, itertools, hashlib, sympy as sp, mpmath as mp
from surface import Th, Tv, E
mp.mp.dps=60
t=sp.symbols('t')
def tolist(M): return [[int(x) for x in row] for row in M.tolist()]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def mpow(A,k):
    n=len(A); R=[[int(i==j) for j in range(n)] for i in range(n)]
    if k<0: A=inv8[id(A)]; k=-k
    for _ in range(k): R=mm(R,A)
    return R
TH=tolist(Th); TV=tolist(Tv); THi=tolist(Th.inv()); TVi=tolist(Tv.inv())
inv8={id(TH):THi,id(TV):TVi}
# Z[phi] arithmetic: (a,b)=a+b*phi, phi^2=phi+1
def zm(x,y): a,b=x;c,d=y; return (a*c+b*d, a*d+b*c+b*d)
def za(x,y): return (x[0]+y[0],x[1]+y[1])
def dmm(A,B): return [[za(zm(A[i][0],B[0][j]),zm(A[i][1],B[1][j])) for j in range(2)] for i in range(2)]
phi4=(2,3)
def Dlet(c,k):
    if c=='H': return [[(1,0),(k*phi4[0],k*phi4[1])],[(0,0),(1,0)]]
    return [[(1,0),(0,0)],[(-k,0),(1,0)]]
def Hlet(c,k): return mpow(TH if c=='H' else TV,k)
Bp=[[-1,0,1,0,0,0,0,0],[0,-1,0,1,0,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,0,-1,1]]  # columns of H1+ basis
def restrict(M):
    # M*b_j = sum_i X_ij b_i ; basis vectors have disjoint supports with entries +-1, 2 nonzeros each
    X=[[0]*4 for _ in range(4)]
    for j,b in enumerate(Bp):
        img=[sum(M[r][c]*b[c] for c in range(8)) for r in range(8)]
        for i,bi in enumerate(Bp):
            idx=[c for c in range(8) if bi[c]!=0][1]  # coordinate where bi has +1
            X[i][j]=img[idx]*bi[idx]
        # verify
        rec=[sum(X[i][j]*Bp[i][c] for i in range(4)) for c in range(8)]
        assert rec==img, 'H1+ not invariant'
    return X
ks=[-4,-3,-2,-1,1,2,3,4]
words=[]
for L in range(1,5):
    for start in 'HV':
        types=[('H' if (i%2==0)==(start=='H') else 'V') for i in range(L)]
        for exps in itertools.product(ks,repeat=L):
            words.append(tuple(zip(types,exps)))
assert len(words)==9360
def phi_float(x): return x[0]+x[1]*(1+5**0.5)/2
def cyc_reduce(w):
    w=list(w)
    while len(w)>1 and w[0][0]==w[-1][0]:
        c=w[0][0]; k=w[0][1]+w[-1][1]; w=w[1:-1] if len(w)>2 else []
        if k!=0: w=[(c,k)]+w
        # re-merge interior if needed (alternating guaranteed otherwise)
    return tuple(w)
def rots(w): return [w[i:]+w[:i] for i in range(len(w))] if w else [w]
def winv(w): return tuple((c,-k) for c,k in reversed(w))
def canon_rot(w): return min(rots(cyc_reduce(w)))
def canon_rr(w): c=cyc_reduce(w); return min(rots(c)+rots(winv(c)))
def refl(w): return tuple((('V' if c=='H' else 'H'),-k) for c,k in w)
def canon_sym(w): return min(canon_rr(w),canon_rr(refl(w)))
rows=[]; nhyp=0; hits=[]; near=[]
for w in words:
    D=[[(1,0),(0,0)],[(0,0),(1,0)]]; M=[[int(i==j) for j in range(8)] for i in range(8)]
    for c,k in w: D=dmm(D,Dlet(c,k)); M=mm(M,Hlet(c,k))
    tr=za(D[0][0],D[1][1]); trf=phi_float(tr)
    if abs(trf)<=2: continue
    nhyp+=1
    Mp=restrict(M)
    p=sp.Matrix(Mp).charpoly(t).as_expr()
    coeffs=[int(c) for c in sp.Poly(p,t).all_coeffs()]
    roots=mp.polyroots(coeffs,maxsteps=200,extraprec=200)
    mods=sorted([abs(z) for z in roots],reverse=True)
    l1,l2=mp.log(mods[0]),mp.log(mods[1])
    if l2<=mp.mpf('1e-30'): continue
    ratio=l1/l2
    if abs(ratio-2)<mp.mpf('1e-40'):
        q=sp.Matrix(Mp)**2; q=q.charpoly(t).as_expr()
        exact=False
        for f,_ in sp.factor_list(p)[1]:
            fr=[abs(z) for z in mp.polyroots([int(c) for c in sp.Poly(f,t).all_coeffs()],maxsteps=200,extraprec=200)] if sp.degree(f,t)>0 else []
            if fr and abs(max(fr)-mods[0])<mp.mpf('1e-40'):
                if sp.degree(sp.gcd(f,q),t)>0 or sp.degree(sp.gcd(f,q.subs(t,-t)),t)>0: exact=True
        rho=(abs(trf)+(trf*trf-4)**0.5)/2
        full=sp.Matrix(M).charpoly(t).as_expr()
        rec={'word':[list(x) for x in w],'trace_Zphi':list(tr),'trace':trf,'rho':rho,'M_plus':Mp,'charpoly_M_plus':str(sp.factor(p)),
             'charpoly_8x8':str(sp.factor(full)),'raw':[1.0,float(l2/l1)],'per_log_rho':[float(l1)/float(mp.log(rho)),float(l2)/float(mp.log(rho))],
             'canon_rot':canon_rot(w),'canon_rr':canon_rr(w),'canon_sym':canon_sym(w)}
        (hits if exact else near).append(rec)
out={'n_words':len(words),'n_hyperbolic':nhyp,'n_directed_hits':len(hits),'n_numeric_only':len(near),
     'n_up_to_rotation':len({h['canon_rot'] for h in hits}),'n_up_to_rotation_reversal':len({h['canon_rr'] for h in hits}),
     'n_up_to_symmetry':len({h['canon_sym'] for h in hits}),
     'n_distinct_abs_trace':len({tuple(abs(x) for x in h['trace_Zphi']) if h['trace']>0 else tuple(-x for x in h['trace_Zphi']) for h in hits}),
     'n_distinct_charpoly8':len({h['charpoly_8x8'] for h in hits})}
# consistency: merged classes share |trace| and 8x8 charpoly
bad=0
from collections import defaultdict
g=defaultdict(set)
for h in hits: g[h['canon_rr']].add((round(abs(h['trace']),9),h['charpoly_8x8']))
bad=sum(1 for v in g.values() if len(v)>1); out['merge_consistency_violations']=bad
json.dump({'summary':out,'hits':hits,'numeric_only':near},open('survey_out.json','w'),indent=1,default=str)
print(json.dumps(out,indent=1))
