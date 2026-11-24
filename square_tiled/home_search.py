import json,glob,itertools,sympy as sp,numpy as np
from fractions import Fraction
from ohom import Origami
t=sp.symbols('t')
target=sp.expand((t**2+14*t+1)*(t**2+4*t+1)); target2=sp.expand((t**2-14*t+1)*(t**2-4*t+1))
cats=[]
for f in sorted(glob.glob('cat_*.json')): cats+=json.load(open(f))
def tofr(M): return np.array([[Fraction(int(sp.fraction(x)[0]),int(sp.fraction(x)[1])) for x in row] for row in M.tolist()],dtype=object)
res=[];summary=[]
for idx,c in enumerate(cats):
    O=Origami(c['r'],c['u'],c['sigma'])
    if O.g2!=8: continue
    Bp,Bm=O.blocks(); pinv=(Bp.T*Bp).inv()*Bp.T
    G={('H',1):tofr(O.Th),('H',-1):tofr(O.Thi),('V',1):tofr(O.Tv),('V',-1):tofr(O.Tvi)}
    def D(c,k): return np.array([[1,k*O.sh],[0,1]],dtype=object) if c=='H' else np.array([[1,0],[-k*O.sv,1]],dtype=object)
    words=[]
    for L,ks in ((2,[-4,-3,-2,-1,1,2,3,4]),(4,[-4,-3,-2,-1,1,2,3,4]),(6,[-2,-1,1,2])):
        for ex in itertools.product(ks,repeat=L):
            words.append(tuple(('H' if i%2==0 else 'V',k) for i,k in enumerate(ex)))
    n52=0;found=[];plus_polys=set()
    for w in words:
        Dm=np.eye(2,dtype=object)
        for cc,k in w: Dm=Dm.dot(D(cc,k))
        if abs(Dm[0,0]+Dm[1,1])!=52: continue
        n52+=1
        M=np.eye(8,dtype=object)*Fraction(1)
        for cc,k in w:
            g=G[(cc,1 if k>0 else -1)]
            for _ in range(abs(k)): M=M.dot(g)
        Mp=pinv*sp.Matrix(M.tolist())*Bp
        p=sp.expand(Mp.charpoly(t).as_expr()); plus_polys.add(str(sp.factor(p)))
        if p==target or p==target2: found.append([list(x) for x in w])
    summary.append({'idx':idx,'N':O.N,'sh':O.sh,'sv':O.sv,'n_trace52':n52,'H1plus_charpolys_at_trace52':sorted(plus_polys),'matches':found})
    print(idx,O.N,O.sh,O.sv,n52,len(found),sorted(plus_polys)[:3],flush=True)
json.dump({'n_origamis':len(summary),'results':summary},open('home_search_out.json','w'),indent=1)
