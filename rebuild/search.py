import random,sys,collections
from origami import *
import sympy as sp
random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
N=int(sys.argv[1]);found={}
x=sp.symbols('x')
for it in range(int(sys.argv[3]) if len(sys.argv)>3 else 200000):
    r=list(range(N));random.shuffle(r);u=list(range(N));random.shuffle(u)
    if stratum(r,u)!=[6]: continue
    if not transitive(r,u): continue
    sigs=[s for s in prym_sigma(r,u) if fixed_points(r,u,s)==2]
    if not sigs: continue
    H=cycles(r);V=cycles(u)
    E=sp.zeros(len(H),len(V))
    for i,h in enumerate(H):
        for j,v in enumerate(V): E[i,j]=len(set(h)&set(v))
    cp=sp.factor_list((E*E.T).charpoly(x).as_expr())
    key=(len(H),len(V),str(cp))
    if key not in found:
        found[key]=(r,u);print(key,r,u,flush=True)
