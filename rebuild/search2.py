import random,sys
from origami import *
import sympy as sp
seed=int(sys.argv[2]);random.seed(seed)
N=int(sys.argv[1]);found={}
x=sp.symbols('x')
def rand_inv(N):
    p=list(range(N));idx=list(range(N));random.shuffle(idx)
    k=random.randint(0,N//2)
    for i in range(k):
        a,b=idx[2*i],idx[2*i+1];p[a]=b;p[b]=a
    return p
def comp(a,b): return [a[b[i]] for i in range(len(a))]
for it in range(int(sys.argv[3])):
    s=rand_inv(N);r=comp(s,rand_inv(N));u=comp(s,rand_inv(N))
    if stratum(r,u)!=[6] or not transitive(r,u): continue
    if fixed_points(r,u,s)!=2: continue
    H=cycles(r);V=cycles(u)
    E=sp.zeros(len(H),len(V))
    for i,h in enumerate(H):
        for j,v in enumerate(V): E[i,j]=len(set(h)&set(v))
    cp=sp.factor_list((E*E.T).charpoly(x).as_expr())
    key=(len(H),len(V),str(cp))
    if key not in found:
        found[key]=(r,u);print(key,r,u,s,flush=True)
