import random,sys,json
from origami import *
seed=int(sys.argv[2]);random.seed(seed);N=int(sys.argv[1]);its=int(sys.argv[3])
def rand_inv(N):
    p=list(range(N));idx=list(range(N));random.shuffle(idx)
    for i in range(random.randint(0,N//2)):
        a,b=idx[2*i],idx[2*i+1];p[a]=b;p[b]=a
    return p
def comp(a,b): return [a[b[i]] for i in range(len(a))]
seen=set();out=[]
for it in range(its):
    s=rand_inv(N);r=comp(s,rand_inv(N));u=comp(s,rand_inv(N))
    if stratum(r,u)!=[6] or not transitive(r,u) or fixed_points(r,u,s)!=2: continue
    key=(tuple(sorted(len(c) for c in cycles(r))),tuple(sorted(len(c) for c in cycles(u))),
         tuple(sorted(len(c) for c in cycles([r[u[x]] for x in range(N)]))),tuple(sorted(len(c) for c in cycles([r[inv(u)[x]] for x in range(N)]))))
    if key in seen: continue
    seen.add(key);out.append({'r':r,'u':u,'sigma':s})
json.dump(out,open(f'cat_{N}_{seed}.json','w'))
print(N,len(out))
