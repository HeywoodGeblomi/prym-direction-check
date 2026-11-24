import random, itertools
from fractions import Fraction
def inv(p):
    q=[0]*len(p)
    for i,x in enumerate(p): q[x]=i
    return q
def cycles(p):
    seen=set();out=[]
    for i in range(len(p)):
        if i in seen: continue
        c=[];j=i
        while j not in seen: seen.add(j);c.append(j);j=p[j]
        out.append(c)
    return out
def transitive(r,u):
    N=len(r);seen={0};st=[0]
    while st:
        x=st.pop()
        for y in (r[x],u[x]):
            if y not in seen: seen.add(y);st.append(y)
    return len(seen)==N
class UF:
    def __init__(s,n): s.p=list(range(n))
    def f(s,x):
        while s.p[x]!=x: s.p[x]=s.p[s.p[x]]; x=s.p[x]
        return x
    def u(s,a,b): s.p[s.f(a)]=s.f(b)
BL,BR,TL,TR=0,1,2,3
def vertices(r,u):
    N=len(r);uf=UF(4*N)
    for s in range(N):
        uf.u(4*s+BR,4*r[s]+BL); uf.u(4*s+TR,4*r[s]+TL)
        uf.u(4*s+TL,4*u[s]+BL); uf.u(4*s+TR,4*u[s]+BR)
    cls={}
    for k in range(4*N): cls.setdefault(uf.f(k),[]).append(k)
    return uf,list(cls.values())
def stratum(r,u):
    uf,cl=vertices(r,u)
    return sorted([len(c)//4-1 for c in cl if len(c)>4],reverse=True)
def prym_sigma(r,u):
    N=len(r);ri=inv(r);ui=inv(u);res=[]
    for t in range(N):
        sig=[-1]*N;sig[0]=t;st=[0];ok=True
        while st and ok:
            x=st.pop()
            for a,b in ((r[x],ri[sig[x]]),(u[x],ui[sig[x]]),(ri[x],r[sig[x]]),(ui[x],u[sig[x]])):
                if sig[a]==-1: sig[a]=b;st.append(a)
                elif sig[a]!=b: ok=False;break
        if ok and all(sig[sig[i]]==i for i in range(N)) and sorted(sig)==list(range(N)):
            res.append(sig)
    return res
def fixed_points(r,u,sig):
    N=len(r);uf,cl=vertices(r,u)
    f=sum(1 for s in range(N) if sig[s]==s)
    f+=sum(1 for s in range(N) if sig[s]==u[s])
    f+=sum(1 for s in range(N) if sig[s]==r[s])
    img={BL:TR,BR:TL,TL:BR,TR:BL}
    for c in cl:
        k=c[0];s,cn=divmod(k,4)
        if uf.f(4*sig[s]+img[cn])==uf.f(k): f+=1
    return f
