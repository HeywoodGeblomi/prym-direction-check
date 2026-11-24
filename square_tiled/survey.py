import itertools, random, math, json
def mul(X,Y): return (X[0]*Y[0]+X[1]*Y[2], X[0]*Y[1]+X[1]*Y[3], X[2]*Y[0]+X[3]*Y[2], X[2]*Y[1]+X[3]*Y[3])
I=(1,0,0,1)
def hmat(k): return (1,-k,0,1)   # transvection on horizontal-type core in its block
def vmat(k): return (1,0,k,1)
def step(state,letter):
    P,Q,A=state;i,k=letter
    if i==0: P=mul(hmat(k),P)
    elif i==2: P=mul(vmat(k),P)
    elif i==1: Q=mul(hmat(k),Q)
    else: Q=mul(vmat(k),Q)
    A=mul((1,0,k,1) if i in (0,1) else (1,k,0,1),A)
    return (P,Q,A)
def Mplus(P,Q):  # assemble 4x4 in basis c0,c1,c2,c3
    return ((P[0],0,P[1],0),(0,Q[0],0,Q[1]),(P[2],0,P[3],0),(0,Q[2],0,Q[3]))
# sanity: recorded example
s=(I,I,I)
for l in [(1,3),(2,4),(2,4),(3,2),(0,2)]: s=step(s,l)
assert Mplus(s[0],s[1])==((-15,0,-2,0),(0,1,0,-3),(8,0,1,0),(0,2,0,-5)) and s[2]==(31,10,65,21), 'example not reproduced'
def ell(t): t=abs(t); return math.acosh(t/2)
def rots(w): return [w[j:]+w[:j] for j in range(len(w))]
def winv(w): return tuple((i,-k) for i,k in reversed(w))
SYMS=[{0:0,1:1,2:2,3:3},{0:1,1:0,2:3,3:2},{0:2,1:3,2:0,3:1},{0:3,1:2,2:1,3:0}]
def c_rot(w): return min(rots(w))
def c_rr(w): return min(rots(w)+rots(winv(w)))
def c_sym(w): return min(c_rr(tuple((m[i],k) for i,k in w)) for m in SYMS)
def analyze(words_iter):
    c=dict(n=0,hypA=0,hypM=0,two=0,band=0,exact=0);hits=[];mats=set()
    for w,(P,Q,A) in words_iter:
        c['n']+=1
        tA=A[0]+A[3]; tP=P[0]+P[3]; tQ=Q[0]+Q[3]
        if abs(tA)>2: c['hypA']+=1
        hp,hq=abs(tP)>2,abs(tQ)>2
        if hp or hq: c['hypM']+=1
        if hp and hq:
            c['two']+=1
            big,small=(tP,tQ) if abs(tP)>=abs(tQ) else (tQ,tP)
            r=ell(big)/ell(small)
            if 1.9<=r<=2.1: c['band']+=1
            if abs(big)==small*small-2:
                c['exact']+=1;hits.append((w,tA,Mplus(P,Q)));mats.add(Mplus(P,Q))
    c['n_matrices_21']=len(mats)
    c['up_to_rotation']=len({c_rot(w) for w,_,_ in hits})
    c['up_to_rotation_reversal']=len({c_rr(w) for w,_,_ in hits})
    c['up_to_symmetry']=len({c_sym(w) for w,_,_ in hits})
    c['exact_and_hypA']=sum(1 for w,tA,_ in hits if abs(tA)>2)
    return c,hits
def exhaustive(L,ks):
    letters=[(i,k) for i in range(4) for k in ks]
    def rec(prefix,state):
        if len(prefix)==L: yield tuple(prefix),state; return
        for l in letters: yield from rec(prefix+[l],step(state,l))
    yield from rec([],(I,I,I))
def lengths(Ls,ks):
    for L in Ls: yield from exhaustive(L,ks)
def rnd(n,L,ks,seed):
    random.seed(seed);letters=[(i,k) for i in range(4) for k in ks]
    for _ in range(n):
        w=tuple(random.choice(letters) for _ in range(L));s=(I,I,I)
        for l in w: s=step(s,l)
        yield w,s
confs={'C1':lambda:exhaustive(5,[1,2,3,4]),'C2':lambda:exhaustive(5,[2,3,4]),'C3':lambda:lengths(range(1,6),[1,2,3,4]),
       'C4':lambda:lengths(range(1,5),[1,2,3,4]),'C5':lambda:rnd(10000,5,[1,2,3,4],0),'C6':lambda:rnd(10000,5,[2,3,4],0),
       'C7':lambda:exhaustive(4,[-4,-3,-2,-1,1,2,3,4])}
out={}
for name,f in confs.items():
    c,hits=analyze(f())
    # consistency of merges
    from collections import defaultdict
    g=defaultdict(set)
    for w,tA,M in hits:
        P=(M[0][0],M[0][2],M[2][0],M[2][2]);Q=(M[1][1],M[1][3],M[3][1],M[3][3])
        g[c_rr(w)].add((tuple(sorted([abs(P[0]+P[3]),abs(Q[0]+Q[3])])),abs(tA)))
    c['merge_violations']=sum(1 for v in g.values() if len(v)>1)
    c['example_hits']=[[list(map(list,w)),tA] for w,tA,_ in hits[:5]]
    out[name]=c; print(name,{k:v for k,v in c.items() if k!='example_hits'},flush=True)
json.dump(out,open('survey_out.json','w'),indent=1)
