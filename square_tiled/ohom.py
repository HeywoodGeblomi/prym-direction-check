"""Exact homology action of horizontal/vertical multitwists on a square-tiled surface (origami)."""
import sympy as sp
from math import lcm
from origami import *
class Origami:
    def __init__(s,r,u,sig):
        s.r,s.u,s.sig=r,u,sig;N=s.N=len(r)
        uf,cl=vertices(r,u);root={uf.f(c[0]):k for k,c in enumerate(cl)}
        vid=lambda sq,cn: root[uf.f(4*sq+cn)]
        nV=len(cl)
        d1=sp.zeros(nV,2*N)
        for x in range(N):
            d1[vid(x,1),x]+=1;d1[vid(x,0),x]-=1      # h_x: BL->BR
            d1[vid(x,2),N+x]+=1;d1[vid(x,0),N+x]-=1  # v_x: BL->TL
        Z=d1.nullspace()
        B=[]
        for x in range(N):
            v=sp.zeros(2*N,1);v[x]+=1;v[N+r[x]]+=1;v[u[x]]-=1;v[N+x]-=1;B.append(v)
        Bm=sp.Matrix.hstack(*B);rb=Bm.rank()
        basis=[];cur=Bm
        for z in Z:
            t=sp.Matrix.hstack(cur,z)
            if t.rank()>cur.rank(): basis.append(z);cur=t
        s.g2=len(basis);s.Hb=sp.Matrix.hstack(*basis);s.full=sp.Matrix.hstack(s.Hb,Bm)
        s.Rc=cycles(r);s.Uc=cycles(u)
        s.sh=lcm(*[len(c) for c in s.Rc]);s.sv=lcm(*[len(c) for c in s.Uc])
        s.Th=s._act(s._twist_h);s.Tv=s._act(s._twist_v);s.P=s._act(s._tau)
        s.Thi=s.Th.inv();s.Tvi=s.Tv.inv()
    def coords(s,ch):
        sol=s.full.solve_least_squares(ch) if False else s.full.gauss_jordan_solve(ch)[0]
        sol=sol.subs({x:0 for x in sol.free_symbols})
        return sol[:s.g2,0]
    def _act(s,f): return sp.Matrix.hstack(*[s.coords(f(s.Hb[:,j])) for j in range(s.g2)])
    def _twist_h(s,ch):
        N=s.N;out=ch.copy()
        for R in s.Rc:
            I=sum(ch[N+y] for y in R); n=s.sh//len(R)
            for x in R: out[x]+=n*I
        return out
    def _twist_v(s,ch):
        N=s.N;out=ch.copy()
        for U in s.Uc:
            J=sum(ch[y] for y in U); n=s.sv//len(U)
            for x in U: out[N+x]-=n*J
        return out
    def _tau(s,ch):
        N=s.N;out=sp.zeros(2*N,1)
        for x in range(N):
            out[s.u[s.sig[x]]]-=ch[x]; out[N+s.r[s.sig[x]]]-=ch[N+x]
        return out
    def hol(s):
        N=s.N; return sp.Matrix([[1]*N+[0]*N,[0]*N+[1]*N])*s.Hb
    def blocks(s):
        Bp=(s.P-sp.eye(s.g2)).nullspace();Bm=(s.P+sp.eye(s.g2)).nullspace()
        return sp.Matrix.hstack(*Bp),sp.Matrix.hstack(*Bm)
