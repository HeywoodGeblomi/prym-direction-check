from surface import *
import sympy as sp
ri=inv(r); ui=inv(u)
hgt=[None]*N; wid=[None]*N
for i,c in enumerate(Hc):
    for s in c: hgt[s]=hvec[i]
for j,c in enumerate(Vc):
    for s in c: wid[s]=wvec[j]
def find_iso(r2,u2,w2,h2):
    """bijections f with f∘r2 = r∘f, f∘u2 = u∘f and rectangle sizes preserved (w2,h2 are sizes on source)"""
    res=[]
    for t in range(N):
        f=[-1]*N;f[0]=t;st=[0];ok=True
        while st and ok:
            x=st.pop()
            for a,b in ((r2[x],r[f[x]]),(u2[x],u[f[x]])):
                if f[a]==-1: f[a]=b;st.append(a)
                elif f[a]!=b: ok=False;break
        if ok and sorted(f)==list(range(N)) and all(sp.simplify(w2[s]-wid[f[s]])==0 and sp.simplify(h2[s]-hgt[f[s]])==0 for s in range(N)):
            res.append(f)
    return res
# candidate derivative G, transformed origami (r2,u2), sizes (w2,h2)
cands={
 'identity (translation automorphisms), G=I':(r,u,wid,hgt),
 'rotation by pi (Prym involution), G=-I':(ri,ui,wid,hgt),
 'rotation by +90deg, G=[[0,-1],[1,0]]':(ui,r,hgt,wid),
 'rotation by -90deg, G=[[0,1],[-1,0]]':(u,ri,hgt,wid),
 'reflection x->-x, G=diag(-1,1)':(ri,u,wid,hgt),
 'reflection y->-y, G=diag(1,-1)':(r,ui,wid,hgt),
 'reflection swap x,y, G=[[0,1],[1,0]]':(u,r,hgt,wid),
 'reflection x,y->-y,-x, G=[[0,-1],[-1,0]]':(ui,ri,hgt,wid),
}
found={k:len(find_iso(*v)) for k,v in cands.items()}
if __name__=='__main__':
    for k,v in found.items(): print(v,k)
    # maximal cylinders: does each r-cycle's top boundary contain the singular vertex?
    uf,cl=vertices(r,u); sing=[c for c in cl if len(c)>4][0]; singroot=uf.f(sing[0])
    for i,c in enumerate(Hc):
        top=any(uf.f(4*s+2)==singroot for s in c); bot=any(uf.f(4*s+0)==singroot for s in c)
        print('H cyl',i,'singular on top',top,'bottom',bot)
    for j,c in enumerate(Vc):
        L=any(uf.f(4*s+0)==singroot for s in c); R=any(uf.f(4*s+1)==singroot for s in c)
        print('V cyl',j,'singular on left',L,'right',R)
