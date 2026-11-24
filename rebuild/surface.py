"""D=5 Prym eigenform in H(6) built by the Thurston-Veech construction on a 10-rectangle grid."""
import sympy as sp
from origami import *
s5=sp.sqrt(5); phi=(1+s5)/2
r=[7,9,1,8,6,4,5,0,3,2]; u=[1,5,6,4,3,0,8,9,2,7]; sig=[8,6,5,7,9,2,1,3,0,4]
N=len(r)
Hc=cycles(r); Vc=cycles(u)
nh,nv=len(Hc),len(Vc)
E=sp.Matrix(nh,nv,lambda i,j: len(set(Hc[i])&set(Vc[j])))
def cyl_of(cyls,s):
    for i,c in enumerate(cyls):
        if s in c: return i
# Prym involution permutes cylinders (with reversal)
piH=[cyl_of(Hc,sig[c[0]]) for c in Hc]; piV=[cyl_of(Vc,sig[c[0]]) for c in Vc]
# PF data: E w = mu h, E^T h = mu w
M=E*E.T
ev=M.eigenvals()
mu2=max(ev,key=lambda e: float(e)); mu2=sp.nsimplify(sp.simplify(mu2),[s5])
mu=sp.sqrt(mu2); mu=sp.nsimplify(sp.simplify(mu),[s5])
hvec=(M-mu2*sp.eye(nh)).nullspace()[0]; hvec=sp.simplify(hvec/hvec[0])
wvec=sp.simplify(E.T*hvec/mu)
# homology basis: 8 cores, order [h_0..h_3, v_0..v_3]
Omega=sp.BlockMatrix([[sp.zeros(nh),E],[-E.T,sp.zeros(nv)]]).as_explicit()
Th=sp.BlockMatrix([[sp.eye(nh),E],[sp.zeros(nv,nh),sp.eye(nv)]]).as_explicit()
Tv=sp.BlockMatrix([[sp.eye(nh),sp.zeros(nh,nv)],[-E.T,sp.eye(nv)]]).as_explicit()
# NB columns = images of basis vectors: T_h(c^v_j)=c^v_j+sum_i E_ij c^h_i -> column j of block (0,1) is E[:,j]
P=sp.zeros(nh+nv)
for i in range(nh): P[piH[i],i]=-1
for j in range(nv): P[nh+piV[j],nh+j]=-1
# periods (Re, Im) of basis
per=sp.Matrix([[mu*hvec[i] for i in range(nh)]+[0]*nv,[0]*nh+[mu*wvec[j] for j in range(nv)]])
