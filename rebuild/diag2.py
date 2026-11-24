import itertools, numpy as np, collections, json, sympy as sp
from surface import Th,Tv,Omega
TH=np.array(Th.tolist(),dtype=float);TV=np.array(Tv.tolist(),dtype=float)
THi=np.linalg.inv(TH);TVi=np.linalg.inv(TV)
Om=np.array(Omega.tolist(),dtype=float)
Bp=np.array([[-1,0,1,0,0,0,0,0],[0,-1,0,1,0,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,0,-1,1]],float).T
Bm=np.array([[1,0,1,0,0,0,0,0],[0,1,0,1,0,0,0,0],[0,0,0,0,1,1,0,0],[0,0,0,0,0,0,1,1]],float).T
phi=(1+5**.5)/2; mu=phi**2
DH=lambda k:np.array([[1,k*mu],[0,1]]); DV=lambda k:np.array([[1,0],[-k*mu,1]])
ks=[-4,-3,-2,-1,1,2,3,4]; ratios=[]; maxerr=0; sympl=0; worst=None
for L in range(1,5):
  for st in 'HV':
    for ex in itertools.product(ks,repeat=L):
      M=np.eye(8); D=np.eye(2)
      for i,k in enumerate(ex):
        c='H' if (i%2==0)==(st=='H') else 'V'
        M=M@np.linalg.matrix_power(TH if c=='H' else TV, k) if k>0 else M@np.linalg.matrix_power(THi if c=='H' else TVi,-k)
        D=D@(DH(k) if c=='H' else DV(k))
      tr=abs(np.trace(D))
      if tr<=2: continue
      rho=(tr+np.sqrt(tr*tr-4))/2
      Mm=np.linalg.lstsq(Bm,M@Bm,rcond=None)[0]; Mp=np.linalg.lstsq(Bp,M@Bp,rcond=None)[0]
      em=max(abs(np.linalg.eigvals(Mm))); err=abs(np.log(em)/np.log(rho)-1)
      if err>maxerr: maxerr=err
      if np.abs(M.T@Om@M-Om).max()>1e-6: sympl+=1
      m=sorted(abs(np.linalg.eigvals(Mp)),reverse=True)
      if np.log(m[1])>1e-9: ratios.append(np.log(m[0])/np.log(m[1]))
r=np.array(ratios)
print('H1- top exponent equals log rho(A): max rel err',maxerr)
print('8x8 symplectic violations',sympl)
print('words with two positive H1+ exponents',len(r))
print('ratio l1/l2 quantiles',np.quantile(r,[0,.05,.25,.5,.75,.95,1]).round(4).tolist())
print('closest to 2:',sorted(r,key=lambda x:abs(x-2))[:5])
print('count in [1.9,2.1]',int(((r>=1.9)&(r<=2.1)).sum()))
