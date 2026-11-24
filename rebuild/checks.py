from surface import *
import sympy as sp, json
out={}
out['squares']=N; out['r']=r; out['u']=u; out['sigma']=sig
out['stratum_zero_orders']=stratum(r,u)
out['prym_fixed_points']=fixed_points(r,u,sig)
out['horizontal_cylinders']=Hc; out['vertical_cylinders']=Vc
out['E']=E.tolist(); out['detE']=int(E.det())
out['mu2']=str(mu2); out['mu']=str(mu); out['h']=[str(x) for x in hvec]; out['w']=[str(x) for x in wvec]
out['piH']=piH; out['piV']=piV
I8=sp.eye(8)
out['Th_symplectic']=sp.simplify(Th.T*Omega*Th-Omega)==sp.zeros(8)
out['Tv_symplectic']=sp.simplify(Tv.T*Omega*Tv-Omega)==sp.zeros(8)
out['P_involution']=P*P==I8
out['P_symplectic']=P.T*Omega*P==Omega
out['P_commutes_Th']=P*Th==Th*P
out['P_commutes_Tv']=P*Tv==Tv*P
# derivative compatibility: per*T = D*per
Dh=sp.Matrix([[1,mu],[0,1]]); Dv=sp.Matrix([[1,0],[-mu,1]])
out['Th_matches_derivative']=sp.simplify(per*Th-Dh*per)==sp.zeros(2,8)
out['Tv_matches_derivative']=sp.simplify(per*Tv-Dv*per)==sp.zeros(2,8)
out['P_matches_minusI']=sp.simplify(per*P+per)==sp.zeros(2,8)
# normalized derivatives: conjugate by diag(mu,1)
C=sp.diag(mu,1)
out['Dh_normalized']=str(sp.simplify(C*Dh*C.inv())); out['Dv_normalized']=str(sp.simplify(C*Dv*C.inv()))
# cylinder geometry: circumference, height, modulus
out['horiz_cyl_circ_height']=[(str(sp.nsimplify(mu*hvec[i])),str(hvec[i])) for i in range(nh)]
out['vert_cyl_circ_width']=[(str(sp.nsimplify(mu*wvec[j])),str(wvec[j])) for j in range(nv)]
print(json.dumps(out,indent=1,default=str))
