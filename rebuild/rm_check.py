from surface import *
import sympy as sp, itertools
# integral bases of H+ and H- (saturated) from P
def int_basis(Mat):
    ns=Mat.nullspace()
    B=sp.Matrix.hstack(*[v*sp.ilcm(*[sp.fraction(x)[1] for x in v]) for v in ns])
    # saturate: HNF of lattice (Q-span ∩ Z^8)
    from sympy.matrices.normalforms import hermite_normal_form
    # compute saturation via kernel of complement: use rational row reduction
    return B
Bplus=int_basis(P-sp.eye(8)); Bminus=int_basis(P+sp.eye(8))
print('B+ cols',Bplus.T.tolist()); print('B- cols',Bminus.T.tolist())
# check saturation: gcd of maximal minors ==1
def sat(B):
    import functools,math
    ms=[int(B.extract(list(rows),list(range(B.shape[1]))).det()) for rows in itertools.combinations(range(8),B.shape[1])]
    return functools.reduce(math.gcd,[abs(m) for m in ms])
print('saturation index + :',sat(Bplus),' - :',sat(Bminus))
# Omega restricted
print('Omega on H+',(Bplus.T*Omega*Bplus).tolist()); print('Omega on H-',(Bminus.T*Omega*Bminus).tolist())
# periods on H-
Pm=sp.simplify(per*Bminus); Pp=sp.simplify(per*Bplus)
print('periods on H+ (should be 0):',Pp.tolist())
# find integer X (4x4) on H- with Pm X = lam Pm, lam=a+b*phi
X=sp.Matrix(4,4,sp.symbols('x0:16')); a,b=sp.symbols('a b')
lam=a+b*phi
eqs=sp.expand(Pm*X-lam*Pm)
lin=[]
for e in eqs:
    e=sp.expand(e); e0=e.subs(s5,0); e1=sp.expand((e-e0)/s5)
    lin+= [e0,e1]
sol=sp.linsolve(lin,list(X)+[a,b])
print(sol)
