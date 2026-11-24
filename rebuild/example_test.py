import sympy as sp, json
phi=(1+sp.sqrt(5))/2; mu=phi**2
H=lambda k: sp.Matrix([[1,k*mu**2],[0,1]]); V=lambda k: sp.Matrix([[1,0],[-k,1]])  # normalized TV generators
w=[(1,3),(2,4),(2,4),(3,2),(0,2)]
out={}
readings={'even index=H, odd=V':lambda i:'H' if i%2==0 else 'V','even index=V, odd=H':lambda i:'V' if i%2==0 else 'H',
          'index<2=H, else V':lambda i:'H' if i<2 else 'V','index<2=V, else H':lambda i:'V' if i<2 else 'H'}
for name,f in readings.items():
    D=sp.eye(2)
    for i,k in w: D=D*(H(k) if f(i)=='H' else V(k))
    out[name]=str(sp.nsimplify(sp.expand(D.trace())))
A=sp.Matrix([[31,10],[65,21]]); Mp=sp.Matrix([[-15,0,-2,0],[0,1,0,-3],[8,0,1,0],[0,2,0,-5]])
t=sp.symbols('t')
out['A_det']=int(A.det()); out['A_trace']=int(A.trace()); out['A_eigs']=[str(e) for e in A.eigenvals()]
out['Mplus_charpoly']=str(sp.factor(Mp.charpoly(t).as_expr())); out['Mplus_eigs']=[str(sp.simplify(e)) for e in Mp.eigenvals()]
out['rho_A_equals_(2+sqrt3)^3']=sp.simplify((26+15*sp.sqrt(3))-(2+sp.sqrt(3))**3)==0
print(json.dumps(out,indent=1))
