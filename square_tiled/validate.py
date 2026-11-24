import sympy as sp, json
from ohom import Origami
O=Origami([7,9,1,8,6,4,5,0,3,2],[1,5,6,4,3,0,8,9,2,7],[8,6,5,7,9,2,1,3,0,4])
H=O.hol();print('genus*2',O.g2,'sh',O.sh,'sv',O.sv)
Dh=sp.Matrix([[1,O.sh],[0,1]]);Dv=sp.Matrix([[1,0],[-O.sv,1]])
print('Th hol-consistent',H*O.Th==Dh*H,'Tv hol-consistent',H*O.Tv==Dv*H,'tau hol=-I',H*O.P==-H)
print('tau^2=I',O.P*O.P==sp.eye(8),'tau commutes',O.P*O.Th==O.Th*O.P,O.P*O.Tv==O.Tv*O.P)
Bp,Bm=O.blocks();print('dim H+',Bp.shape[1],'dim H-',Bm.shape[1],'hol on H+ zero',H*Bp==sp.zeros(2,Bp.shape[1]))
