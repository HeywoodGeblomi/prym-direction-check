"""Toy model that reproduces the recorded (A, M_plus) exactly.
H1+ basis c0,c1 (horizontal-type cores), c2,c3 (vertical-type cores), form Om with Om(c0,c2)=Om(c1,c3)=a (a=-1).
Letter (i,k): homology = k-fold transvection x -> x + k*Om(c_i,x)*c_i ; derivative = A_h(k) if i in {0,1}, A_v(k) if i in {2,3}.
Word applied left to right (first letter acts first)."""
import numpy as np
def Om(a=-1):
    O=np.zeros((4,4),dtype=object); O[0,2]=a;O[2,0]=-a;O[1,3]=a;O[3,1]=-a; return O
def transvection(i,k,O):
    M=np.eye(4,dtype=object).astype(object)
    for j in range(4):  # image of basis vector j
        M[i,j]+=k*O[i,j]
    return M
Ah=lambda k:np.array([[1,0],[k,1]],dtype=object); Av=lambda k:np.array([[1,k],[0,1]],dtype=object)
def word_action(w,O):
    M=np.eye(4,dtype=object); D=np.eye(2,dtype=object)
    for i,k in w:
        M=transvection(i,k,O).dot(M); D=(Ah(k) if i in (0,1) else Av(k)).dot(D)
    return M,D
if __name__=='__main__':
    w=[(1,3),(2,4),(2,4),(3,2),(0,2)]
    M,D=word_action(w,Om(-1))
    print('M_plus from toy:',M.tolist()); print('A from toy:',D.tolist())
    print('matches recorded M_plus:',M.tolist()==[[-15,0,-2,0],[0,1,0,-3],[8,0,1,0],[0,2,0,-5]])
    print('matches recorded A:',D.tolist()==[[31,10],[65,21]])
    # sign-consistency test: on a single square torus a twist acts on H1 exactly as its derivative.
    # Block (c1,c3) is a torus seeing letters (1,3) then (3,2): derivative-consistent action would be Av(2)Ah(3)
    b1=np.array([[M[1,1],M[1,3]],[M[3,1],M[3,3]]]); b0=np.array([[M[0,0],M[0,2]],[M[2,0],M[2,2]]])
    print('block(c0,c2)',b0.tolist(),'trace',b0.trace(),' block(c1,c3)',b1.tolist(),'trace',b1.trace())
    c0=Ah(2).dot(Av(8)); c1=Av(2).dot(Ah(3))
    print('derivative-consistent torus blocks: ',c0.tolist(),'trace',c0.trace(),';',c1.tolist(),'trace',c1.trace())
    import math
    l=lambda t:math.log((abs(t)+math.sqrt(t*t-4))/2)
    print('ratio with recorded blocks',l(b0.trace())/l(b1.trace()),' ratio with sign-consistent blocks',l(c0.trace())/l(c1.trace()))
    # with a=+1
    M2,_=word_action(w,Om(1)); print('a=+1 gives',M2.tolist())
