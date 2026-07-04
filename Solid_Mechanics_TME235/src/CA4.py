# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 22:36:13 2024

@author: Sankalp
"""

import numpy as np
from sympy import *
import matplotlib.pyplot as plt

#%% Given Data

L0 = 1000 #mm
h = 30 #mm
b = 30 #mm
E_steel = 210000 #MPa
nu_steel = 0.3
E_rubber = 20 #MPa
nu_rubber = 0.45
p0 = 10000 #N
q0 = 0
MI0 = b*h**3/12

#%% Pure stress vs stretch for rubber

G = E_rubber/(2*(1+nu_rubber))

lamb = np.linspace(0.1,1.9,100)

stress = []
for i in lamb:
    st = G*(i**2-1)
    stress.append(st)

plt.figure()
plt.plot(lamb,stress)
plt.xlabel("Stretch in Rubber ($\lambda$)", fontsize =16)
plt.ylabel("True Stress in MPa ($\sigma$)", fontsize =16)
plt.title("True Stress v/s Stretch", fontsize =16)
plt.grid(True)
plt.savefig("True_stress_v_stretch.png", dpi=300)
plt.show()

#%% Abaqus Data
import pandas as pd

Data = pd.read_csv('../data/Ca4_def_U2_data.csv')

nodes = Data['Node']
U_down = Data['U_down']
U_up = Data['U_up']

plt.figure()
plt.plot(nodes,U_up,'b',label='Deflection for 10kN force')
plt.plot(nodes,np.abs(U_down),'r',label='Deflection for -10kN force')
plt.xlabel("Distance along the length in mm", fontsize = 16)
plt.ylabel("Absolute Deflection in mm", fontsize = 16)
plt.legend()
plt.grid(True)
plt.title("Deflection in Beam", fontsize = 16)
plt.tight_layout()
plt.savefig("Deflection.png", dpi=300)
plt.show()

#%% Euler-Bernoulli beam model

#Variables
p,q,x,L,Emod,I,C1,C2,C3,C4,z = symbols ('p q x L Emod I C1 C2 C3 C4 z',real=True)

# Equation
w1 = (integrate(-q,x,x,x,x)+C1*x**3/6+C2*x**2/2+C3*x+C4)/(Emod*I)
w1prime = diff(w1,x)
M1 = -simplify(Emod*I*diff(w1prime,x))
V1 = simplify (diff(M1,x))

#Boundary Conditions upward force of 10kN
eqns1 = [ w1.subs(x,0), w1prime.subs(x,0),M1.subs(x,L), V1.subs(x,L)-p] 

#Boundary Conditions downward force of 10kN
eqns2 = [ w1.subs(x,0), w1prime.subs(x,0),M1.subs(x,L), V1.subs(x,L)+p] 

#Solution 1
unknowns = (C1,C2,C3,C4)
sol1 = solve (eqns1,unknowns)
w_total_up = simplify (w1.subs(sol1))

#Solution 2
unknowns = (C1,C2,C3,C4)
sol2 = solve (eqns2,unknowns)
w_total_down = simplify (w1.subs(sol2))

x_values = np.linspace(0, L0)
def_up_EB = []
def_down_EB = []

for i in x_values:
    w_up= w_total_up.subs({p:p0,L:L0,q:q0,I:MI0,Emod:E_steel,x:i})
    def_up_EB.append(w_up)
    w_down = w_total_down.subs({p:p0,L:L0,q:q0,I:MI0,Emod:E_steel,x:i})
    def_down_EB.append(w_down)


plt.figure()
plt.plot(nodes,U_up,color='r',label='Beam Deflection with Rubber Support')
plt.plot(x_values,def_up_EB,color='b',label='Euler-Bernoulli beam model w/o Rubber support')
plt.xlabel('Distance along the length of Beam (mm)',fontsize=16)
plt.ylabel('Deflection in mm',fontsize=16)
plt.title('Deflection comparison with and without rubber support for P = 10kN',fontsize=16)
plt.legend()
plt.tight_layout()
plt.grid(True)
plt.savefig('Deflection_up.png',dpi=300,bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(nodes,U_down,color='r',label='Beam Deflection with Rubber Support')
plt.plot(x_values,def_down_EB,color='b',label='Euler-Bernoulli beam model w/o Rubber support')
plt.xlabel('Distance along the length of Beam (mm)',fontsize=16)
plt.ylabel('Deflection in mm',fontsize=16)
plt.title('Deflection comparison with and without rubber support for P = -10kN',fontsize=16)
plt.legend()
plt.tight_layout()
plt.grid(True)
plt.savefig('Deflection_down.png',dpi=300,bbox_inches='tight')
plt.show()