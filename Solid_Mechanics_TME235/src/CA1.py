import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import calfem.vis_mpl as cfv
import calfem.core as cfc
from scipy.sparse import csr_matrix

#-------------------------------------------------------------
#Given Data
#-------------------------------------------------------------

E = 200000000000 #in Pascal
nu = 0.3
yield_stress = 450000000 #in Pascal
density = 7800 #kg/m3
g = 9.81 #m/s2
m = 120 #kg
h= 0.05 #m
b = 0.05 #m
K = 5/6
L1 = 2 #m
L2 = 0.25 #m
q0 = h*b*density*g
p_load = m*g
MI = (b*h**3)/12
G = E/(2*(1+nu))
ptype = 1
Dmat = cfc.hooke(ptype,E,nu)
ep = [ptype,h]

#-------------------------------------------------------------
#%%Euler Bernouli Theorem
#-------------------------------------------------------------
#Variables
p,q,x,L,Emod,I,C1,C2,C3,C4,z = symbols ('p q x L Emod I C1 C2 C3 C4 z',real=True)

# Equation
w1 = (integrate(-q,x,x,x,x)+C1*x**3/6+C2*x**2/2+C3*x+C4)/(Emod*I)
w1prime = diff(w1,x)
M1 = -simplify(Emod*I*diff(w1prime,x))
V1 = simplify (diff(M1,x))

#Boundary Conditions
eqns1 = [ w1.subs(x,0), w1prime.subs(x,0),M1.subs(x,L), V1.subs(x,L)+p] #Boundary Conditions for distributed load

#Solution
unknowns1 = (C1,C2,C3,C4)
sol1 = solve (eqns1,unknowns1)
w_total = simplify (w1.subs(sol1))

#Normal stress σxx

Normal_stress = -simplify(Emod*z*(diff(w1prime.subs(sol1),x)))

#Solving for Beam 1 where Length of the beam is 2 m

x1_values = np.linspace(0,L1,20)
deflection_for_L1 = []
Normal_stresses_for_L1 = []

for xi in x1_values:
    w_for_L1 = w_total.subs({p:p_load,L:L1,q:q0,I:MI,Emod:E,x:xi})
    deflection_for_L1.append(w_for_L1)
    stress= Normal_stress.subs({p:p_load,L:L1,q:q0,I:MI,Emod:E,x:xi,z:-h/2})
    Normal_stresses_for_L1.append(stress)

#Von Mises Stress
Von_Mises_stress_EB1 = np.abs(Normal_stresses_for_L1)

EB_yielding_L1 = False

for i, xi in enumerate(x1_values):
    if Von_Mises_stress_EB1[i]>yield_stress:
        print(f"The beam starts yielding at position x = {x1_values[i]:.4f} meters with von Mises stress = {Von_Mises_stress_EB1[i]:.2f}Pa")
        EB_yielding_L1 = True
        break
if not EB_yielding_L1:
    print('Beam 1 is not yielding as per Euler-Bernoulli Beam Theory')

#Solving for Beam 2 where Length of the beam is 0.25 m

x2_values_EB = np.linspace(0,L2,20)
deflection_for_L2 = []
Normal_stresses_for_L2 = []

for i in x2_values_EB:
    w_for_L2 = w_total.subs({p:p_load,L:L2,q:q0,I:MI,Emod:E,x:i})
    deflection_for_L2.append(w_for_L2)
    stress2= Normal_stress.subs({p:p_load,L:L2,q:q0,I:MI,Emod:E,x:i,z:-h/2})
    Normal_stresses_for_L2.append(stress2)

#Von Mises Stress
Von_Mises_stress_EB_2 = np.abs(Normal_stresses_for_L2)

EB_yielding_L2 = False

for i, xi in enumerate(x2_values_EB):
    if Von_Mises_stress_EB_2[i]>yield_stress:
        print(f"The beam 2 starts yielding at position x = {x2_values[i]:.4f} meters with von Mises stress = {Von_Mises_stress_EB_2[i]:.2f}Pa")
        EB_yielding_L2 = True
        break
if not EB_yielding_L2:
    print('Beam 2 is not yielding as per Euler-Bernoulli Beam Theory')
    
#Graphs
plt.figure()
plt.plot(x1_values,deflection_for_L1,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Euler-Bernoulli: Deflection of the Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli Deflection B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x1_values,Normal_stresses_for_L1,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress in Pa",fontsize =16)
plt.title("Euler-Bernoulli: Normal Stress in Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli NS B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x1_values,Von_Mises_stress_EB1,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Von-Mises Stress in Pa",fontsize =16)
plt.title("Euler-Bernoulli: Von-Mises Stress in Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli VM B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_EB,deflection_for_L2,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Euler-Bernoulli: Deflection of the Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli Deflection B2.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_EB,Normal_stresses_for_L2,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress in Pa",fontsize =16)
plt.title("Euler-Bernoulli: Normal Stress in Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli NS B2.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_EB,Von_Mises_stress_EB_2,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Von-Mises Stress in Pa",fontsize =16)
plt.title("Euler-Bernoulli: Von-Mises Stress in Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Euler-Bernoulli VM B2.png', dpi=300,bbox_inches='tight')

#-------------------------------------------------------------
#%%Timoshenko Theorem
#-------------------------------------------------------------
#Variables
p,q,x,L,Emod,I,C1,C2,C3,C4,z,Gmod,phi,K_s,A = symbols ('p q x L Emod I C1 C2 C3 C4 z Gmod phi K_s A',real=True)

#Body Weight Equation
phi1 = (-q*x**3/6+C1*x**2/2+C2*x+C3)/(Emod*I)
M1 = -simplify(Emod*I*diff(phi1,x))
V1 = simplify(diff(M1,x))
w1 = (M1/(Gmod*K_s*A))+integrate(phi1,x)+C4

#Boundary Conditions
eqns1 = [w1.subs(x,0),phi1.subs(x,0),M1.subs(x,L),V1.subs(x,L)+p]

#Solution
unknowns1 = (C1,C2,C3,C4)
sol1 = solve(eqns1,unknowns1)
w_total_TT = simplify(w1.subs(sol1))

#Normal Stress
Normal_stress_TT = -Emod*z*diff(phi1.subs(sol1),x)

#Shear Stress
shear_stress = G*(-phi1.subs(sol1)+diff(w1.subs(sol1),x))

#Solving for Beam 1 where length is 2m
x1_values_TT = np.linspace(0,L1,20)
deflection_for_L1_TT = []
Normal_stress_for_L1_TT = []
Shear_stress_for_L1_TT = []

for xi in x1_values_TT:
    w_for_L1_TT = w_total_TT.subs({p:p_load,L:L1,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b})
    deflection_for_L1_TT.append(w_for_L1_TT)
    sigma_for_L1_TT = Normal_stress_TT.subs({p:p_load,L:L1,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b,z:-h/2})
    Normal_stress_for_L1_TT.append(sigma_for_L1_TT)
    s_stress_for_L1_TT = shear_stress.subs({p:p_load,L:L1,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b,z:-h/2})
    Shear_stress_for_L1_TT.append(s_stress_for_L1_TT)

#Von-Mises Stress
Von_Mises_stress_TT_1=[]

for i in range (len(x1_values_TT)):
    Von_Mises_TT_1 = ((Normal_stress_for_L1_TT[i])**2+(3*(Shear_stress_for_L1_TT[i])**2))**0.5
    Von_Mises_stress_TT_1.append(Von_Mises_TT_1)

TT_yielding_L1 = False

for i, xi in enumerate(x1_values_TT):
    if Von_Mises_stress_TT_1[i]>yield_stress:
        print(f"The beam 1 starts yielding at position x = {x1_values_TT[i]:.4f} meters with von Mises stress = {Von_Mises_stress_TT_1[i]:.2f}Pa as per Timoshenko Beam Theory")
        TT_yielding_L1 = True
        break
if not TT_yielding_L1:
    print('Beam 1 is not yielding as per Timoshenko Beam Theory')

#Solving for Beam 2 where length is 0.25m
x2_values_TT = np.linspace(0,L2,20)
deflection_for_L2_TT = []
Normal_stress_for_L2_TT = []
Shear_stress_for_L2_TT=[]

for xi in x2_values_TT:
    w_for_L2_TT = w_total_TT.subs({p:p_load,L:L2,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b})
    deflection_for_L2_TT.append(w_for_L2_TT)
    sigma_for_L2_TT = Normal_stress_TT.subs({p:p_load,L:L2,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b,z:-h/2})
    Normal_stress_for_L2_TT.append(sigma_for_L2_TT)
    s_stress_for_L2_TT = shear_stress.subs({p:p_load,L:L2,q:q0,I:MI,Emod:E,x:xi,Gmod:G,K_s:K,A:h*b,z:-h/2})
    Shear_stress_for_L2_TT.append(s_stress_for_L2_TT)

#Von-Mises Stress
Von_Mises_stress_TT_2=[]

for i in range (len(x2_values_TT)):
    Von_Mises_TT_2 = ((Normal_stress_for_L2_TT[i])**2+(3*(Shear_stress_for_L2_TT[i])**2))**0.5
    Von_Mises_stress_TT_2.append(Von_Mises_TT_2)

TT_yielding_L2 = False

for i, xi in enumerate(x2_values_TT):
    if Von_Mises_stress_TT_2[i]>yield_stress:
        print(f"The beam 2 starts yielding at position x = {x2_values_TT[i]:.4f} meters with von Mises stress = {Von_Mises_stress_TT_2[i]:.2f}Pa as per Timoshenko Beam Theory")
        TT_yielding_L2 = True
        break
if not TT_yielding_L2:
    print('Beam 2 is not yielding as per Timoshenko Beam Theory')
    
#Graphs
plt.figure()
plt.plot(x1_values_TT,deflection_for_L1_TT,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Deflection of the Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model Deflection B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x1_values_TT,Normal_stress_for_L1_TT,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress in Pa",fontsize =16)
plt.title("Normal Stress in Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model NS B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x1_values_TT,Von_Mises_stress_TT_1,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Von-Mises Stress in Pa",fontsize =16)
plt.title("Von-Mises Stress in Beam 1 v/s Position along the beam",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model VM B1.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_TT,deflection_for_L2_TT,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Deflection of the Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model Deflection B2.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_TT,Normal_stress_for_L2_TT,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress in Pa",fontsize =16)
plt.title("Normal Stress in Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model NS B2.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(x2_values_TT,Von_Mises_stress_TT_2,color='b',label='Euler-Bernouli Beam Model',linewidth=2)
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Von-Mises Stress in Pa",fontsize =16)
plt.title("Von-Mises Stress in Beam 2 v/s Position along the beam",fontsize =16)
plt.grid(True)
plt.tight_layout()
plt.savefig('Timoshenko Beam Model VM B2.png', dpi=300,bbox_inches='tight')
plt.show()

#%% FEM Method
# Quad Mesh
def ex_ey_quadmesh(p1,p2,nelx,nely,ndofs):
    xv=np.linspace(p1[0],p2[0],nelx+1)
    yv=np.linspace(p1[1],p2[1],nely+1)

    nel  = nelx*nely;
    Ex=np.zeros((nel,4))
    Ey=np.zeros((nel,4))
    for m in range(0,nely):
        for n in range(0,nelx):
            Ex[n+m*nelx,0]=xv[n]
            Ex[n+m*nelx,1]=xv[n+1]
            Ex[n+m*nelx,2]=xv[n+1]
            Ex[n+m*nelx,3]=xv[n]
        #
            Ey[n+m*nelx,0]=yv[m]
            Ey[n+m*nelx,1]=yv[m]
            Ey[n+m*nelx,2]=yv[m+1]
            Ey[n+m*nelx,3]=yv[m+1]
            
    return Ex, Ey

def edof_quadmesh(nelx,nely,ndofs):
    Edof=np.zeros((nelx*nely,4*ndofs),'i')
    for m in range(0,nely):
        for n in range(0,nelx):
            Edof[n+m*nelx,0]=n*ndofs+1+m*(nelx+1)*ndofs
            Edof[n+m*nelx,1]=n*ndofs+2+m*(nelx+1)*ndofs
            Edof[n+m*nelx,2]=(n+1)*ndofs+1+m*(nelx+1)*ndofs
            Edof[n+m*nelx,3]=(n+1)*ndofs+2+m*(nelx+1)*ndofs
        #
            Edof[n+m*nelx,4]=(n+1)*ndofs+1+(m+1)*(nelx+1)*ndofs
            Edof[n+m*nelx,5]=(n+1)*ndofs+2+(m+1)*(nelx+1)*ndofs
            Edof[n+m*nelx,6]=n*ndofs+1+(m+1)*(nelx+1)*ndofs      
            Edof[n+m*nelx,7]=n*ndofs+2+(m+1)*(nelx+1)*ndofs
    return Edof

def B1B2B3B4_quadmesh(nelx,nely,ndofs):
    #lower boundary, dofs
    B1=np.linspace(1,(nelx+1)*ndofs,(nelx+1)*ndofs)
    B1=B1.astype(int)
    B2=np.zeros(((nely+1)*ndofs),'i')
    nn=0
    for n in range(0,nely+1):
        B2[nn]=(nelx+1)*ndofs*(n+1)-1
        if ndofs>1:
            B2[nn+1]=(nelx+1)*ndofs*(n+1)+0
        nn=nn+ndofs

    B3=np.linspace(1,(nelx+1)*ndofs,(nelx+1)*ndofs)+(nelx+1)*ndofs*nely
    B3=B3.astype(int)

    B4=np.zeros(((nely+1)*ndofs),'i')
    nn=0
    for n in range(0,nely+1):
        B4[nn]=(nelx+1)*ndofs*n+1
        if ndofs>1:
            B4[nn+1]=(nelx+1)*ndofs*n+2
        nn=nn+ndofs
   
    P1=np.zeros((2),'i'); P2=np.zeros((2),'i'); P3=np.zeros((2),'i'); P4=np.zeros((2),'i')
    for m in range(0,2):
        P1[m]=B1[m]
        P2[m]=B2[m]
        P4[m]=B3[m]
    P3[0]=B3[-1]-1
    P3[1]=B3[-1]
    return B1,B2,B3,B4,P1,P2,P3,P4

def quadmesh(p1,p2,nelx,nely,ndofs):
    Ex, Ey=ex_ey_quadmesh(p1,p2,nelx,nely,ndofs)
    Edof=edof_quadmesh(nelx,nely,ndofs)
    B1,B2,B3,B4,P1,P2,P3,P4=B1B2B3B4_quadmesh(nelx,nely,ndofs)
    return Ex,Ey,Edof,B1,B2,B3,B4,P1,P2,P3,P4

# %%Beam 1 Calfem

p1=np.array([0.,0.])
p2 = np.array([L1,h])
nelx=640;nely=16;
ndof_per_node = 2
nnode = (nelx+1)*(nely+1)
nDofs = ndof_per_node*nnode
Ex,Ey,Edof,B1,B2,B3,B4,P1,P2,P3,P4=quadmesh(p1,p2,nelx,nely,ndof_per_node)
cfv.eldraw2(Ex,Ey)
K= np.zeros([nDofs,nDofs])
F_global = np.zeros((nDofs,1))

##Boundary Conditions
bc = B4
bcVal=0*np.ones(np.size(bc))
force_per_node = -p_load/(int(len(B2)/2))
F_global[B2[::2]]+=force_per_node
eq = np.array([[0],[-density*g]])

for eltopo, elx,ely in zip(Edof,Ex,Ey):
    Ke,fe= cfc.planqe(elx,ely,ep,Dmat,eq)
    cfc.assem(eltopo,K,Ke,F_global,fe)
    
Ks = csr_matrix(K,shape=(nDofs,nDofs))
a,r=cfc.spsolveq(Ks,F_global,bc,bcVal)

deflection1 = a[1::2][:nelx+1]
Beam_length = np.linspace(0,L1,len(deflection1))
length_1 = np.linspace(0,L1,nelx)

Lower_Beam_Stress_FEM = []
Von_mises_stress_FEM1 = []
    
for eltopo, elx, ely in zip(Edof, Ex, Ey):
    lower_node_element = np.any(np.isin(eltopo, B1))
    
    if lower_node_element:
        elx = np.array(elx).flatten()
        ely = np.array(ely).flatten()
        element = np.array(a[eltopo - 1]).flatten()
        es, et = cfc.planqs(elx, ely, ep, Dmat, element)
        sigma_xx = es[0]
        sigma_yy = es[1]
        tau_xy = es[2]
        Lower_Beam_Stress_FEM.append(sigma_xx)
        VM_stress = np.sqrt((sigma_xx**2)-(sigma_xx*sigma_yy)+(sigma_yy**2)+(3*tau_xy**2))
        Von_mises_stress_FEM1.append(VM_stress)

#%% Beam 2 Calfem

p1=np.array([0.,0.])
p2 = np.array([L2,h])
nelx=160;nely=32;
ndof_per_node = 2
nnode = (nelx+1)*(nely+1)
nDofs = ndof_per_node*nnode
Ex,Ey,Edof,B1,B2,B3,B4,P1,P2,P3,P4=quadmesh(p1,p2,nelx,nely,ndof_per_node)
cfv.eldraw2(Ex,Ey)
K= np.zeros([nDofs,nDofs])
F_global = np.zeros((nDofs,1))

##Boundary Conditions
bc = B4
bcVal=0*np.ones(np.size(bc))
force_per_node = -p_load/(int(len(B2)/2))
F_global[B2[::2]]+=force_per_node
eq = np.array([[0],[-density*g]])

for eltopo, elx,ely in zip(Edof,Ex,Ey):
    Ke,fe= cfc.planqe(elx,ely,ep,Dmat,eq)
    cfc.assem(eltopo,K,Ke,F_global,fe)
    
Ks = csr_matrix(K,shape=(nDofs,nDofs))
a,r=cfc.spsolveq(Ks,F_global,bc,bcVal)

deflection_Beam2 = a[1::2][:nelx+1]
Beam_length_2 = np.linspace(0,L2,len(deflection_Beam2))
length_2 = np.linspace(0,L2,nelx)

Lower_Beam_Stress_FEM_2 = []
Von_mises_stress_FEM_2 = []
    
for eltopo, elx, ely in zip(Edof, Ex, Ey):
    lower_node_element2 = np.any(np.isin(eltopo, B1))
    
    if lower_node_element2:
        elx = np.array(elx).flatten()
        ely = np.array(ely).flatten()
        element = np.array(a[eltopo - 1]).flatten()
        es, et = cfc.planqs(elx, ely, ep, Dmat, element)
        sigma_xx = es[0]
        sigma_yy = es[1]
        tau_xy = es[2]
        Lower_Beam_Stress_FEM_2.append(sigma_xx)
        VM_stress2 = np.sqrt((sigma_xx**2)-(sigma_xx*sigma_yy)+(sigma_yy**2)+(3*tau_xy**2))
        Von_mises_stress_FEM_2.append(VM_stress2)
        

#%% Abaqus Data

import pandas as pd
##Beam1
#Deflection,Normal Stress,Von-Mises Stress
data_Beam_1 = pd.read_csv('../data/Consolidated_data_Beam1.csv')
Nodes_FEM_abaqus_B1 = data_Beam_1['Node']
Deflection_FEM_abaqus_B1 = data_Beam_1['U2']
Normal_stress_FEM_abaqus_B1 = data_Beam_1['Normal Stress']
Von_Mises_stress_FEM_abaqus_B1 = data_Beam_1['Von-Mises Stress']

#Horizontal Displacement Beam 1
Data_Beam_1_disp = pd.read_csv('../data/U1_wrt_z_beam1.csv')
Nodes_Beam1_h = Data_Beam_1_disp['Node']
Horizontal_disp_Beam_1 = Data_Beam_1_disp['U1']

##Beam2
#Deflection,Normal Stress,Von-Mises Stress
data_Beam_2=  pd.read_csv('../data/Consolidated Table Beam2.csv')
Nodes_FEM_abaqus_B2 = data_Beam_2['Node']
Deflection_FEM_abaqus_B2 = data_Beam_2['U2']
Normal_stress_FEM_abaqus_B2 = data_Beam_2['Normal Stress']
Von_Mises_stress_FEM_abaqus_B2 = data_Beam_2['Von-Mises Stress']

#Horizontal Displacement Beam 1
Data_Beam_2_disp = pd.read_csv('../data/U1_wrt_z_beam2.csv')
Nodes_Beam2_h = Data_Beam_2_disp['Node']
Horizontal_disp_Beam_2 = Data_Beam_2_disp['U1']


 
#%%Combining Graphs
#For deflection in Beam 1
plt.figure()
plt.plot(x1_values,deflection_for_L1,color='b',label='Euler-Bernouli Beam Model')
plt.plot(x1_values_TT,deflection_for_L1_TT,color='r',linestyle='--',label='Timoshenko Beam Model')
plt.plot(Beam_length,deflection1,color='g',label = 'FEM Model using CalFEM')
plt.plot(Nodes_FEM_abaqus_B1,Deflection_FEM_abaqus_B1,color='m',label='FEM Model using Abaqus')
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Deflection of the Beam 1 v/s Position along the beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Deflection EBvTT B1.png', dpi=300,bbox_inches='tight')

#For Normal Stress in Beam 1
plt.figure()    
plt.plot(x1_values,Normal_stresses_for_L1,color='b',label='Euler-Bernouli Beam Model')
plt.plot(x1_values_TT,Normal_stress_for_L1_TT,color='r',linestyle='--',label='Timoshenko Beam Model')
plt.plot(length_1,Lower_Beam_Stress_FEM,color='g',label='FEM Model using CalFEM')
plt.plot(Nodes_FEM_abaqus_B1,Normal_stress_FEM_abaqus_B1,color='m',label='FEM Model using Abaqus')
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress along \nthe beam in Pascal",fontsize =16)
plt.title("Normal Stress of the Beam 1 v/s Position along the beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('NS EBvTT B1.png', dpi=300,bbox_inches='tight')

#For deflection in Beam 2
plt.figure()
plt.plot(x2_values_EB,deflection_for_L2,color='b',label='Euler-Bernouli Beam Model')
plt.plot(x2_values_TT,deflection_for_L2_TT,color='r',linestyle='--',label='Timoshenko Beam Model')
plt.plot(Beam_length_2,deflection_Beam2,color='g',label='FEM Model using CalFEM')
plt.plot(Nodes_FEM_abaqus_B2,Deflection_FEM_abaqus_B2 ,color='m',label='FEM Model using Abaqus')
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Total Deflection in m",fontsize =16)
plt.title("Deflection of the Beam 2 v/s Position along the beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Deflection EBvTT B2.png', dpi=300,bbox_inches='tight')


#For Normal stress in Beam 2
plt.figure()    
plt.plot(x2_values_EB,Normal_stresses_for_L2,color='b',label='Euler-Bernouli Beam Model')
plt.plot(x2_values_TT,Normal_stress_for_L2_TT,color='r',linestyle='--',label='Timoshenko Beam Model')
plt.plot(length_2,Lower_Beam_Stress_FEM_2,color='g',label='FEM Model using CalFEM')
plt.plot(Nodes_FEM_abaqus_B2,Normal_stress_FEM_abaqus_B2,color='m',label='FEM Model using Abaqus')
plt.xlabel("Position along the Beam in m",fontsize =16)
plt.ylabel("Normal Stress along \nthe beam in Pascal",fontsize =16)
plt.title("Normal Stress of the Beam 2 v/s Position along the beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Ns EBvTT B2.png', dpi=300,bbox_inches='tight')
plt.show()

# #Horizontal Displacement Beam 1
plt.figure()
plt.plot(Nodes_Beam1_h,Horizontal_disp_Beam_1,label='FEM Model using Abaqus')
plt.xlabel("Height of the Beam in m",fontsize =16)
plt.ylabel("Horizontal Displacement \nin m",fontsize =16)
plt.title("Horizontal Displacement in Beam 1 v/s Height of the Beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Horizontal Displacement B1.png', dpi=300,bbox_inches='tight')

# #Horizontal Displacement Beam 2
plt.figure()
plt.plot(Nodes_Beam2_h,Horizontal_disp_Beam_2,label='FEM Model using Abaqus')
plt.xlabel("Height of the Beam in m",fontsize =16)
plt.ylabel("Horizontal Displacement \nin m",fontsize =16)
plt.title("Horizontal Displacement in Beam 2 v/s Height of the Beam",fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Horizontal Displacement B2.png', dpi=300,bbox_inches='tight')
plt.show()
