import numpy as np
from sympy import *
import sympy as sp
import matplotlib.pyplot as plt

#-------------------------------------------------------------
#Given Data
#-------------------------------------------------------------
a= 100 #mm
b = 200 #mm
h = 20 #mm
p = 100 #Mpa
E = 210000 #Mpa
nu = 0.3

#-------------------------------------------------------------
#Analytical Method
#-------------------------------------------------------------

#Variables
fr,r,a0,b0,h0,Emod,nu0,C1,C2, P = symbols ('fr r a0 b0 h0 Emod nu0 C1 C2 P',real=True)

#Equation
ur_ = -((1-nu0**2)/(Emod*r))*integrate(r*integrate(fr,r),r)+C1*r/2+C2/r
Strain_rr = diff(ur_,r)
Strain_phiphi = ur_/r
Sigma_rr_ = (Emod/(1-nu0**2))*(Strain_rr+nu0*Strain_phiphi)
Sigma_phiphi  = (Emod/(1-nu0**2))*(nu0*Strain_rr+Strain_phiphi)
unknowns = (C1,C2)

#Boundary Conditions
eqns = [ur_.subs(r,a),Sigma_rr_.subs(r,b)+P]

#Solution
sol = solve(eqns,unknowns)
ur = ur_.subs(sol)
Sigma_rr = Sigma_rr_.subs(sol)

r_values = np.linspace(a,b,20)

deflection = []
Normal_stress = []

for ri in r_values:
    ur_for_disc = ur.subs({fr:0,r:ri,Emod:E,nu0:nu,P:p})
    deflection.append(ur_for_disc)
    Normal_stress_for_disc = Sigma_rr.subs({fr:0,r:ri,Emod:E,nu0:nu,P:p})
    Normal_stress.append(Normal_stress_for_disc)

plt.figure()
plt.plot(r_values,deflection)
plt.xlabel("Radius Values in mm",fontsize =16)
plt.ylabel("Radial Displacement in mm",fontsize =16)
plt.title("Radial Displacement in an Axisymmetric Disc",fontsize =16)
plt.grid()


plt.figure()
plt.plot(r_values,Normal_stress)
plt.xlabel("Radius Values in mm",fontsize =16)
plt.ylabel("Radial Stress in MPa",fontsize =16)
plt.title("Radial Stress in an Axisymmetric Disc",fontsize =16)
plt.grid()

#-------------------------------------------------------------
#FEM Method
#-------------------------------------------------------------

#Constant Thickness
#-------------------------------------------------------------

r,r1,r2,h0 = symbols('r r1 r2 h0',real=True)
num_elements = 30
num_nodes = num_elements+1
r_values_FEM = np.linspace(a,b,num_nodes)

K_global = np.zeros((num_nodes,num_nodes))
F_global = np.zeros(num_nodes)

D = (E/(1-nu**2))*np.array([[1,nu],[nu,1]])

Ke = []

for i in range(num_elements):
    r_1 = r_values_FEM[i]
    r_2 = r_values_FEM[i+1]
    Ne1 = (r2-r)/(r2-r1)
    Ne2 = (r-r1)/(r2-r1)
    B = Matrix([[diff(Ne1,r),diff(Ne2,r)],[Ne1/r,Ne2/r]])
    K = B.T*D*B*h0*r

    Ke_Matrix = 2*np.pi*integrate(K,(r,r1,r2))
    Ke = Ke_Matrix.subs({r1:r_1,r2:r_2,h0:h})
    K_global[i:i+2,i:i+2]+=Ke


#Boundary Condition

K_global[0,:] = 0
K_global[0,0] = 1
F_global[0] = 0
F_global[-1] = -p*h*2*np.pi*b

u_r = np.linalg.solve(K_global,F_global)

Normal_stress_rr = np.zeros(num_nodes)

for i in range (num_elements):
    du_r_dr = (u_r[i+1]-u_r[i])/(r_values_FEM[i+1]-r_values_FEM[i])

    sig_rr = (E/(1-nu**2))*(du_r_dr+(nu*u_r[i]/r_values_FEM[i]))
    Normal_stress_rr[i]= sig_rr

Normal_stress_rr[-1]= -p

plt.figure()
plt.plot(r_values_FEM,u_r,'*', color='blue',label='Radial Displacement FEM Method')
plt.plot(r_values,deflection, color= 'red', label = 'Radial Displacement Analytical Method')
plt.xlabel("Radius Values in mm",fontsize =16)
plt.ylabel("Radial Displacement \nin mm",fontsize =16)
plt.title('Radial Displacement:FEM Method v/s Analytical Method',fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Rdispanalytical.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(r_values_FEM,Normal_stress_rr,'*',color='blue',label='Normal Stress FEM Method')
plt.plot(r_values,Normal_stress,color= 'red', label = 'Normal Stress Analytical Method')
plt.xlabel("Radius Values in mm",fontsize =16)
plt.ylabel("Normal Stress in Mpa",fontsize =16)
plt.title('Normal Stress:FEM Method v/s Analytical Method',fontsize =16)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('Rstressanalytical.png', dpi=300,bbox_inches='tight')
plt.show()

#%%Varying Thickness
#-------------------------------------------------------------

K_global_2 = np.zeros((num_nodes,num_nodes))
F_global_2 = np.zeros(num_nodes)

D_2 = (E/(1-nu**2))*np.array([[1,nu],[nu,1]])


Ke_2 = []

for i in range(num_elements):
    r_1 = r_values_FEM[i]
    r_2 = r_values_FEM[i+1]
    Ne1 = (r2-r)/(r2-r1)
    Ne2 = (r-r1)/(r2-r1)
    B_2 = Matrix([[diff(Ne1,r),diff(Ne2,r)],[Ne1/r,Ne2/r]])
    K_2 = B.T*D*B*((h0*(r-a)/(b-a))+h0)*r

    Ke_Matrix_2 = 2*np.pi*integrate(K_2,(r,r1,r2))
    Ke_2 = Ke_Matrix_2.subs({r1:r_1,r2:r_2,h0:h})
    K_global_2[i:i+2,i:i+2]+=Ke_2


#Boundary Condition

K_global_2[0,:] = 0
K_global_2[0,0] = 1
F_global_2[0] = 0
F_global_2[-1] = -p*2*h*2*np.pi*b

u_r_2 = np.linalg.solve(K_global_2,F_global_2)

Normal_stress_rr_2 = np.zeros(num_nodes)

for i in range (num_elements):
    du_r_dr_2 = (u_r_2[i+1]-u_r_2[i])/(r_values_FEM[i+1]-r_values_FEM[i])

    sig_rr = (E/(1-nu**2))*(du_r_dr_2+(nu*u_r_2[i]/r_values_FEM[i]))
    Normal_stress_rr_2[i]= sig_rr

Normal_stress_rr_2[-1]= -p
    
plt.figure()
plt.plot(r_values_FEM,u_r_2)
plt.xlabel("Radius Values",fontsize =16)
plt.ylabel("Radial Displacement",fontsize =16)
plt.title('Radial Displacement $u_r$ for Axisymmteric \nDisc for Varying Thickness',fontsize =16)
plt.grid()
plt.tight_layout()
plt.savefig('Vthickness.png', dpi=300,bbox_inches='tight')

plt.figure()
plt.plot(r_values_FEM,Normal_stress_rr_2)
plt.xlabel("radius Values")
plt.ylabel("Normal Stress")
plt.grid()
plt.title('Normal Stress $sigma_{rr}$ for Axisymmteric Disc for Varying Thickness')
plt.show()





    


