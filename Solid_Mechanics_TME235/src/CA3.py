import numpy as np
import matplotlib.pyplot as plt
from sympy import *

r,A1,A2,A3,A4=symbols('r,A1,A2,A3,A4',real=True)


b=0.2
a=0.1
h=0.02
E=210e9
v=0.3
D=E*h**3/(12*(1-v**2))
print(D)
q= 10e6

w=integrate(1/r*integrate(r*integrate( 1/r*integrate(q*r/D,r),r),r),r)+A1*r**2*log(r/b) + A2*r**2 + A3*log(r/b) + A4
wprime=diff(w,r)

Mr=D*(-diff(wprime,r)-v/r*wprime ) #radial bending moment
Mphi=D*(-1/r*wprime -v*diff(wprime,r))#circumf bending moment

V=diff(Mr,r)+1/r*(Mr-Mphi)


eqns = [
    w.subs(r,b),
    Mr.subs(r,b),
    Mr.subs(r,a),
    V.subs(r,a)] #bc’s

unknowns=(A1,A2,A3,A4)

sol=solve(eqns,unknowns)
w_=simplify(w.subs(sol))

r_vals = np.linspace(a, b, 200)
w_vals=[]

for vals in r_vals: 
    w_vals.append(w_.subs(r,vals))
#%%
print(w_vals[0])
# plt.figure(figsize=(3.94,3.5), dpi=600)  
plt.figure()
plt.plot(r_vals, w_vals,'m--', label='deflection')
plt.xlabel('$r$ (m)')
plt.ylabel('$w(r)$ (m)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('Deflection.png', dpi=300,bbox_inches='tight')

#-------------------------------------------------------------------------
#%%Normal Stresses
#-------------------------------------------------------------------------
z = h / 2
w_double_prime=diff(wprime,r)
sigma_rr = (E/(1-v**2)) * (-z*w_double_prime.subs(sol) -z*v*wprime.subs(sol)/r)

sigma_rr_vals= []

for val in r_vals:
    sigma_rr_vals.append(sigma_rr.subs(r, val).evalf())  

print(np.max(sigma_rr_vals))
# Plot normal stress
# plt.figure(figsize=(3.94,3), dpi=600) 
plt.figure()
plt.plot(r_vals, sigma_rr_vals,'r--', label='$\sigma_{rr}$')
plt.xlabel('$r$ (m)')
plt.ylabel('$\sigma_{rr}$ (Pa)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('Normal_Stress.png', dpi=300,bbox_inches='tight')

#-------------------------------------------------------------------------
#%%Von Misses Stresses
#-------------------------------------------------------------------------

z = h / 2
w_double_prime=diff(wprime,r)
sigma_pp = (E/(1-v**2)) * (-v*z*w_double_prime.subs(sol) -z*wprime.subs(sol)/r)


sigma_pp_vals= []

for val in r_vals:
    sigma_pp_vals.append(sigma_pp.subs(r, val).evalf())  


sigma_em=[]
for i in range(len(sigma_pp_vals)):
    value = (1.5 * (1/9) * ((2 * sigma_rr_vals[i] - sigma_pp_vals[i])**2 +
        (2 * sigma_pp_vals[i] - sigma_rr_vals[i])**2 +
        (sigma_rr_vals[i] + sigma_pp_vals[i])**2
    ))**0.5
    sigma_em.append(value)

print(sigma_em[0])
plt.figure()
plt.plot(r_vals, sigma_em,'r--', label='$\sigma_{eM}$')
plt.xlabel('radius')
plt.ylabel('Von Misses (Pa)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('VM_stress.png', dpi=300,bbox_inches='tight')


#-------------------------------------------------------------------------
#Comparison
#-------------------------------------------------------------------------

#%%**********************deflection***********************************

d1 = np.genfromtxt("../data/h_2_w.txt", comments="%", delimiter="\t")
x_1=np.linspace(a, b, 101)
#x_1=d1[:,0] 
#x_1=x_1+0.1
w_1=d1[:,1]
# plt.figure(figsize=(3.94,3.5), dpi=600)  
plt.figure()
plt.plot(r_vals, w_vals,'r--', label='Analytical_Results')
plt.plot(x_1, w_1,'b--', label='FEM Results')
plt.xlabel('$r$ (m)')
plt.ylabel('$w(r)$ (m)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('Deflection_Analytical_v_FEM.png', dpi=300,bbox_inches='tight')


#%%*********************Normal_Stress****************************

d2 = np.genfromtxt("../data/h_2_s_11.txt", comments="%", delimiter="\t")
x_2=np.linspace(a, b, 101)
#x_2=d2[:,0] 
#x_2=x_2+0.1
s_11=d2[:,1]

# plt.figure(figsize=(3.94,3.1), dpi=600) 
plt.figure()
plt.plot(r_vals, sigma_rr_vals,'r--', label='Analytical_Results')
plt.plot(x_2,s_11,'b--', label='FEM Results')
plt.xlabel('$r$ (m)')
plt.ylabel('$\sigma_{rr}$ (Pa)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('NS_Analytical_v_FEM.png', dpi=300,bbox_inches='tight')



#%%*********************Von_misses_Stress****************************

d3 = np.genfromtxt("../data/h_2_von_misses.txt", comments="%", delimiter="\t")
x_3=np.linspace(a, b, 101)
#x_3=d3[:,0] 
#x_3=x_3+0.1
s_em=d3[:,1]

plt.figure()
plt.plot(r_vals, sigma_em,'m--', label='$\sigma_{eM}$')
plt.plot(x_3, s_em,'b--', label='$\sigma_{eM}$_abaqus')
plt.xlabel('radius')
plt.ylabel('Von Misses (Pa)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('VM_Analytical_v_FEM.png', dpi=300,bbox_inches='tight')



#************************************************************************
#------------------------------------------------------------------------
# %%Reduction of Thickness:
#------------------------------------------------------------------------

B1,B2,B3,B4=symbols('B1,B2,B3,B4',real=True)

h_new = 0.01  # Reduced thickness
D_new = E * h_new**3 / (12 * (1 - v**2))

# New deflection w for reduced thickness
w_new = integrate(1/r * integrate(r * integrate(1/r * integrate(q * r / D_new, r), r), r), r) + B1 * r**2 * log(r/b) + B2 * r**2 + B3 * log(r/b) + B4
wprime_new = diff(w_new, r)

# New moments for reduced thickness
Mr_new = D_new * (-diff(wprime_new, r) - v / r * wprime_new)  
Mphi_new = D_new * (-1 / r * wprime_new - v * diff(wprime_new, r))
V_new = diff(Mr_new, r) + 1 / r * (Mr_new - Mphi_new)

# Solve for constants again
unknowns_n=(B1,B2,B3,B4)
sol_new = solve([w_new.subs(r, b), Mr_new.subs(r, b), Mr_new.subs(r, a), V_new.subs(r, a)], unknowns_n)
w_new_ = simplify(w_new.subs(sol_new))

# Compute deflection for reduced thickness
w_vals_new = [w_new_.subs(r, vals) for vals in r_vals]

# Plot comparison of deflections (original vs reduced thickness)
# plt.figure(figsize=(3.94, 3.5), dpi=600)  
plt.figure()
plt.plot(r_vals, w_vals, 'm--', label='Deflection (h=0.02m)')
plt.plot(r_vals, w_vals_new, 'b--', label='Deflection (h=0.01m)')
plt.xlabel('$r$ (m)')
plt.ylabel('$w(r)$ (m)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('Deflection_original_v_reduced_thickness.png', dpi=300,bbox_inches='tight')


# Compute new Von Mises stresses for reduced thickness
z = h_new / 2
w_double_prime_new=diff(wprime_new,r)

sigma_pp_new = (E/(1-v**2)) * (-v*z*w_double_prime_new.subs(sol_new) -z*wprime_new.subs(sol_new)/r)
sigma_pp_vals_new= []
for val in r_vals:
    sigma_pp_vals_new.append(sigma_pp_new.subs(r, val).evalf())  

sigma_rr_new = (E/(1-v**2)) * (-z*w_double_prime_new.subs(sol_new) -z*v*wprime_new.subs(sol_new)/r)
sigma_rr_vals_new= []
for val in r_vals:
    sigma_rr_vals_new.append(sigma_rr_new.subs(r, val).evalf())  

sigma_em_new=[]
for i in range(len(sigma_pp_vals_new)):
    value = (1.5 * (1/9) * ((2 * sigma_rr_vals_new[i] - sigma_pp_vals_new[i])**2 +
        (2 * sigma_pp_vals_new[i] - sigma_rr_vals_new[i])**2 +
        (sigma_rr_vals_new[i] + sigma_pp_vals_new[i])**2))**0.5
    sigma_em_new.append(value)

# Plot comparison of Von Mises stresses (original vs reduced thickness)
# plt.figure(figsize=(3.94, 3.5), dpi=600)  
plt.figure()
plt.plot(r_vals, sigma_em, 'm--', label='Von Mises (h=0.02m)')
plt.plot(r_vals, sigma_em_new, 'b--', label='Von Mises (h=0.01m)')
plt.xlabel('radius')
plt.ylabel('Von Mises (Pa)')
plt.legend()
plt.grid()
plt.show()
plt.savefig('VM_original_v_reduced_thickness.png', dpi=300,bbox_inches='tight')


# %%
