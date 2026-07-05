clear all
close all
clc

xi=sym('xi',[2,1],'real');
N1_u=(1-xi(1)-xi(2))*(1-2*xi(1)-2*xi(2));
N3_u=xi(1)*(2*xi(1)-1);
N5_u=xi(2)*(2*xi(2)-1);
N2_u=4*xi(1)*(1-xi(1)-xi(2));
N4_u=4*xi(1)*xi(2);
N6_u=4*xi(2)*(1-xi(1)-xi(2));

dN1_dxi=gradient(N1_u,xi);
dN2_dxi=gradient(N2_u,xi);
dN3_dxi=gradient(N3_u,xi);
dN4_dxi=gradient(N4_u,xi);
dN5_dxi= gradient(N5_u,xi);
dN6_dxi = gradient(N6_u,xi);

N1_p = 1-xi(1)-xi(2);
N2_p = xi(1);
N3_p = xi(2);

Xe1=sym('xe1',[2,1],'real');
Xe2=sym('xe2',[2,1],'real');
Xe3=sym('xe3',[2,1],'real');
Xe4=sym('xe4',[2,1],'real');
Xe5=sym('xe5',[2,1],'real');
Xe6=sym('xe6',[2,1],'real');

X=N1_u*Xe1+N2_u*Xe2+N3_u*Xe3+N4_u*Xe4+N5_u*Xe5+N6_u*Xe6;

Fisop=jacobian(X,xi);
detFisop=det(Fisop);

dN1_dX=simplify(inv(Fisop)'*dN1_dxi);
dN2_dX=simplify(inv(Fisop)'*dN2_dxi);
dN3_dX=simplify(inv(Fisop)'*dN3_dxi);
dN4_dX=simplify(inv(Fisop)'*dN4_dxi);
dN5_dX=simplify(inv(Fisop)'*dN5_dxi);
dN6_dX=simplify(inv(Fisop)'*dN6_dxi);

B_eps = [dN1_dX(1) 0 dN2_dX(1) 0 dN3_dX(1) 0 dN4_dX(1) 0 dN5_dX(1) 0 dN6_dX(1) 0;
      0 dN1_dX(2) 0 dN2_dX(2) 0 dN3_dX(2) 0 dN4_dX(2) 0 dN5_dX(2) 0 dN6_dX(2);
      dN1_dX(2) dN1_dX(1) dN2_dX(2) dN2_dX(1) dN3_dX(2) dN3_dX(1) dN4_dX(2) dN4_dX(1) dN5_dX(2) dN5_dX(1) dN6_dX(2) dN6_dX(1)];

B_div = [dN1_dX(1) dN1_dX(2) dN2_dX(1) dN2_dX(2) dN3_dX(1) dN3_dX(2) dN4_dX(1) dN4_dX(2) dN5_dX(1) dN5_dX(2) dN6_dX(1) dN6_dX(2)];

B_dev = B_eps - 1/3*[1;1;0]*B_div;

N_p = [N1_p N2_p N3_p];

matlabFunction(B_eps,B_div,B_dev,N_p,detFisop,'File','Be_Taylor_Hood_implementation','Vars',{xi,Xe1,Xe2,Xe3,Xe4,Xe5,Xe6})

