clear all
close all
clc

xi = sym('xi',[2,1],'real');

N1=(1-xi(1)-xi(2))*(1-2*xi(1)-2*xi(2));
N3=xi(1)*(2*xi(1)-1);
N5=xi(2)*(2*xi(2)-1);
N2=4*xi(1)*(1-xi(1)-xi(2));
N4=4*xi(1)*xi(2);
N6=4*xi(2)*(1-xi(1)-xi(2));

dN1_dxi = gradient(N1,xi);
dN2_dxi = gradient(N2,xi);
dN3_dxi = gradient(N3,xi);
dN4_dxi = gradient(N4,xi);
dN5_dxi = gradient(N5,xi);
dN6_dxi = gradient(N6,xi);

Xe1=sym('Xe1',[2,1],'real');
Xe2=sym('Xe2',[2,1],'real');
Xe3=sym('Xe3',[2,1],'real');
Xe4=sym('Xe4',[2,1],'real');
Xe5=sym('Xe5',[2,1],'real');
Xe6=sym('Xe6',[2,1],'real');

X=N1*Xe1+N2*Xe2+N3*Xe3+N4*Xe4+N5*Xe5+N6*Xe6;

Fisotr=jacobian(X,xi);

detJ = det(Fisotr);

matlabFunction(X,detJ,'File','Gauss_quad','Vars',{xi,Xe1,Xe2,Xe3,Xe4,Xe5,Xe6});