clear all
close all
clc

xi = sym('xi',[2,1],'real');

N1 = 1-xi(1)-xi(2);
N2 = xi(1);
N3 = xi(2);

dN1_dxi = gradient(N1,xi);
dN2_dxi = gradient(N2,xi);
dN3_dxi = gradient(N3,xi);

Xe1=sym('Xe1',[2,1],'real');
Xe2=sym('Xe2',[2,1],'real');
Xe3=sym('Xe3',[2,1],'real');

X=N1*Xe1+N2*Xe2+N3*Xe3;

Fisotr=jacobian(X,xi);

detJ = det(Fisotr);

matlabFunction(X,detJ,'File','Gauss','Vars',{xi,Xe1,Xe2,Xe3});