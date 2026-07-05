%define 2d deformation gradient in Voigt format:
Fv = sym('Fv',[4,1],'real')
F=[Fv(1) Fv(3) 0; Fv(4) Fv(2) 0; 0 0 1];
C=F'*F
invC=simplify(inv(C))
J=det(F)
I = eye(3)

syms Emod nu 'real'
Gmod = Emod / (2 * (1 + nu));
c1 = Gmod/2;
c2 = -Gmod/10;
c3 = Gmod/30;
D1 = 0.02;
D2 = 0.01;
D3 = 0.01;

% Compute 2nd Piola-Kirchhoff stress tensor
S = 2 * (c1 + 2 * c2 * (J^(-2/3) * trace(C) - 3) + 3 * c3 * (J^(-2/3) * trace(C) - 3)^2) ...
    * J^(-2/3) * (I - (1/3) * trace(C) * inv(C)) ...
    + 2 * ( (1/D1) + (2/D2) * (J - 1)^2 + (3/D3) * (J - 1)^4 ) * (J - 1) * J * inv(C);

% Compute 1st Piola-Kirchhoff stress tensor
P=F*S
Pv=[P(1,1) P(2,2) P(1,2) P(2,1)]'

% Compute Cauchy stress tensor
Sigma = (1 / J) * (F * S * F');
Sv = [Sigma(1,1) Sigma(2,2) Sigma(1,2) Sigma(2,1) Sigma(3,3)]';

% Compute Material Tangent Stiffness Matrix dP/dF
dPvdFv=sym('dPvFv',[4,4],'real')
for i=1:4
    dPvdFv(i,:)=gradient(Pv(i),Fv)
end


matlabFunction(Pv,dPvdFv,Sv,'File','Yeoh_func_S33','Vars',{Fv,Emod, nu});