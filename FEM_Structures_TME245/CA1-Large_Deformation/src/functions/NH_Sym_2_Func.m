%define 2d deformation gradient in Voigt format:
Fv = sym('Fv',[4,1],'real')
F=[Fv(1) Fv(3) 0; Fv(4) Fv(2) 0; 0 0 1];
C=F'*F
invC=simplify(inv(C))
J=det(F)
syms Emod nu 'real'
Gmod = Emod / (2 * (1 + nu));
lambda = (Emod*nu)/((1+nu)*(1-2*nu))
S=Gmod*( eye(3)-invC )+ lambda*log(J)*invC;
P=F*S
Pv=[P(1,1) P(2,2) P(1,2) P(2,1)]'

% Compute Cauchy stress tensor
Sigma = (1 / J) * (F * S * F');

matlabFunction(Pv, Sigma,'File','NH_func','Vars',{Fv,Emod,nu});