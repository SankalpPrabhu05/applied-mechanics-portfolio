xi = sym('xi',[2,1],'real');

% Define shape functions
N1 = (1 - xi(1) - xi(2))*(1- 2*xi(1) - 2*xi(2)); 
N2 = xi(1)*(2*xi(1) - 1);
N3 = xi(2)*(2*xi(2) - 1);
N4 = 4*xi(1)*(1 - xi(1) - xi(2));
N5 = 4*xi(1)*xi(2);
N6 = 4*xi(2)*(1 - xi(1) - xi(2));

% Define N-matrix for the element
Ne = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0;
      0 N1 0 N2 0 N3 0 N4 0 N5 0 N6];

% Differentiate shape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,xi);
dN2_dxi=gradient(N2,xi);
dN3_dxi=gradient(N3,xi);
dN4_dxi=gradient(N4,xi);
dN5_dxi=gradient(N5,xi);
dN6_dxi=gradient(N6,xi);

% Introduce initial node positions
Xe1 = sym('xe1',[2,1],'real');
Xe2 = sym('xe2',[2,1],'real');
Xe3 = sym('xe3',[2,1],'real');
Xe4 = sym('xe4',[2,1],'real');
Xe5 = sym('xe5',[2,1],'real');
Xe6 = sym('xe6',[2,1],'real');

% Introduce spatial coordinate as fcn of isoparam. coord.
X = N1*Xe1 + N2*Xe2 + N3*Xe3 + ...
    N4*Xe4 + N5*Xe5 + N6*Xe6;

% Compute Jacobian
Fisop = jacobian(X,xi);
detFisop = det(Fisop);

% Use chain rule to compute material derivatives
dN1_dX=simplify( inv(Fisop)'*dN1_dxi );
dN2_dX=simplify( inv(Fisop)'*dN2_dxi );
dN3_dX=simplify( inv(Fisop)'*dN3_dxi );
dN4_dX=simplify( inv(Fisop)'*dN4_dxi );
dN5_dX=simplify( inv(Fisop)'*dN5_dxi );
dN6_dX=simplify( inv(Fisop)'*dN6_dxi );

% Define B0-matrix of element
Be0 = [dN1_dX(1)        0       dN2_dX(1)       0       dN3_dX(1)       0       dN4_dX(1)       0       dN5_dX(1)       0       dN6_dX(1)       0
        0           dN1_dX(2)      0        dN2_dX(2)       0       dN3_dX(2)       0       dN4_dX(2)       0       dN5_dX(2)       0       dN6_dX(2)
       dN1_dX(2)      0        dN2_dX(2)       0       dN3_dX(2)       0       dN4_dX(2)       0       dN5_dX(2)       0       dN6_dX(2)        0
        0           dN1_dX(1)        0       dN2_dX(1)       0       dN3_dX(1)       0       dN4_dX(1)       0       dN5_dX(1)       0       dN6_dX(1)];

matlabFunction(Be0,detFisop,Ne,'File','Be0_TRIA6_func','Vars',...
{xi,Xe1,Xe2,Xe3,Xe4, Xe5, Xe6});