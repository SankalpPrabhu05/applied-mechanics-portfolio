%% Compute Kirchhoff Plate Element N_oop function
xi = sym('xi',[2,1],'real');

% Define base functions for isoparam quadrilateral element
N1 = 0.25*(xi(1) - 1)*(xi(2) - 1); N2 = -0.25*(xi(1) + 1)*(xi(2) - 1);
N3 = 0.25*(xi(1) +1)*(xi(2) + 1); N4 = -0.25*(xi(1) - 1)*(xi(2) + 1);

% Introduce node positions
xe1 = sym('xe1',[2,1],'real');
xe2 = sym('xe2',[2,1],'real');
xe3 = sym('xe3',[2,1],'real');
xe4 = sym('xe4',[2,1],'real');

% Introduce spatial coordinate as fcn of isoparam. coord.
x=N1*xe1+N2*xe2+N3*xe3+N4*xe4;

% Compute Jacobian
Fisop=jacobian(x,xi);
invFisop=simplify(inv(Fisop));
detFisop = simplify(det(Fisop));

% Create Matlab functions for later use
matlabFunction(invFisop,'File','invFisop_4node_func','Vars',{xi,xe1,xe2,xe3,xe4});
matlabFunction(detFisop,'File','detFisop_4node_func','Vars',{xi,xe1,xe2,xe3,xe4});

% Define P matrix, alpha coefficients and w(xi)
P(xi)=[1 xi(1) xi(2) xi(1)^2 xi(1)*xi(2) xi(2)^2 ... %6
xi(1)^3 xi(1)^2*xi(2) xi(1)*xi(2)^2 xi(2)^3 ...
xi(1)^3*xi(2) xi(1)*xi(2)^3];
alpha = sym('alpha',[12,1],'real');
w(xi)=P(xi(1),xi(2))*alpha;

% Differentiate P with respect to x, y. Use chain rule.
dPdxi = sym('dPxi',[12,2],'real')
dPdxi1(xi)=diff(P(xi(1),xi(2)),xi(1))
dPdxi2(xi)=diff(P(xi(1),xi(2)),xi(2))
dPdx1(xi)=jacobian(P(xi(1),xi(2)),xi)*invFisop(:,1)
dPdx2(xi)=jacobian(P(xi(1),xi(2)),xi)*invFisop(:,2)

% Initialize A matrix
A = sym('A',[12,12],'real');

% Give points (corner nodes in isoparam element)
np=[-1 -1; -1 -1; -1 -1;...
1 -1; 1 -1; 1 -1;...
1 1; 1 1; 1 1; ...
-1 1; -1 1; -1 1];

% Define A matrix based on P and derivatives of P
for i=[1 4 7 10]
    A(i,:)=P(np(i,1),np(i,2));
end

for i=[2 5 8 11]
    A(i,:)=dPdx1(np(i,1),np(i,2));
end

for i=[3 6 9 12]
    A(i,:)=dPdx2(np(i,1),np(i,2));
end

% Finally, solve for N and produce matlab function
N_oop=simplify(P*inv(A));
matlabFunction(N_oop,'File','Nk_func','Vars',{xi,xe1,xe2,xe3,xe4});

%% Compute Be_kirchoff and B* matrix

% Initalize dNdx
dNdx=sym('dNdx',[12,2],'real');

% Define dNdx by using the chain rule
dNdx(:,1) = diff(N_oop(xi(1),xi(2)),xi(1))*invFisop(1,1) + ...
diff(N_oop(xi(1),xi(2)),xi(2))*invFisop(2,1);
dNdx(:,2) = diff(N_oop(xi(1),xi(2)),xi(1))*invFisop(1,2) + ...
diff(N_oop(xi(1),xi(2)),xi(2))*invFisop(2,2);


% Now Bast can be be computed
Bast(1,:) = diff(dNdx(:,1),xi(1))*invFisop(1,1) + ...
diff(dNdx(:,1),xi(2))*invFisop(2,1);
Bast(2,:) = diff(dNdx(:,2),xi(1))*invFisop(1,2) + ...
diff(dNdx(:,2),xi(2))*invFisop(2,2);
Bast(3,:) = 2*(diff(dNdx(:,1),xi(1))*invFisop(1,2) + ...
diff(dNdx(:,1),xi(2))*invFisop(2,2));

matlabFunction(dNdx,'File','Be_kirch_func','Vars',{xi,xe1,xe2,xe3,xe4});
matlabFunction(Bast,'File','Bast_kirchoff_func','Vars',{xi,xe1,xe2,xe3,xe4});