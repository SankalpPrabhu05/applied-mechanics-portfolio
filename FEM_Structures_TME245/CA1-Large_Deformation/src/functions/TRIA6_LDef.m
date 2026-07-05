function [fe_int,Ke_int,Fv_m,Pv_m,Sv_m, dPv_dFv_m]=...
TRIA6_LDef(ae, ex, ey, mtrlpar, t)

% TRIA6_LDef Element Routine for 2nd order triangular elements under large deformation
%
% INPUT:  ex = [x1 x2 ... x6]
%         ey = [y1 y2 ... y6]                   Element node coordinates
%         ae = [u1 v1 u2 v2 ... u6 v6]          Nodal displacements u and v     
%         mtrlpar =  matlab structure containing Emod and nu, used as input the Yeoh_func
% 
% OUTPUT: Ke_int : element internal stiffness matrix (12 x 12)
%         fe_int : element internal forces vector (12 x 1)
%         Fv_m   : deformation gradient at each Gauss Points
%         Pv_m   : 1st Piola-Kirchhoff stress at each Gauss Point
%         Sv_m   : Cauchy stress at each Gauss Point
%         dPv_dFv_m : material tangent matrix

E = mtrlpar.Emod;
nu = mtrlpar.nu;
Xe = zeros(2,size(ex,2));
for i=1:size(ex,2)
    Xe(:,i) = [ex(i) ey(i)]';
end

% Define Gauss Integration points and weights (same for all elements)
ngp = 3;
H_v = ones(1,ngp) / 6;         % Equal weight of 1/6 for each gauss point
xi_v = [1/6 1/6 2/3;           % Each column represents the position of a Gauss Point
        1/6 2/3 1/6];         

% Initialize output
Fv_m = zeros(4,ngp);
Pv_m = zeros(4,ngp);
Sv_m = zeros(4,ngp);
dPv_dFv_m = zeros(4,4,ngp);

% Initialize internal element force and stiffness
Ke_int = zeros(12);
fe_int = zeros(12,1);
I_v = [1 1 0 0]';       % Kroenecker's Delta in Voigt's Form

% Computation of deformation gradient, stresses, internal element stiffness and force
for j=1:ngp
    H = H_v(j);
    xi = xi_v(:,j);
    [Be0, detFisop, ~] = Be0_TRIA6_func(xi,Xe(:,1),Xe(:,2),Xe(:,3),Xe(:,4),Xe(:,5),Xe(:,6));
    Fv_m(:,j) = Be0*ae + I_v;    % Deformation gradient at Gauss points
    
    % Piola-K stress and material tangent stiffness from Yeoh Near-Incompressible model
    [Pv_m(:,j), dPv_dFv_m(:,:,j), Sv_m(:,j)] = Yeoh_func(Fv_m(:,j), E, nu);

    Ke_int = Ke_int + Be0' * dPv_dFv_m(:,:,j) * Be0 * detFisop*t*H;               % Internal element stiffness
    fe_int = fe_int + Be0' * Pv_m(:,j) * detFisop*t*H;                            % Internal element force
end

end