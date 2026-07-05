function [Kww, Kaa, fwe_ext, fae_ext] = Full_KQuad_Func(ex,ey,ep,D,q,eq)

 % -------------------------------------------------------------
 % PURPOSE
 %   Calculate the stiffness matrix for a Kirchhoff plane stress quad element
 %  
 % 
 %  INPUT:  ex = [x1 x2 x3 x4]      element coordinates
 %          ey = [y1 y2 y3 y4]
 % 
 %          ep = [ptype t]          ptype: analysis type → 1 for Plane Stres
 %                                  t: thickness
 % 
 %          q = scalar              Uniform pressure - out of plane
 %
 %          D                       constitutive matrix → Obtained by Hooke
 % 
 %          eq = [bx;               bx: body force x-dir
 %                by]               by: body force y-dir
 % 
 %  OUTPUT: Kww : element bending stiffness matrix (12 x 12)
 %          Kaa : element stiffness matrix (8 x 8)
 %			fwe_ext: element external force matrix, out-of-plane (12 x 1)
 % 			fae_ext: element external force matrix, in-plane (8 x 1)
 
ptype = ep(1);
h = ep(2);

Dbar = D*h^3/12;

% Gauss points
H_v = ones(1,4);
xi_v = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3);
        -1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];

Kww = zeros(12,12);
Kaa = zeros(8);
fwe_ext = zeros(12,1);
fae_ext = zeros(8,1);

for gp=1:4
    Hgp = H_v(gp);
    xin = xi_v(:,gp);
    detFisop = detFisop_4node_func(xin, [ ex(1) ey(1) ]', ...
                                        [ ex(2) ey(2) ]', ...
                                        [ ex(3) ey(3) ]', ...
                                        [ ex(4) ey(4) ]');
    
    [Be, Ne] = Be_Quad_func(xin, [ ex(1) ey(1) ]', ...
                                 [ ex(2) ey(2) ]', ...
                                 [ ex(3) ey(3) ]',...
                                 [ ex(4) ey(4) ]');
    
    Bastn = Bast_kirchoff_func(xin, [ ex(1) ey(1) ]', ...
                                    [ ex(2) ey(2) ]', ...
                                    [ ex(3) ey(3) ]',...
                                    [ ex(4) ey(4) ]');    

    N_oop = Nk_func(xin,[ ex(1) ey(1) ]',[ ex(2) ey(2) ]',...
                        [ ex(3) ey(3) ]',[ ex(4) ey(4) ]');

    fwe_ext = fwe_ext + N_oop'*q*detFisop*Hgp;
    fae_ext = fae_ext + Ne'*eq*detFisop*h*Hgp;

    Kaa = Kaa + Be'*D*Be*detFisop*h*Hgp;
    Kww = Kww + Bastn'*Dbar*Bastn*detFisop*Hgp;
end