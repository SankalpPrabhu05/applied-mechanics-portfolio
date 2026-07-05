function [Gr] = Kirch_G(ex, ey, t, eSS)
 % -------------------------------------------------------------
 % PURPOSE
 %   Calculate the Geometric matrix for a Kirchhoff Plate element
 %  
 % 
 %  INPUT:  ex = [x1 x2 x3 x4]      element coordinates
 %          ey = [y1 y2 y3 y4]
 % 
 %          t = scalar              t: thickness
 %
 %          D                       constitutive matrix → Obtained by Hooke
 % 
 %          eSS =  [Sxx;            Element Stress State (output from K_inplane_S function) 
 %                  Syy; 
 %                  Sxy];             
 % 
 %  OUTPUT: Gr : element geometric matrix (12 x 12)

Nsec = t*[eSS(1) eSS(3);
          eSS(3) eSS(2)];

Gr = zeros(12,12);


% Gauss points
ngp = 9;
H_v = [25/81 25/81 25/81 25/81 40/81 40/81 40/81 40/81 64/81];
xi_v = [-sqrt(3/5)  sqrt(3/5) -sqrt(3/5) sqrt(3/5)      0        0      -sqrt(3/5) sqrt(3/5) 0;
        -sqrt(3/5) -sqrt(3/5)  sqrt(3/5) sqrt(3/5) -sqrt(3/5) sqrt(3/5)     0           0    0];

for gp=1:ngp
    Hgp = H_v(gp);
    xin = xi_v(:,gp);

    detFisop = detFisop_4node_func(xin, [ ex(1) ey(1) ]', ...
                                        [ ex(2) ey(2) ]', ...
                                        [ ex(3) ey(3) ]', ...
                                        [ ex(4) ey(4) ]');
    Be_kirch = Be_kirch_func(xin,[ ex(1) ey(1) ]', ...
                                 [ ex(2) ey(2) ]', ...
                                 [ ex(3) ey(3) ]',...
                                 [ ex(4) ey(4) ]');

    Gr = Gr + Be_kirch*Nsec*Be_kirch'*Hgp*detFisop;
end