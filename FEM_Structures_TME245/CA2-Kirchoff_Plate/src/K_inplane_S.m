function [e_ip_sig] = K_inplane_S(aue, awe, ex, ey, z, D)
% Calculates In-Plane Stress for Kirchhoff plate element in Voigt Form

% Gauss points
ngp = 4;
xi_v = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3);
        -1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];

Sig_ip_e = zeros(3,ngp); %       GP1 GP2 GP3 GP4        - Sigma_inplane_element
                         % S11
                         % S22
                         % S12

for j=1:ngp
    xin = xi_v(:,j);               % Gauss-Point Coordinate

    [Be, ~] =  Be_Quad_func(xin, [ ex(1) ey(1) ]', ...
                                 [ ex(2) ey(2) ]', ...
                                 [ ex(3) ey(3) ]',...
                                 [ ex(4) ey(4) ]');
    
    Bastn = Bast_kirchoff_func(xin, [ ex(1) ey(1) ]', ...
                                    [ ex(2) ey(2) ]', ...
                                    [ ex(3) ey(3) ]',...
                                    [ ex(4) ey(4) ]');

    Sig_ip_e(:,j) = D * (Be*aue' - z*Bastn*awe');
end

e_ip_sig = mean(Sig_ip_e,2);        % Averages Gauss Point stresses - element_inplane_sigma (averaged out) 3x1 vector