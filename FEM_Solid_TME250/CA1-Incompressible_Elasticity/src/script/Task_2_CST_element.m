clear all
close all
clc
%%

ny_list = [1,3,5,7,10,15,20,25];
nu0 = [0.3,0.4,0.45,0.48,0.49,0.499,0.4999];

E = 200e9; %   Pa
thickness = 1; % m
H = 20e-3; %m

%% Load Cases

%Case 1
Tval = 1e8; % Pa
ueval = zeros(size(nu0,2),size(ny_list,2)); 
%Case 2

Uval = 1e-7; %m
Teval = zeros(size(nu0,2),size(ny_list,2));
%% Problem 1


mu = linspace(0.3,0.4999,100);
U_1 = zeros(size(mu,2),1);
U_2 = zeros(size(mu,2),1);
T_1 = zeros(size(mu,2),1);
T_2 = zeros(size(mu,2),1);

for i = 1:size(mu,2)

    U_1(i,1) = Tval/((E*(1-mu(i)))/((1+mu(i))*(1-2*mu(i))*H));
    U_2(i,1) = Tval/((E)/((1-mu(i)^2)*H));
    T_1(i,1) = Uval*((E*(1-mu(i)))/((1+mu(i))*(1-2*mu(i))*H));
    T_2(i,1) = Uval*((E)/((1-mu(i)^2)*H));
end

%% Problem 2

for iny = 1:size(ny_list,2)

ny = ny_list(iny);
order = 1;
[coords, edges, triangle] = mesh(ny,order,false);

Ex = coords(1,:);
Ey = coords(2,:);


nnodes = size(coords,2);
nelem = size(triangle,2);
Edof = zeros(nelem, size(triangle,1)*2+1);

for i = 1:nelem
    
    dofx = 2*triangle(:,i)-1;
    dofy = 2*triangle(:,i);
    Edof(i,1)=i;
    Edof(i,2:2:end) = dofx;
    Edof(i,3:2:end) = dofy;

end
%% Material Properties


ptype = 2; 
ep = [ptype thickness];

ndof = 2*nnodes;

%% Case I - Problem 2

for inu = 1:size(nu0,2)
    
    nu = nu0(inu);

    D = hooke(ptype,E,nu);

    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    bound_id = 4;
    
    for ielem = 1:nelem
        
        enodes = triangle(:,ielem);
        ex = coords(1,enodes);
        ey = coords(2,enodes);
        [Ke,fe] = plante(ex,ey,ep,D);

        edof = zeros(1,6);

        for k = 1:3
            
            n = enodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end

        K(edof,edof) = K(edof,edof)+Ke;
        F(edof) = F(edof)+fe;
        
    end

    [funit,gamma_L] = Loadvector_interpolation(coords,edges,bound_id);

    F_T = F + (-Tval)*funit;
    
    fixed_dofs = [];
    left_edge = find(edges(end,:)==5);

    for le = 1:size(left_edge,2)
        n1 = edges(1,left_edge(1,le));
        n2 = edges(2,left_edge(1,le));

        fixed_dofs = [fixed_dofs,2*n1-1,2*n2-1];
    end
    
    bottom_edge = find(edges(end,:)==1);

    for be = 1:size(bottom_edge,2)
        n1 = edges(1,bottom_edge(1,be));
        n2 = edges(2,bottom_edge(1,be));

        fixed_dofs = [fixed_dofs,2*n1,2*n2];
    end
    
    fdof = unique(fixed_dofs);
    free_dof = setdiff(1:ndof,fdof);

    u= zeros(ndof,1);

    K_free = K(free_dof,free_dof);
    F_free = F_T(free_dof);

    u(free_dof) = K_free\F_free;
    
    ueval(inu,iny) = -(1/(gamma_L))*(funit'*u);
    % 
    
end


%% Case II - Problem 2

for inu = 1:size(nu0,2)

    nu = nu0(inu);

    D = hooke(ptype,E,nu);

    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    bound_id = 4;

    for ielem = 1:nelem
        enodes = triangle(:,ielem);
        ex = coords(1,enodes);
        ey = coords(2,enodes);
        [Ke,fe] = plante(ex,ey,ep,D);

        edof = zeros(1,6);

        for k = 1:3

            n = enodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end

        K(edof,edof) = K(edof,edof)+Ke;
        F(edof) = F(edof)+fe;

    end

    fixed_dofs2 = [];
    left_edge = find(edges(end,:)==5);

    for le = 1:size(left_edge,2)
        n1 = edges(1,left_edge(1,le));
        n2 = edges(2,left_edge(1,le));

        fixed_dofs2 = [fixed_dofs2,2*n1-1,2*n2-1];
    end
    
    bottom_edge = find(edges(end,:)==1);

    for be = 1:size(bottom_edge,2)
        n1 = edges(1,bottom_edge(1,be));
        n2 = edges(2,bottom_edge(1,be));

        fixed_dofs2 = [fixed_dofs2,2*n1,2*n2];
    end
    
    top_left_edge = find(edges(end,:)==4);
    pres_dofs = [];
    for tle = 1:size(top_left_edge,2)
        n1 = edges(1,top_left_edge(1,tle));
        n2 = edges(2,top_left_edge(1,tle));

        fixed_dofs2 = [fixed_dofs2,2*n1,2*n2];
        pres_dofs = [pres_dofs, 2*n1, 2*n2];
    end

    fixdof2 = unique(fixed_dofs2);
    pdof = unique(pres_dofs);
    free_dof2 = setdiff(1:ndof,fixdof2);

    u2 = zeros(ndof,1);
    u2(pdof)= -Uval;
    K_f = K(free_dof2,free_dof2);
    F_f = F(free_dof2,1);
    K_fc = K(free_dof2,fixdof2);

    u_f = K_f\(F_f-K_fc*u2(fixdof2));

    u2(free_dof2)=u_f;

    RF = K*u2;

    T = 0;

    for i = 1: size(pdof,2)
        T = T+ RF(pdof(i));
    end

    Teval(inu,iny) = -T/gamma_L ; 
   
end


end

figure; hold on; grid on;
colors = lines(numel(nu0));

for j = 1:numel(nu0)
    plot(ny_list, ueval(j,:), '-o', ...
        'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('\\nu = %.4f', nu0(j)));
end
set(gca, 'FontSize');
xlabel('Number of Elements in Y-direction');
ylabel('Average Displacement');
title('Mesh Convergence for Uniform Traction Load');
legend(Location='best');

figure; hold on; grid on;
colors = lines(numel(nu0));

for j = 1:numel(nu0)
    plot(ny_list, Teval(j,:), '-o', ...
        'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('\\nu = %.4f', nu0(j)));
end
set(gca, 'FontSize');
xlabel('Number of Elements in Y-direction');
ylabel('Average Traction');
title('Mesh Convergence for Pres. Displacement Load');
legend(Location='best');






%% Pressure Distribution

ny = 20;
order = 1;
[coords, edges, triangle] = mesh(ny,order,false);

Ex = coords(1,:);
Ey = coords(2,:);

ueval_c = zeros(size(nu0,2),1);
Teval_c = zeros(size(nu0,2),1);

nnodes = size(coords,2);
nelem = size(triangle,2);
Edof = zeros(nelem, size(triangle,1)*2+1);

for i = 1:nelem
    
    dofx = 2*triangle(:,i)-1;
    dofy = 2*triangle(:,i);
    Edof(i,1)=i;
    Edof(i,2:2:end) = dofx;
    Edof(i,3:2:end) = dofy;

end
%% Material Properties


ptype = 2; 
ep = [ptype thickness];

ndof = 2*nnodes;

%% Case I - Problem 2

for inu = 1:size(nu0,2)
    
    nu = nu0(inu);

    D = hooke(ptype,E,nu);

    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    bound_id = 4;
    
    for ielem = 1:nelem
        
        enodes = triangle(:,ielem);
        ex = coords(1,enodes);
        ey = coords(2,enodes);
        [Ke,fe] = plante(ex,ey,ep,D);

        edof = zeros(1,6);

        for k = 1:3
            
            n = enodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end

        K(edof,edof) = K(edof,edof)+Ke;
        F(edof) = F(edof)+fe;
        
    end

    [funit,gamma_L] = Loadvector_interpolation(coords,edges,bound_id);

    F_T = F + (-Tval)*funit;
    
    fixed_dofs = [];
    left_edge = find(edges(end,:)==5);

    for le = 1:size(left_edge,2)
        n1 = edges(1,left_edge(1,le));
        n2 = edges(2,left_edge(1,le));

        fixed_dofs = [fixed_dofs,2*n1-1,2*n2-1];
    end
    
    bottom_edge = find(edges(end,:)==1);

    for be = 1:size(bottom_edge,2)
        n1 = edges(1,bottom_edge(1,be));
        n2 = edges(2,bottom_edge(1,be));

        fixed_dofs = [fixed_dofs,2*n1,2*n2];
    end
    
    fdof = unique(fixed_dofs);
    free_dof = setdiff(1:ndof,fdof);

    u= zeros(ndof,1);

    K_free = K(free_dof,free_dof);
    F_free = F_T(free_dof);

    u(free_dof) = K_free\F_free;
    
    ueval_c(inu,1) = -(1/(gamma_L))*(funit'*u);
    % 
    p1 = zeros(nelem,1);

    for ie = 1:nelem

        enodes = triangle(:,ie);
        ex1 = coords (1,enodes);
        ey1 = coords (2,enodes);

        ed = zeros(6,1);

        for k = 1:size(ed,1)/2

            n = enodes(k);
            ed(2*k-1)= u(2*n-1);
            ed(2*k)= u(2*n);

        end

        [es,et] = plants(ex1,ey1,ep,D,ed');
        p1(ie) = (-1/3)*(es(1)+es(2));

    end
    

    figure; hold on; axis equal; view(2);
    patch('Faces', triangle', 'Vertices', coords', ...
      'FaceColor', 'flat', 'CData', p1, 'EdgeColor', 'none');
    colormap(jet); colorbar;
    title(['Pressure for Uniform traction Load at \nu=',num2str(nu)]);

end



figure
semilogy(mu,U_1,LineWidth=1.5,DisplayName='\epsilon_{11}=0')
hold on
semilogy(mu,U_2,LineWidth=1.5,DisplayName='\sigma_{11}=0')
semilogy(nu0,ueval_c, LineWidth=1.5,DisplayName='Uniform Traction load')
hold off
set(gca, 'FontSize');
xlabel ('Poisson Ratio')
ylabel ('Displacement (m)')
title('Average Displacement For CASE I')
legend (Location="best")

figure

semilogy(nu0,ueval_c, LineWidth=1.5,DisplayName='Uniform Traction load')
set(gca, 'FontSize');
xlabel ('Poisson Ratio')
ylabel ('Displacement (m)')
title('Average Displacement For CASE I')
legend (Location="best")

%% Case II - Problem 2

for inu = 1:size(nu0,2)

    nu = nu0(inu);

    D = hooke(ptype,E,nu);

    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    bound_id = 4;

    for ielem = 1:nelem
        enodes = triangle(:,ielem);
        ex = coords(1,enodes);
        ey = coords(2,enodes);
        [Ke,fe] = plante(ex,ey,ep,D);

        edof = zeros(1,6);

        for k = 1:3

            n = enodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end

        K(edof,edof) = K(edof,edof)+Ke;
        F(edof) = F(edof)+fe;

    end

    fixed_dofs2 = [];
    left_edge = find(edges(end,:)==5);

    for le = 1:size(left_edge,2)
        n1 = edges(1,left_edge(1,le));
        n2 = edges(2,left_edge(1,le));

        fixed_dofs2 = [fixed_dofs2,2*n1-1,2*n2-1];
    end
    
    bottom_edge = find(edges(end,:)==1);

    for be = 1:size(bottom_edge,2)
        n1 = edges(1,bottom_edge(1,be));
        n2 = edges(2,bottom_edge(1,be));

        fixed_dofs2 = [fixed_dofs2,2*n1,2*n2];
    end
    
    top_left_edge = find(edges(end,:)==4);
    pres_dofs = [];
    for tle = 1:size(top_left_edge,2)
        n1 = edges(1,top_left_edge(1,tle));
        n2 = edges(2,top_left_edge(1,tle));

        fixed_dofs2 = [fixed_dofs2,2*n1,2*n2];
        pres_dofs = [pres_dofs, 2*n1, 2*n2];
    end

    fixdof2 = unique(fixed_dofs2);
    pdof = unique(pres_dofs);
    free_dof2 = setdiff(1:ndof,fixdof2);

    u2 = zeros(ndof,1);
    u2(pdof)= -Uval;
    K_f = K(free_dof2,free_dof2);
    F_f = F(free_dof2,1);
    K_fc = K(free_dof2,fixdof2);

    u_f = K_f\(F_f-K_fc*u2(fixdof2));

    u2(free_dof2)=u_f;

    RF = K*u2;

    T = 0;

    for i = 1: size(pdof,2)
        T = T+ RF(pdof(i));
    end

    Teval_c(inu,1) = -T/gamma_L ; 

    p2 = zeros(nelem,1);

    for ie = 1:nelem

        enodes = triangle(:,ie);
        ex2 = coords (1,enodes);
        ey2 = coords (2,enodes);

        ed = zeros(6,1);

        for k = 1:size(ed,1)/2

            n = enodes(k);
            ed(2*k-1)= u2(2*n-1);
            ed(2*k)= u2(2*n);

        end

        [es,et] = plants(ex2,ey2,ep,D,ed');
        p2(ie) = (-1/3)*(es(1)+es(2));

    end

    figure; hold on; axis equal; view(2);
    patch('Faces', triangle', 'Vertices', coords', ...
      'FaceColor', 'flat', 'CData', p2, 'EdgeColor', 'none');
    colormap(jet); colorbar;
    title(['Pressure for Prescribed Displacement Load at \nu=',num2str(nu)]);

end

%%
figure

semilogy(mu,T_1,LineWidth=1.5,DisplayName='\epsilon_{11}=0')
hold on
semilogy(mu,T_2,LineWidth=1.5,DisplayName='\sigma_{11}=0')
semilogy(nu0,Teval_c, LineWidth=1.5,DisplayName='Prescribed Displacement load')
hold off
set(gca, 'FontSize');
xlabel ('Poisson Ratio')
ylabel ('Traction (Pa)')
title('Average Traction For CASE II')
legend (Location="best")

figure

semilogy(nu0,Teval_c, LineWidth=1.5,DisplayName='Prescribed Displacement Load')
set(gca, 'FontSize');
xlabel ('Poisson Ratio')
ylabel ('Traction (Pa)')
title('Average Traction For CASE II')
legend (Location="best")