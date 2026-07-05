clear all
close all
clc

%%

ny_list = [2,4,8,10,15,20,25,30];
nu0 = [0.3,0.4,0.45,0.48,0.49,0.5];

E = 200e9; %   Pa
t = 1; % m
H = 20e-3; %m

%% Case 1
Tval = 1e8; % Pa
ueval = zeros(size(nu0,2),size(ny_list,2));
ueval_c = zeros(size(nu0,2),1);
% Case 2

Uval = 1e-7; %m
Teval = zeros(size(nu0,2),size(ny_list,2));
Teval_c = zeros(size(nu0,2),1);
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

%% Problem 3 

for iny = 1:size(ny_list,2)
    ny = ny_list(iny);
    order = 2;
    [coords, edges, triangle] = mesh(ny,order,false);
    
    Ex = coords(1,:);
    Ey = coords(2,:);
    
    nnodes_u = size(coords,2);
    vertices = [triangle(1,:),triangle(3,:),triangle(5,:)];
    unique_vertices = unique(vertices);
    nnodes_p = size(unique_vertices,2);
    
    ndofs = nnodes_u*2+nnodes_p;
    nelem = size(triangle,2);
    
    bound_id = 4;
    
    %% Case I
    
    for inu=1:size(nu0,2)
    
        nu = nu0(inu);
        K = sparse(ndofs,ndofs);
        F = zeros(ndofs,1);
        
        for ielem = 1:nelem
            unodes = triangle(:,ielem);
            pnodes = triangle(1:2:end,ielem);
        
            ex = coords(1,unodes);
            ey = coords(2,unodes);
            [Ke]= Taylor_hood(ex,ey,E,nu,t);
        
            edof = zeros(1,15);
        
            for k = 1:6
                n = unodes(k);
                edof(2*k-1:2*k) = [2*n-1,2*n];
            end
        
            for k= 1:3
                n=pnodes(k);
                np= find(unique_vertices==n)+nnodes_u*2;
                edof(2*size(triangle,1)+k)= np;
            end
            Edof(ielem,:) = edof;
            K(edof,edof)=K(edof,edof)+Ke;
        end
    
        [funit,gamma_L] = Loadvector_interpolation_Taylor_hood(coords,edges,bound_id,ndofs);
    
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
        free_dof = setdiff(1:ndofs,fdof);
    
        u= zeros(ndofs,1);
    
        K_free = K(free_dof,free_dof);
        F_free = F_T(free_dof);
    
        u(free_dof) = K_free\F_free;
    
        ueval(inu,iny) = -(1/(gamma_L))*(funit'*u);
    end
        
    
    %% Case II
    
    for inu=1:size(nu0,2)
    
        nu = nu0(inu);
        K = sparse(ndofs,ndofs);
        F = zeros(ndofs,1);
    
        for ielem = 1:nelem
            unodes = triangle(:,ielem);
            pnodes = triangle(1:2:end,ielem);
    
            ex = coords(1,unodes);
            ey = coords(2,unodes);
            [Ke]= Taylor_hood(ex,ey,E,nu,t);
    
            edof = zeros(1,15);
    
            for k = 1:6
                n = unodes(k);
                edof(2*k-1:2*k) = [2*n-1,2*n];
            end
    
            for k= 1:3
                n=pnodes(k);
                np= find(unique_vertices==n)+nnodes_u*2;
                edof(2*size(triangle,1)+k)= np;
            end
    
            K(edof,edof)=K(edof,edof)+Ke;
        end
    
        [funit,gamma_L] = Loadvector_interpolation_Taylor_hood(coords,edges,bound_id,ndofs);
    
        
    
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
        free_dof2 = setdiff(1:ndofs,fixdof2);
        
        u2= zeros(ndofs,1);
    
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
        'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('\\nu = %.3f', nu0(j)));
end

xlabel('Number of Elements in Y-direction');
ylabel('Average Displacement (m)');
title('Mesh Convergence for Uniform Traction Load Case');
legend(Location='best');

figure; hold on; grid on;
colors = lines(numel(nu0));

for j = 1:numel(nu0)
    plot(ny_list, Teval(j,:), '-o', ...
        'Color', colors(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('\\nu = %.3f', nu0(j)));
end

xlabel('Number of Elements in Y-direction');
ylabel('Average Traction (Pa)');
title('Mesh Convergence for Pres. Displacement Load Case');
legend(Location='best');


%% Pressure Distribution

ny = 15;
order = 2;
[coords, edges, triangle] = mesh(ny,order,false);

Ex = coords(1,:);
Ey = coords(2,:);

nnodes_u = size(coords,2);
vertices = [triangle(1,:),triangle(3,:),triangle(5,:)];
unique_vertices = unique(vertices);
nnodes_p = size(unique_vertices,2);

ndofs = nnodes_u*2+nnodes_p;
nelem = size(triangle,2);

bound_id = 4;

%% Case I

for inu=1:size(nu0,2)

    nu = nu0(inu);
    K = sparse(ndofs,ndofs);
    F = zeros(ndofs,1);
    
    for ielem = 1:nelem
        unodes = triangle(:,ielem);
        pnodes = triangle(1:2:end,ielem);
    
        ex = coords(1,unodes);
        ey = coords(2,unodes);
        [Ke]= Taylor_hood(ex,ey,E,nu,t);
    
        edof = zeros(1,15);
    
        for k = 1:6
            n = unodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end
    
        for k= 1:3
            n=pnodes(k);
            np= find(unique_vertices==n)+nnodes_u*2;
            edof(2*size(triangle,1)+k)= np;
        end

        K(edof,edof)=K(edof,edof)+Ke;
    end

    [funit,gamma_L] = Loadvector_interpolation_Taylor_hood(coords,edges,bound_id,ndofs);

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
    free_dof = setdiff(1:ndofs,fdof);

    u= zeros(ndofs,1);

    K_free = K(free_dof,free_dof);
    F_free = F_T(free_dof);

    u(free_dof) = K_free\F_free;

    ueval_c(inu,1) = -(1/(gamma_L))*(funit'*u);

    p1 = zeros(nelem,1);

    for ip = 1:nelem

        pnodes = triangle(1:2:end,ip);
        press_on_nodes_1 = zeros(numel(pnodes),1);
        for k = 1:numel(press_on_nodes_1)
            idx= find(unique_vertices==pnodes(k));
            press_on_nodes_1(k) = u(nnodes_u*2+idx);
        end
        p1(ip,1)= mean(press_on_nodes_1);
    end
    
    figure; hold on; axis equal; view(2);
    patch('Faces', triangle(1:2:end,:)', 'Vertices', coords', ...
      'FaceColor', 'flat', 'CData', p1, 'EdgeColor', 'none');
    colormap(jet); colorbar;
    title(['Pressure for Uniform Traction Load at \nu=',num2str(nu)]);

end

figure
semilogy(mu,U_1,LineWidth=1.5,DisplayName='\epsilon_{11}=0')
hold on
semilogy(mu,U_2,LineWidth=1.5,DisplayName='\sigma_{11}=0')
semilogy(nu0,ueval_c, LineWidth=1.5,DisplayName='Uniform Traction load')
hold off

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

%% Case II

for inu=1:size(nu0,2)

    nu = nu0(inu);
    K = sparse(ndofs,ndofs);
    F = zeros(ndofs,1);

    for ielem = 1:nelem
        unodes = triangle(:,ielem);
        pnodes = triangle(1:2:end,ielem);

        ex = coords(1,unodes);
        ey = coords(2,unodes);
        [Ke]= Taylor_hood(ex,ey,E,nu,t);

        edof = zeros(1,15);

        for k = 1:6
            n = unodes(k);
            edof(2*k-1:2*k) = [2*n-1,2*n];
        end

        for k= 1:3
            n=pnodes(k);
            np= find(unique_vertices==n)+nnodes_u*2;
            edof(2*size(triangle,1)+k)= np;
        end

        K(edof,edof)=K(edof,edof)+Ke;
    end

    [funit,gamma_L] = Loadvector_interpolation_Taylor_hood(coords,edges,bound_id,ndofs);


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
    free_dof2 = setdiff(1:ndofs,fixdof2);
    
    u2= zeros(ndofs,1);

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
    % p2 = zeros(nnodes_u,1);
    % p2(unique_vertices) = u2(nnodes_u*2+1:end);

    for ip = 1:nelem

        pnodes = triangle(1:2:end,ip);
        press_on_nodes_2 = zeros(numel(pnodes),1);
        for k = 1:numel(press_on_nodes_2)
            idx= find(unique_vertices==pnodes(k));
            press_on_nodes_2(k) = u2(nnodes_u*2+idx);
        end
        p2(ip,1)= mean(press_on_nodes_2);
    end
    
    figure; hold on; axis equal; view(2);
    patch('Faces', triangle(1:2:end,:)', 'Vertices', coords', ...
      'FaceColor', 'flat', 'CData', p2, 'EdgeColor', 'none');
    colormap(jet); colorbar;
    title(['Pressure for Prescribed Displacement Load at \nu=',num2str(nu)]);

end

figure
semilogy(mu,T_1,LineWidth=1.5,DisplayName='\epsilon_{11}=0')
hold on
semilogy(mu,T_2,LineWidth=1.5,DisplayName='\sigma_{11}=0')
semilogy(nu0,Teval_c, LineWidth=1.5,DisplayName='Prescribed Displacement load')
hold off

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