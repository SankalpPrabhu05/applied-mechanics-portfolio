clear all
close all
clc

%% Material Properties
L = 40e-3; %m
H = 20e-3; %m
E = 50e9; %Pa
nu = 0.3;
thickness = 1;
mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));
D = hooke(2,E,nu);
Dq = [D(1,1),D(1,2),D(1,4);
    D(2,1),D(2,2),D(2,4);
    D(4,1),D(4,2),D(4,4)];
grid0 = load("grid0.mat");

%% Initiation
num_refinement = 10;
Ndof_plot = zeros(num_refinement,1);
eta = zeros(num_refinement,1); 
current_grid = grid0.grid;
load("exact_energy_norm.mat")

%% Refinement Loop
for j = 1:num_refinement

    %% Linear Element
    p= 1;
    [Edof,Ed_dof,Dof,Ex,Ey,Ed_x,Ed_y,Ed_l,Ed_no,Coord]=grid2CALFEM(current_grid,p);
    nnode = size(Coord,1);
    ndof = 2*nnode;
    Ndof_plot(j)=ndof;

    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    
    for iel = 1:size(Edof,1)

        ex = Ex(iel,:);
        ey = Ey(iel,:);
        edof = Edof(iel,2:end);
        [Ke,~] = plante(ex,ey,[2, thickness],D);
        K(edof,edof)=K(edof,edof)+Ke;
        [Fe] = Force_exact_solution(ex,ey,mu,lambda,L,H,thickness);
        F(edof)=F(edof)+Fe;
    end

    edges = current_grid.edges;
    fixed_dofs = [];

    for edge_id = 1:5

        fix_edge = find(edges(:,3)==edge_id);

        for le = 1:length(fix_edge)

            edge_index = fix_edge(le);
            n1 = edges(edge_index,1);
            n2 = edges(edge_index,2);        
            fixed_dofs = [fixed_dofs, 2*n1-1, 2*n1, 2*n2-1, 2*n2];
        end
    end

    fdof = unique(fixed_dofs);
    free_dof = setdiff(1:ndof,fdof);
    u= zeros(ndof,1);
    K_free = K(free_dof,free_dof);
    F_free = F(free_dof);
    u(free_dof) = K_free\F_free;
    norm_uh_square = 0;

    for il = 1:size(Edof,1)

        ex = Ex(il,:);
        ey = Ey(il,:);
        [Ke,~] = plante(ex,ey,[2, thickness],D);
        ae = u(Edof(il,2:end));
        norm_e = ae'*Ke*ae;
        norm_uh_square = norm_uh_square + norm_e;
    end

    norm_uh = sqrt(norm_uh_square);
    norm_error_square = energy_norm^2-norm_uh^2;
    exact_error = sqrt(abs(norm_error_square));

    %% Quadratic Element
    p = 2;
    [Edofq,Ed_dofq,Dofq,Exq,Eyq,Ed_xq,Ed_yq,Ed_lq,Ed_noq,Coordq]=grid2CALFEM(current_grid,p);
    num_elem_q = size(Edofq,1);
    Ex_quad = zeros(num_elem_q,6);
    Ey_quad = zeros(num_elem_q,6);

    for i = 1:num_elem_q
    
        edof = Edofq(i,2:end);
        xdofs = edof(1:2:end);
        node = (xdofs+1)/2;
        Ex_quad(i,:) = Coordq(node,1);
        Ey_quad(i,:) = Coordq(node,2);
    end

    nnodeq = size(Coordq,1);
    ndofq = 2*nnodeq;
    Kq = sparse(ndofq,ndofq);
    Fq = zeros(ndofq,1);

    for iel = 1:num_elem_q

        ex = Ex_quad(iel,:);
        ey = Ey_quad(iel,:);
        edof = Edofq(iel,2:end);
        [Ke,~] = plant6e(ex,ey,[thickness, 2],Dq);
        Kq(edof,edof)=Kq(edof,edof)+Ke;
        [Fe] = Force_exact_quad_solution(ex,ey,mu,lambda,L,H,thickness);
        Fq(edof)=Fq(edof)+Fe;
    end

    edges_quad = current_grid.edges;
    nlnodes = length(current_grid.nodes(:,1)); 
    fixed_dofs_quad = [];

    for edge_id = 1:5

        fix_edge = find(edges_quad(:,3)==edge_id);

        for le = 1:length(fix_edge)

            edge_index = fix_edge(le);
            n1 = edges_quad(edge_index,1);
            n2 = edges_quad(edge_index,2);
            n_mid = nlnodes+edge_index;
            fixed_dofs_quad = [fixed_dofs_quad, 2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n_mid-1, 2*n_mid];
        end
    end

    fdofq = unique(fixed_dofs_quad);
    free_dof_quad = setdiff(1:ndofq, fdofq);
    uq = zeros(ndofq,1);
    K_free_quad = Kq(free_dof_quad,free_dof_quad);
    F_free_quad = Fq(free_dof_quad,1);
    uq(free_dof_quad) = K_free_quad\F_free_quad;

    %% Approximate error
    a_lin_quad = lin2quad(u,Ed_dofq);
    e_tilde = uq-a_lin_quad;
    norm_e_square = 0;

    for il = 1:size(Edof,1)

        ex = Ex_quad(il,:);
        ey = Ey_quad(il,:);
        [Ke,~] = plant6e(ex,ey,[1,2],Dq);
        ae = e_tilde(Edofq(il,2:end));
        norm_e = ae'*Ke*ae;
        norm_e_square = norm_e_square + norm_e;
    end

    estimated_error = sqrt(norm_e_square);
    eta(j) = estimated_error/exact_error;

    if j <num_refinement

        current_grid = refinegrid(current_grid);
    end
end

%% Plotting

figure
loglog(Ndof_plot,eta,'-o','Color', [1 0 0], 'LineWidth', 1.5, 'DisplayName', 'For Linear Elements')
xlabel('Number of Degree of Freedom');
ylabel('Effective Index');
title('Effective Index v/s Number of Degree of Freedom');
grid on