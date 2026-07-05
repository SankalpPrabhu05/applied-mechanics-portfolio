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
num_refinement = 12;
grid0 = load("grid0.mat");

%% Relative Error: Linear Elements
Ndof_plot = zeros(num_refinement,1);
rel_error = zeros(num_refinement,1);
load("exact_energy_norm.mat");
current_grid = grid0.grid;

for j = 1:num_refinement

    p = 1;
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
    norm_error = sqrt(abs(norm_error_square));
    relative_error = norm_error/energy_norm;
    rel_error(j)=relative_error;

    if j <num_refinement

        current_grid = refinegrid(current_grid);
    end
end

fprintf('\n========== SUMMARY ==========\n');
fprintf('Level | N_dof | ||u||_a | ||u_h||_a | Rel Error\n');
fprintf('------|-------|---------|-----------|----------\n');

for j = 1:num_refinement

    uh_computed = sqrt(abs(energy_norm^2 - (rel_error(j)*energy_norm)^2));
    fprintf('%5d | %5d | %7.2f | %9.2f | %.6f\n', ...
        j, Ndof_plot(j), energy_norm, uh_computed, rel_error(j));
end


%% Relative Error: Quadratic Elements
Ndof_plot_quad = zeros(num_refinement,1);
rel_error_quad = zeros(num_refinement,1);
energy_norm_quad = zeros(num_refinement,1);
current_grid_quad = grid0.grid;

for j = 1:num_refinement
    
    p = 2;   
    [Edofq,Ed_dofq,Dofq,Exq,Eyq,Ed_xq,Ed_yq,Ed_lq,Ed_noq,Coordq]=grid2CALFEM(current_grid_quad,p);
    num_elem_q = size(Edofq,1);
    Ex_quad = zeros(num_elem_q,6);
    Ey_quad = zeros(num_elem_q,6);
    nlnodes = length(current_grid_quad.nodes(:,1));
    
    for i = 1:num_elem_q
    
        edof = Edofq(i,2:end);
        xdofs = edof(1:2:end);
        node = (xdofs+1)/2;
        Ex_quad(i,:) = Coordq(node,1);
        Ey_quad(i,:) = Coordq(node,2);
    end
    
    nnodeq = size(Coordq,1);
    ndofq = 2*nnodeq;
    Ndof_plot_quad(j) = ndofq;

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

    edges_quad = current_grid_quad.edges;
    nlnodes = length(current_grid_quad.nodes(:,1)); 
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
    norm_uh_square_quad = 0;

    for il = 1:num_elem_q

        ex = Ex_quad(il,:);
        ey = Ey_quad(il,:);
        [Ke,~] = plant6e(ex,ey,[1, 2],Dq);
        ae = uq(Edofq(il,2:end));
        norm_e = ae'*Ke*ae;
        norm_uh_square_quad = norm_uh_square_quad + norm_e;
    end

    norm_uh_quad = sqrt(norm_uh_square_quad);
    norm_error_square_quad = energy_norm^2-norm_uh_quad^2;
    norm_error_quad = sqrt(abs(norm_error_square_quad));
    relative_error_quad = norm_error_quad/energy_norm;
    rel_error_quad(j)=relative_error_quad;

    if j <num_refinement

        current_grid_quad = refinegrid(current_grid_quad);
    end
end

fprintf('\n========== SUMMARY ==========\n');
fprintf('Level | N_dof | ||u||_a | ||u_h||_a | Rel Error\n');
fprintf('------|-------|---------|-----------|----------\n');

for j = 1:num_refinement

    uh_computed = sqrt(abs(energy_norm^2 - (rel_error_quad(j)*energy_norm)^2));
    fprintf('%5d | %5d | %7.2f | %9.2f | %.6f\n', ...
        j, Ndof_plot_quad(j), energy_norm, uh_computed, rel_error_quad(j));
end

%% Plotting
figure
loglog(Ndof_plot,rel_error,'-o','Color', [0 0 1], 'LineWidth', 1.5, 'DisplayName', 'For Linear Elements')
hold on
loglog(Ndof_plot_quad,rel_error_quad,'-o','Color',[1 0 0] , 'LineWidth', 1.5, 'DisplayName', 'For Quadratic Elements')
xlabel('Number of Degree of Freedom');
ylabel('Relative Error');
title('Relative Error v/s Number of Degree of Freedom');
grid on
legend(Location='best');

p_fit = polyfit(log(Ndof_plot_quad), log(rel_error_quad), 1);
slope = p_fit(1);
q_quad = -2 * slope;

p_fit = polyfit(log(Ndof_plot), log(rel_error), 1);
slope = p_fit(1);
q = -2 * slope;

fprintf('Convergence rate q in linear element : %.3f\n', q);
fprintf('Convergence rate q in quadratic element : %.3f\n', q_quad);
