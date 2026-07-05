clear all
close all
clc

%% Material Properties

E = 50e3; %Mpa
nu = 0.3;
rho = 8e-6;%kg/mm3
T = 5e-3;%sec
dt = 20e-6;%sec
v_ini = 10e3;%mm/sec
edge_length = 2;%mm
thickness = 1; % Assumed
R = 40;
H = 20;%mm
L = 2*H;%mm

%% Meshing
ny = H/edge_length;
order = 1;
[coords, edges, triangle] = mesh(ny,order,false);
coords = 1000.*coords;%mm
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

%% Element Routine

ptype = 2; 
ep = [ptype thickness];
ndof = 2*nnodes;
D = hooke(ptype,E,nu);
S = sparse(ndof,ndof);
F = zeros(ndof,1);

for ielem = 1:size(Edof,1)
    
    enodes = triangle(:,ielem);
    ex = Ex(enodes);
    ey = Ey(enodes);
    [Ke,fe] = plante(ex,ey,ep,D);
    edof = Edof(ielem,2:end);
    S(edof,edof) = S(edof,edof)+Ke;
    F(edof) = F(edof)+fe;

end

%% Boundary Conditions

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

fixdof = unique(fixed_dofs);

%% Cylinder movement

y_ini = 60;
y_final = 55;
n_s = 250; %Number of steps
y = linspace(y_ini,y_final,n_s);
X_bar = zeros(2,n_s);
X_bar(1,:) = 0;
X_bar(2,:) = y;

%% Defining Contact Nodes

Contact_id = [3,4];
top_right_edge = find(edges(end,:)==Contact_id(1));
top_left_edge = find(edges(end,:)==Contact_id(2));
top_right_nodes = reshape(edges(1:end-1,top_right_edge),1,[]);
top_left_nodes = reshape(edges(1:end-1,top_left_edge),1,[]);
C_nodes = unique([top_right_nodes,top_left_nodes]); 
num_candidate = numel(C_nodes);
active_set = false(num_candidate,1);

%% Pre-Initiation before the step loop

u0 = zeros(ndof,1);
force_hist = zeros(n_s,1);
indentation = linspace(0,y_ini-y_final,n_s);
g0 = zeros(num_candidate,1);
N0 = zeros(2,num_candidate);

for ic = 1:num_candidate
    
    ni = C_nodes(ic);
    xi = coords(:,ni);
    u_prev = [u0(2*ni-1);u0(2*ni)];
    [g0(ic),N0(:,ic)]= gap_lin(X_bar(:,1),xi,R,u_prev);

end

%% Loading Loop

for n= 2:n_s+1
    x_c = X_bar(:,n-1);
    active_local=active_set;
    bool = true;
    k= 0 ;
    while k ==0 || bool

        nC = sum(active_local);
        index = find(active_local);
        C_c = zeros(nC,ndof);
        G_c = zeros(nC,1);
        node_list = zeros(nC,1);
        if nC~=0
            for i = 1:nC 

                idx = index(i);
                ni = C_nodes(idx);
                C_i = zeros(1,ndof);
                C_i(1,2*ni-1)= N0(1,idx);
                C_i(1,2*ni) = N0(2,idx);
                C_c(i,:) = C_i; 
                G_c(i) = g0(idx);
                node_list(i)= ni;
     
            end
            
            total_dof = ndof+nC;
            u = zeros(total_dof,1);
            freedof = setdiff((1:total_dof)',fixdof);
            z = sparse(nC,nC);
            K = [S , C_c';
                C_c, z];
            Force = [F;-G_c+C_c*u0];
            u(freedof)= K(freedof,freedof)\Force(freedof);
            lambda = u(ndof+1:end);

        else
           u = zeros(ndof,1);
           freedof = setdiff((1:ndof)',fixdof);
           u(freedof)= S(freedof,freedof)\F(freedof);
           lambda = []; 
        end
        k= k+1;
        new_active_local = active_local;
        if ~isempty(lambda)
            
            for i = 1:nC 
                idx = index(i);
                
                if lambda(i)<=0
            
                    new_active_local(idx)=true;
                else
                    new_active_local(idx)= false;
                end   
            end
        end
        g_new = zeros(num_candidate,1);
        N_new = zeros(2,num_candidate);
        for j= 1:num_candidate
            ni = C_nodes(j);
            xi = coords(:,ni);
            u_g = [u(2*ni-1);u(2*ni)];
            u_prev = [u0(2*ni-1);u0(2*ni)];
            [gx,Nx]= gap_lin(x_c,xi,R,u_prev);
            g_new(j)= gx+Nx'*u_g;
            N_new(:,j)= Nx;
            if g_new(j) < 0
            
                new_active_local(j)=true;
            end
        end

        bool = any(new_active_local~=active_local);
        active_local= new_active_local;
    end
    g0 = g_new;
    N0 = N_new;
    u0 = u(1:ndof);
    active_set= new_active_local;

    if ~isempty(lambda)

        NL = zeros(numel(lambda),1);

        for il = 1:numel(lambda)

            ni= node_list(il);
            X = coords(:,ni);
            u_node = [u(2*ni-1);u(2*ni)];
            [~, N_node]= gap_lin(x_c,X,R,u_node);
            NL(il)= N_node(2);

        end
        F_contact = sum(lambda.*NL);
    else
        F_contact = 0;
    end

    force_hist(n-1)= F_contact;
end

vm = zeros(size(Edof,1),1);
for ielem= 1:size(Edof,1)

    enodes = triangle(:,ielem);
    ex = Ex(enodes);
    ey = Ey(enodes);
    edof = Edof(ielem,2:end);
    ed = u0(edof);
    [es,et]= plants(ex,ey,ep,D,ed');
    sigma = [es(1),es(4),0;
            es(4),es(2),0;
            0,0,es(3)];
    sij = sigma - (1/3)*trace(sigma)*eye(3);
    vm(ielem) = sqrt(3/2*sum(sum(sij.*sij)));

end

ed1 = extract_ed(Edof,u0);
Ex_plot = zeros(size(Edof,1),3);
Ey_plot = zeros(size(Edof,1),3);
for i = 1:size(Edof,1)

    enodes = triangle(:,i);
    Ex_plot(i,:)= coords(1,enodes);
    Ey_plot(i,:)= coords(2,enodes);

end

%% Plotting

theta = linspace(-pi/6, -pi/2, 200); 
xc = x_c(1); yc = x_c(2);
xarc = xc + R*cos(theta);
yarc = yc + R*sin(theta);

%Force v/s Indentation Plot
figure
plot(indentation,force_hist,LineWidth=1.5)
xlabel ('Indentation [in mm]')
ylabel ('Force [in N]')
title('Total Contact Force v/s Indentation')

%Von Mises Stress Plot
figure; 
fill((Ex_plot(:,1:3)+ed1(:,1:2:6))', (Ey_plot(:,1:3)+ed1(:,2:2:6))', vm); % Factor of Safety
hold on;
plot(xarc, yarc, 'k-', 'LineWidth', 2);
view(2);
axis equal tight;
colormap("turbo");
colorbar;
c = colorbar;
title('Von Mises Stress in Deformed Mesh');

% Displacement Plot
ux = u0(1:2:end);       
uy = u0(2:2:end);       
umag = hypot(ux, uy);   
faces = triangle';      
verts = coords';   
verts_def = verts +  [ux, uy];
figure;
trisurf(faces, verts_def(:,1), verts_def(:,2), umag);
hold on;
plot(xarc, yarc, 'k-', 'LineWidth', 2);
view(2);
axis equal tight;
colormap turbo;
cb = colorbar;
title('Dispalcement in Deformed Mesh');


