clear all
close all
clc

%% Material Properties

E = 50e3; %Mpa
nu = 0.3;
rho = 8e-6; %kg/mm3
T = 5e-3; %sec
dt = 20e-6; %sec
v_ini = 10e3; %mm/sec
edge_length = 2; %mm
thickness = 1; % Assumed
R = 40; %mm
x_c = [0;60];
H = 20; %mm
L = 2*H; %mm

%% Meshing

ny = H/edge_length;
order = 1;
[coords, edges, triangle] = mesh(ny,order,false);
coords = 1000.*coords;
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
M = sparse(ndof,ndof);
F = zeros(ndof,1);

for ielem = 1:size(Edof,1)
    
    enodes = triangle(:,ielem);
    ex = Ex(enodes);
    ey = Ey(enodes);
    tri_coords = [ex;ey;1,1,1]'; 
    Ae = 0.5*abs(det(tri_coords));
    [Ke,~] = plante(ex,ey,ep,D);    
    Me = Ae*rho/3*eye(6);
    edof = Edof(ielem,2:end);
    S(edof,edof) = S(edof,edof)+Ke;
    M(edof,edof) = M(edof,edof) + Me;

end

%% Boundary Conditions

fixed_dofs = [];
left_edge = find(edges(end,:)==5);

for le = 1:size(left_edge,2)
    n1 = edges(1,left_edge(1,le));
    n2 = edges(2,left_edge(1,le));
    fixed_dofs = [fixed_dofs,2*n1-1,2*n2-1];
end

fixdof = unique(fixed_dofs);

%% Defining Contact Nodes

Contact_id = [3,4];
top_right_edge = find(edges(end,:)==Contact_id(1));
top_left_edge = find(edges(end,:)==Contact_id(2));
top_right_nodes = reshape(edges(1:end-1,top_right_edge),1,[]);
top_left_nodes = reshape(edges(1:end-1,top_left_edge),1,[]);
C_nodes = unique([top_right_nodes,top_left_nodes]); 
num_candidate = numel(C_nodes);
active_set = false(num_candidate,1);

%% Initiations before loop

v0 = zeros(ndof,1);
v0(2:2:end)= v_ini;
v0(fixdof)=0;
u0 = zeros(ndof,1);
n_s = 250;
U = zeros(ndof,n_s);
V = zeros(ndof,n_s);
int_st = zeros(n_s,1);
kin_eng = zeros(n_s,1);
T_s = linspace(0,T,n_s);
force_hist = zeros(n_s,1);

%% Initial Gap
g0 = zeros(num_candidate,1);
N0 = zeros(2,num_candidate);

for ic = 1:num_candidate
    
    ni = C_nodes(ic);
    xi = coords(:,ni);
    u_prev = [u0(2*ni-1);u0(2*ni)];
    [g0(ic),N0(:,ic)]= gap_lin(x_c,xi,R,u_prev);

end

%% Time Loop

for t = 2:n_s+1

    u_old = u0;
    v_old = v0;
    F_eff = F+((2/(dt^2))*M-0.5*S)*u_old+((2/dt)*M)*v_old;
    S_eff = ((2/(dt^2))*M+0.5*S);
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
            for j = 1:nC 
                idx = index(j);
                ni = C_nodes(idx);
                C_i = zeros(1,ndof);
                C_i(1,2*ni-1)= N0(1,idx);
                C_i(1,2*ni) = N0(2,idx);
                C_c(j,:) = C_i; 
                G_c(j) = g0(idx);
                node_list(j)= ni;
            end
            
            total_dof = ndof+nC;
            u = zeros(total_dof,1);
            freedof = setdiff((1:total_dof)',fixdof);
            z = sparse(nC,nC);
            K = [S_eff , C_c';
                C_c, z];
            Force = [F_eff;-G_c+C_c*u_old];
            u(freedof)= K(freedof,freedof)\Force(freedof);
            lambda = u(ndof+1:end);
            u_n = u(1:ndof);

        else
           u = zeros(ndof,1);
           freedof = setdiff((1:ndof)',fixdof);
           u(freedof)= S_eff(freedof,freedof)\F_eff(freedof);
           lambda = []; 
           u_n = u;
        end
        k= k+1;
        new_active_local = active_local;
        if ~isempty(lambda) 
            for l = 1:nC 
                idx = index(l);
                if lambda(l)<=0
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
            u_g = [u_n(2*ni-1);u_n(2*ni)];
            u_prev = [u_old(2*ni-1);u_old(2*ni)];
            [gx,Nx]= gap_lin(x_c,xi,R,u_prev);
            g_new(j)= gx+Nx'*u_g;
            N_new(:,j)=Nx;
            if g_new(j) < 0
            
                new_active_local(j)=true;

            end
        end        
        bool = any(new_active_local~=active_local);
        active_local= new_active_local;
    end

    u0 = u_n;
    for j = 1:num_candidate
        ni = C_nodes(j);
        xi = coords(:,ni);
        u_conv = [u0(2*ni-1); u0(2*ni)];
        [gx, N0(:,j)] = gap_lin(x_c, xi, R, u_conv);
        g0(j)= gx+N0(:,j)'*u_conv;
    end
    U(:,t-1)= u0;
    active_set= new_active_local;
    v0= (2/(dt))*(u0-u_old)-v_old;
    V(:,t-1)= v0;
    int_st(t-1,1)=0.5*u0'*S*u0;
    kin_eng(t-1,1)=0.5*v0'*M*v0;
    if ~isempty(lambda)
        Ny = zeros(numel(lambda),1);
        for il = 1:numel(lambda)
            idx = index(il);
            Ny(il)= N0(2,idx);
        end
        F_contact = sum(lambda.*Ny);
    else
        F_contact = 0;
    end
    force_hist(t-1)= F_contact;
end

%% Plotting

% Force History

figure
plot(T_s,force_hist)
xlabel ('Indentation [in mm]')
ylabel ('Force [in N]')
title('Total Contact Force v/s Indentation')

% Simulation, Comment it for Faster run
 
% scale=1;
% frame_pause = 0.02;
% vert = coords';
% Face = triangle';
% figure
% h.patch=patch('Faces', Face, 'Vertices', vert, ...
%                 'FaceColor','none', ...       % no face fill
%                 'EdgeColor',[0 0 0], ...     % black mesh lines
%                 'LineWidth',0.5);
% 
% axis equal;hold on;
% theta = linspace(-pi/6, -pi/2, 200); 
% xc = x_c(1); yc = x_c(2);
% xarc = xc + R*cos(theta);
% yarc = yc + R*sin(theta);
% h.arc = plot(xarc, yarc, 'r', 'LineWidth', 1.5);
% xlabel('x (mm)'); ylabel('y (mm)');
% title('Deformed mesh (edges only)');
% drawnow;
% 
% for i= 1:n_s
%     uk = U(:,i);
%     ux = uk(1:2:end);
%     uy = uk(2:2:end);
%     vertk= coords'+scale.*[ux, uy];
%     set(h.patch,'Vertices',vertk);
%     title(sprintf('Time step %d / %d', i, n_s));
%     drawnow;
%     pause(frame_pause)
% end

% Energy over Time

figure
plot(T_s,int_st,LineWidth=1.5,DisplayName='Internal Strain Energy');
hold on
plot(T_s,kin_eng,LineWidth=1.5,DisplayName='Kinetic Energy');
plot(T_s,int_st+kin_eng,LineWidth=1.5,DisplayName='Total Energy');
xlabel ('Time [in sec]')
ylabel ('Energy [in mJ]')
legend (Location="best")
title('Development of Energy over Time')

%Deformation of Block for Selected Time Step

time_step =round([1e-3/dt,2e-3/dt,3e-3/dt,4e-3/dt,5e-3/dt]);
scale=1;
vert = coords';
Face = triangle';
figure
patch('Faces', Face, 'Vertices', vert,'FaceColor','none','EdgeColor',[0 0 0], 'LineWidth',0.5);
axis equal;hold on;
theta = linspace(-pi/6, -pi/2, 200); 
xc = x_c(1); yc = x_c(2);
xarc = xc + R*cos(theta);
yarc = yc + R*sin(theta);
plot(xarc, yarc, 'r', 'LineWidth', 1.5);
xlabel('x (mm)'); ylabel('y (mm)');
title ('Initial Mesh at Time 0 ms')

for i=1:numel(time_step)

    k = time_step(i);
    uk = U(:,k);
    ux = uk(1:2:end);
    uy = uk(2:2:end);
    vertk= coords'+scale.*[ux, uy];
    figure
    patch('Faces', Face, 'Vertices', vertk,'FaceColor','none','EdgeColor',[0 0 0], 'LineWidth',0.5);
    axis equal;hold on;
    theta = linspace(-pi/6, -pi/2, 200); 
    xc = x_c(1); yc = x_c(2);
    xarc = xc + R*cos(theta);
    yarc = yc + R*sin(theta);
    plot(xarc, yarc, 'r', 'LineWidth', 1.5);
    xlabel('x (mm)'); ylabel('y (mm)');
    title (sprintf('Deformed mesh at Time %d ms',k*dt*1000))

end