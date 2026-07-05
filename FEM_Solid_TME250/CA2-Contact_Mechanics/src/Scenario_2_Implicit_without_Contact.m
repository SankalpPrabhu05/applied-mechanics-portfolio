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
top_right_edge = find(edges(end,:)==3);

for be = 1:size(top_right_edge,2)
    n1 = edges(1,top_right_edge(1,be));
    n2 = edges(2,top_right_edge(1,be));
    fixed_dofs = [fixed_dofs,2*n1,2*n2];
end

top_left_edge = find(edges(end,:)==4);

for be = 1:size(top_left_edge,2)
    n1 = edges(1,top_left_edge(1,be));
    n2 = edges(2,top_left_edge(1,be));
    fixed_dofs = [fixed_dofs,2*n1,2*n2];
end

fixdof = unique(fixed_dofs);
freedof = setdiff((1:ndof)',fixdof);

%% Initiations before loop

v = zeros(ndof,1);
v(2:2:end)= v_ini;
v(fixdof)=0;
u = zeros(ndof,1);
n_s = 250;
U = zeros(ndof,n_s);
V = zeros(ndof,n_s);
int_st = zeros(n_s,1);
kin_eng = zeros(n_s,1);
T_s = linspace(0,T,n_s);

%% Time Loop

for i = 1:n_s

    u_old = u;
    v_old = v;
    F_eff = F+((2/(dt^2))*M-0.5*S)*u_old+((2/dt)*M)*v_old;
    S_eff = ((2/(dt^2))*M+0.5*S);
    u(freedof) = S_eff(freedof,freedof)\F_eff(freedof);
    U(:,i)= u;
    v= (2/(dt))*(u-u_old)-v_old;
    V(:,i)= v;
    int_st(i,1)=0.5*u'*S*u;
    kin_eng(i,1)=0.5*v'*M*v;
end

%% Plotting
% Simulation, Comment it for faster Run
scale=1;
vert = coords';
Face = triangle';
figure
h.patch=patch('Faces', Face, 'Vertices', vert,'FaceColor','none','EdgeColor',[0 0 0], 'LineWidth',0.5);
axis equal;
xlabel('x (mm)'); ylabel('y (mm)');
title('Deformed mesh (edges only)');
drawnow;

for k= 1:n_s

    uk = U(:,k);
    ux = uk(1:2:end);
    uy = uk(2:2:end);
    vertk= coords'+scale.*[ux, uy];
    set(h.patch,'Vertices',vertk);
    title(sprintf('Time step %d / %d', k, n_s));
    drawnow;

end

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
axis equal;
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
    axis equal;
    xlabel('x (mm)'); ylabel('y (mm)');
    title (sprintf('Deformed mesh at Time %d ms',k*dt*1000))

end