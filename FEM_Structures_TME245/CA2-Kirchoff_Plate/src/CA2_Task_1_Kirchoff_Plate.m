clear
close all
clc

% ********************* Unit scheme SI: N m kg s Pa ********************* %

%% Property

E = 210e9;                      % Young's Modulus
nu = 0.3;                       % Poisson's Ration
Sy = 450e6;                     % Yield Strength
G_modulus = E/(2*(1+nu));       % Transverse Modulus
G = G_modulus*eye(2);
D = hooke(1,E,nu);              % Material Matrix
rho = 1e3;                      % Water density for hydrostatic pressure
h = 10e-3;                      % Thickness
ep = [1 h];
%% Meshing

xmin = 2; xmax = 4; ymin = 0; ymax = 2;
nelx = 40; nely = 40;
[mesh, coord, Edof_ip, Edof_op] = rectMesh(xmin, xmax, ymin, ymax, nelx, nely);

% Extract x and y coordinates for each node in mesh
Ex = reshape(coord(mesh, 1), 4, [])'; % Get x-coordinates from coord using mesh indices
Ey = reshape(coord(mesh, 2), 4, [])'; % Get y-coordinates from coord using mesh indices

% Edof arrays for each edge of the square
LS_nodes = find(abs(coord(:,1)-xmin)<eps);
RS_nodes = find(abs(coord(:,1)-xmax)<eps);
BS_nodes = find(abs(coord(:,2)-ymax)<eps);
TS_nodes = find(abs(coord(:,2)-ymin)<eps);
TT_nodes = find(abs(coord(TS_nodes,1) >= 2.5 & abs(coord(TS_nodes,1)) <= 3.5));

Edof_ip = full(Edof_ip);
Edof_op = full(Edof_op);

% Find all DoFs and count it for in-plane and out-of-plane
ndof_ip = numel(unique(Edof_ip));
ndof_op = numel(unique(Edof_op));
nnode = ndof_ip/2;                 % Number of nodes

Dof_ip = zeros(nnode,2);           % Rows are nodes and columns are dofs (x,y)
Dof_op = zeros(nnode,3);           % Rows are nodes and columns are dofs (w,th1,th2) 

% Populate Dof matrices 
for i=1:nnode       
    Dof_ip(i,1) = 2*i-1;
    Dof_ip(i,2) = 2*i;
    Dof_op(i,1) = 3*i-2;
    Dof_op(i,2) = 3*i-1;
    Dof_op(i,3) = 3*i;
end

% Retrieve the dofs of edges of interest for both in-plane and out-of-plane
RS_ip = Dof_ip(RS_nodes,:);
TS_ip = Dof_ip(TS_nodes,:);
LS_ip = Dof_ip(LS_nodes,:);
BS_ip = Dof_ip(BS_nodes,:);
TT_ip = Dof_ip(TT_nodes,:);           % Top traction in-plane DOFs
RS_op = Dof_op(RS_nodes,:);
LS_op = Dof_op(LS_nodes,:);
BS_op = Dof_op(BS_nodes,:);

nelem = size(Edof_ip,1);
%% Load and Boundary Conditions

l = 1;                                % Length of the edge
g = 9.806;                            % Gravity
P = (150*g)/(l*h);                    % Uniform Pressure
eq = g*[0; 0];                        % Body Forces
faa = zeros(ndof_ip,1);               % Global In-Plane Load vector
fww = zeros(ndof_op,1);               % Global Out-of-Plane Load vector

% Calculate Pressure load for In-Plane
for i = 1:size(TT_nodes, 1)-1           % For the number of elements on the Top T edge (# of nodes -1)  
    pos1 = coord(TT_nodes(i),:);
    pos2 = coord(TT_nodes(i+1),:);
    Ep = pressure_load(pos1, pos2, P, h);   % Element Pressure
    dofs = reshape(TT_ip(i:i+1,:)', 1, []); % Flattens the dofs matrix
    faa(dofs) = faa(dofs) + Ep;
end

% Create BC vector for In-Plane
ux_0 = [LS_ip(:,1); BS_ip(:,1); RS_ip(:,1)];    % Creates a vector with all DoFs that will have ux = 0
uy_0 = [LS_ip(:,2); BS_ip(:,2); RS_ip(:,2)];    % Creates a vector with all DoFs that will have uy = 0
pdis_ip = unique([ux_0; uy_0]);                 % Concatenates prescribed DoFS in a vector
pd_0_ip = zeros(size(pdis_ip,1),1);             % Creates a vector of zeroes that correspond to the prescribed displacements
bc_ip = [pdis_ip pd_0_ip];                      % Column 1 is DoF, Column 2 is prescribed value

% Create BC vector for Out-Plane
uz_0 = [LS_op(:,1); BS_op(:,1); RS_op(:,1)];    % w = 0
tx_0 = [LS_op(:,2); BS_op(:,2); RS_op(:,2)];    % theta_x = 0
ty_0 = [LS_op(:,3); BS_op(:,3); RS_op(:,3)];    % theta_y = 0
pdis_op = unique([uz_0; tx_0; ty_0]);           % Concatenates prescribed DoFS in a vector
pd_0_op = zeros(size(pdis_op,1),1);             % Creates a vector of zeroes that correspond to the prescribed displacements
bc_op = [pdis_op pd_0_op];                      % Column 1 is DoF, Column 2 is prescribed value

%% Stifness Assembly + Out-of-Plane Loads

% Assembly Global Kaa and Kww Matrices and fww matrix (out-of-plane pressure)
Kaa = zeros(ndof_ip, ndof_ip);     % Global In-Plane Stiffness Matrix
Kww = zeros(ndof_op, ndof_op);     % Global Bending Stiffness Matrix

for i=1:nelem
    edof_ip = Edof_ip(i,2:end);
    edof_op = Edof_op(i,2:end);
    ex = Ex(i,:);
    ey = Ey(i,:);
    yel = sum(ey)/size(ey,2);
    q = rho*g*yel;
    [Kew, Kea, fwe_ext, ~] = Full_KQuad_Func(ex,ey,ep,D,q,eq);
    Kaa(edof_ip, edof_ip) = Kaa(edof_ip, edof_ip) + Kea;
    Kww(edof_op, edof_op) = Kww(edof_op, edof_op) + Kew;
    fww(edof_op) = fww(edof_op) + fwe_ext;
end

%% Solve In-Plane Problem

% Solution of the [Kaa]{Ua} = {fa} system
fdof_ip = 1:ndof_ip;                        % Vector that contains all degrees of freedom
U_ip = zeros(ndof_ip,1);                    % Displacement vector
pdof_ip = bc_ip(:,1);                       % Retrieves dofs that have a prescribed value
dp_ip = bc_ip(:,2);                         % Retrieves prescribed value of dofs
fdof_ip(pdof_ip) = [];                      % DOFs with U = 0 are crossed-out

% Solve
s_ip = Kaa(fdof_ip,fdof_ip)\(faa(fdof_ip)-Kaa(fdof_ip,pdof_ip)*dp_ip); % Condensed problem to solve for unknown displacements s = [Kff]^-1 * (ff - Kfc*ac) 

% Populate displacements vector and calculate reaction force 
U_ip(pdof_ip) = dp_ip;
U_ip(fdof_ip) = s_ip;
Ra = Kaa*U_ip - faa;

%% Solve Out-of-Plane Problem

% Solution of the [Kww]{Uw} = {fw} system
fdof_op = 1:ndof_op;                        % Vector that contains all degrees of freedom
U_op = zeros(ndof_op,1);                    % Displacement vector
pdof_op = bc_op(:,1);                       % Retrieves dofs that have a prescribed value
dp_op = bc_op(:,2);                         % Retrieves prescribed value of dofs
fdof_op(pdof_op) = [];                      % DOFs with U = 0 are crossed-out

Kww_FF = Kww(fdof_op,fdof_op);
fww_FF = fww(fdof_op);
% Solve
s_op = Kww_FF\(fww_FF-Kww(fdof_op,pdof_op)*dp_op); 

% Populate displacements vector and calculate reaction force 
U_op(pdof_op) = dp_op;
U_op(fdof_op) = s_op;
Ro = Kww*U_op - fww;
%% Extracting Displacements and Stresses

Ux = U_ip(1:2:end);
Uy = U_ip(2:2:end);
Uw = U_op(1:3:end);
Tx = U_op(2:3:end);
Ty = U_op(3:3:end);
E_ip = extract_ed(Edof_ip, U_ip);
E_op = extract_ed(Edof_op, U_op);

%% Plotting Displacements

% Compute deformed nodal coordinates
x_def_ip = Ex + 600000*E_ip(:,1:2:end);
y_def_ip = Ey + 600000*E_ip(:,2:2:end);

E_ip_mag =  sqrt(sum(E_ip.^2, 2));

% In plane
figure('Position', [900, 100, 720, 720]); % [left, bottom, width, height]
fill(x_def_ip', y_def_ip', E_ip_mag); 
colormap(turbo);
c = colorbar;
xlabel('X [m]', 'Interpreter', 'latex');
ylabel('Y [m]', 'Interpreter', 'latex');
ylabel(c, 'Displacements [m]', Interpreter='latex', FontSize=14)
title('In-Plane Displacements',Interpreter='latex', FontSize=16);
set(gca, 'YDir', 'reverse')

% Compute deformed nodal coordinates
x_def = coord(:,1);% + Ux;
y_def = coord(:,2);% + Uy;
z_def = Uw; % out-of-plane displacement

% Reshape the deformed coordinates to a grid (nely+1 by nelx+1)
x_def_grid = reshape(x_def, [nely+1, nelx+1]);
y_def_grid = reshape(y_def, [nely+1, nelx+1]);
z_def_grid = reshape(z_def, [nely+1, nelx+1]);

% Plot the deformed configuration as a 3D surface
figure('Position', [100, 100, 720, 720]); % [left, bottom, width, height]
surf(x_def_grid, y_def_grid, z_def_grid);
colormap turbo;
colorbar_handle = colorbar;
xlabel('X [m]', Interpreter='latex');
ylabel('Y [m]', Interpreter='latex');
ylabel(colorbar_handle, 'Displacement [m]', 'Interpreter', 'latex', FontSize=14);       % Set the label for the colorbar
zlabel('Displacement [m]', Interpreter='latex');
title('Out of Plane Displacements',Interpreter='latex', FontSize=16);
zlim([min(z_def_grid(:)) 0.035]); % Adjust the upper limit of z axis
view(3);