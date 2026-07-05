% Import Variables and Draw Mesh

% ********************* Unit scheme: mm ton s MPa ********************* %   % https://www.dynasupport.com/howtos/general/consistent-units
clear
close all
clc

mesh = load("topology_coarse_3node.mat");

% Get field names
fields = fieldnames(mesh);

% Loop through each field and assign it to a variable
for i = 1:numel(fields)
    field = fields{i};
    assignin('base', field, mesh.(field)); % Assign to base workspace
end

% Plot mesh
figure
eldraw2(Ex,Ey, plotpar);

Emod = 20;
nu = 0.45;
thick = 100;

% Material Matrix
D = (Emod)/((1+nu)*(1-2*nu)) * [1-nu nu 0; 
                                nu 1-nu 0;
                                0 0 1/2*(1-2*nu)];
%%

% Initialize Disps
a = zeros(ndofs,1);                 % Displacements vector 
aold = zeros(ndofs,1);              % Old displacements (stores the displacement of the previous iteration)
da = a-aold;                        % Delta a calculation


% Define Free and Constrained (Prescribed) DoFs
dof_F = 1:ndofs;                                                   % Define all as free
dof_C = [dof_left(1:2:end) dof_right(1:2:end) dof_corner(end)];    % Define the constrained dofs vector - ux of left and right edges and uy of corner
dof_F(dof_C) = [];                                                 % Remove the constrained from the free

% Time Stepping
ntime = 20;                         % Number of Timesteps (load increments)
tend = 100;                         % End time 
t = linspace(0, tend, ntime);       % All time steps

% Displacement control of ux dofs of right side
umax = -15;
uu = umax*(t./tend);

% Post-Processing variables
P = zeros(size(uu));                 % Applied forces vector
RF = zeros(size(uu));                 % Reaction forces vector (?) 
K = spalloc(ndofs, ndofs, 20*ndofs); % Defines K as a sparse matrix

% Initialize internal and external force
fint = zeros(ndofs,1);
fext = zeros(ndofs,1);

% Tolerance value for Newton Iteration
tol = 1e-6;

%% Newton Iterations

for i=1:ntime
    % Guess for unknown displacement field
    a(dof_F) = aold(dof_F) + da(dof_F);     % a(dof_F) since we're only solving for the Free DoFs

    % Upadte prescribed values
    a(dof_right(1:2:end)) = uu(i);          % Odd DoFs of right side receive applied displacement

    % Newton iteration to find unknown displacements
    unbal = 1e10;                           % Define large initial value for unbalance
    niter = 0;                              % Newton iterations counter variable

    while unbal > tol
        % Initialize Tangent stiffness and internal force vectors
        K = K.*0;               
        fint = fint .* 0;

        % Loop over elements
        for iel = 1:nelem
            ed = a(Edof(iel,2:end))';
            ex = Ex(iel,:);                                                    % Sums element current displacements to initial coordinates
            ey = Ey(iel,:);
            [strs, strn] = plants(ex, ey, [2 thick], D, ed);                   % Extracts strains and stresses
            Ke = plante(ex, ey, [2 thick], D);
            K(Edof(iel,2:end), Edof(iel,2:end)) = K(Edof(iel,2:end), Edof(iel,2:end)) + Ke;
            f_int = plantf(ex,ey,[2 thick], strs)';
            fint(Edof(iel,2:end)) = fint(Edof(iel,2:end)) + f_int;         % Sums element internal forces         
        end

        %Unbalance equation
        g_F = fint(dof_F) - fext(dof_F);    % Calculates force balance at free DoFs
        unbal = norm(g_F);                  % Magnitude of unbalance vector
        if unbal > tol
            % Newton update displacements
            a(dof_F) = a(dof_F) - K(dof_F, dof_F) \ g_F;    % Solve [K]^-1 * force unbalance
        end

        niter = niter + 1;                  % Update Newton iteration counter

        % Print loadstep + newton residual information
        fprintf('Step: %d | ', i); fprintf('Global iter: %4.0f | ', niter); fprintf('Global residual: %4.2e \n', unbal);

        if niter>20
            disp('No convergence in Newton iteration')
            break
        end
    end

    % Save data for post-processing when Newton iteration has converged in the time increment
    u_rightside(:,i) = a(dof_right(1:2:end));
    
    % Save how large the displacement has changed during the timestep
    da=a-aold;
    aold = a;
    
    P(i) = sum(fint(dof_right(1:2:end)));
    RF(i) = sum(fint(dof_left(1:2:end)));
    
    plot(uu, RF, '-') % Plotting internal force as a function of applied displacement during simulation
    hold on
    plot(uu, P, '-')  %
    hold off
    drawnow
   
end

plot(uu, P, '-')  %
hold on
% Highlight the maximum value on the plot
plot(uu(end), P(end), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
hold off
legend('Applied Force', 'Interpreter','latex',Location='northwest')

% Add text annotation for the maximum value
text(uu(end), P(end), sprintf('Max: %.2f', P(end)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('Displacement of Right Edge [mm]', Interpreter='latex')
ylabel('Force [N]', Interpreter='latex')
title('External Work',Interpreter='latex')
ax = gca;
ax.GridLineStyle = '--';
grid on;

fprintf('Total Reaction Force is Fx = %.2f [N]. \n', RF(end));

% Plotting the deformed mesh
E_d = extract_ed(Edof,a);
figure
eldisp2(Ex, Ey, E_d, plotpar, 1);