% Import Variables and Draw Mesh

% ********************* Unit scheme: mm ton s MPa ********************* %   % https://www.dynasupport.com/howtos/general/consistent-units
clear
close all
clc

mesh = load("topology_medium_6node.mat");

% Get field names
fields = fieldnames(mesh);

% Loop through each field and assign it to a variable
for i = 1:numel(fields)
    field = fields{i};
    assignin('base', field, mesh.(field)); % Assign to base workspace
end

% Plot mesh
eldraw2_ext(Ex, Ey, plotpar)

mtrlpar.Emod = 20;
mtrlpar.nu = 0.45;
thick = 100;

%%
% Start timer
tic;

% Initialize Disps
a = zeros(ndofs,1);                 % Displacements vector 
aold = zeros(ndofs,1);              % Old displacements (stores the displacement of the previous iteration)
da = a-aold;                        % Delta a calculation

% Define Free and Constrained (Prescribed) DoFs
dof_F = 1:ndofs;                                                   % Define all as free
dof_C = [dof_left(1:2:end) dof_right(1:2:end) dof_corner(end)];    % Define the constrained dofs vector - ux of left and right edges and uy of corner
dof_F(dof_C) = [];                                                 % Remove the constrained from the free

% Time Stepping
ntime = 19;                         % Number of Timesteps (load increments)
tend = 1;                         % End time 
t = linspace(0, tend, ntime);       % All time steps

% Displacement control of ux dofs of right side
umax = -15;
uu = umax*(t./tend);

% Post-Processing variables
P = zeros(size(uu));                  % Applied forces vector
RF = zeros(size(uu));                 % Reaction forces vector 
K = spalloc(ndofs, ndofs, 20*ndofs);  % Defines K as a sparse matrix
Sv = zeros(5,3,nelem);                % 3D Matrix to store Cauchy stress of each element (1st is stress components, 2nd is gauss point, 3rd dimension is element )
Pv = zeros(4,3,nelem);                % 3D Matrix to store 1st Piola-K stress of each element
Mises = zeros(1,3,nelem);                % 3D Matrix to store von mises at each gauss point

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
            ae = a(Edof(iel,2:end));
            ex = Ex(iel,:);             % Sums element current displacements to initial coordinates
            ey = Ey(iel,:);
            [fe_int, Ke_int, ~, ~, ~, ~] = TRIA6_LDef(ae, ex, ey, mtrlpar, thick);
            K(Edof(iel,2:end), Edof(iel,2:end)) = K(Edof(iel,2:end), Edof(iel,2:end)) + Ke_int;
            fint(Edof(iel,2:end)) = fint(Edof(iel,2:end)) + fe_int;         % Sums element internal forces         
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

% Stop the timer and display the elapsed time
elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);

% Loop over elements to obtain stresses
for iel = 1:nelem
    ae = a(Edof(iel,2:end));
    ex = Ex(iel,:);
    ey = Ey(iel,:);
    [~, ~, ~, Pv_m, Sv_m,~] = TRIA6_LDef(ae, ex, ey, mtrlpar, thick);
    Sv(:,:,iel) = Sv_m;
    Pv(:,:,iel) = Pv_m;
    
    % Loop over gauss points to obtain their von Mises equivalent stress
    for gp = 1:size(Sv,2)
        Sigma = [Sv_m(1,gp) Sv_m(3,gp) 0; Sv_m(4,gp) Sv_m(2,gp) 0; 0 0 Sv_m(5,gp)];          % Cauchy-Stress Tensor 3x3
        sij = Sigma - 1/3 * trace(Sigma) * eye(3);
        vm = sqrt(3/2*sum(sum(sij.*sij)));
        Mises(1,gp,iel) = vm;
    end
end

Mises_avg = mean(Mises, 2);     % Compute the mean along the second dimension
Mises_avg = squeeze(Mises_avg); % Remove the singleton dimension

%plot(uu, RF, '-') % Plotting internal force as a function of applied displacement during simulation

figure
plot(uu, P, '-')  %
hold on
% Find the index of the maximum value
max_value = P(end);

% Highlight the maximum value on the plot
plot(uu(end), P(end), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
hold off
legend('Applied Force', 'Interpreter','latex', Location='northwest')

% Add text annotation for the maximum value
text(uu(end), P(end), sprintf('Max: %.2f', max_value), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
xlabel('Displacement of Right Edge [mm]', Interpreter='latex')
ylabel('Force [N]', Interpreter='latex')
title('External Work',Interpreter='latex')

fprintf('Total Reaction Force is Fx = %.2f [N]. \n', RF(end));
fprintf('Max Von Mises Stress is = %.2f [MPa]. \n', max(Mises_avg));

% Plot mesh
ed = extract_ed(Edof,a);
%eldisp2_ext(Ex, Ey, ed, plotpar, 1)

figure('Position', [100, 100, 600*0.85, 600]); % [left, bottom, width, height]
fill((Ex(:,1:3)+ed(:,1:2:6))', (Ey(:,1:3)+ed(:,2:2:6))', Mises_avg); % Factor of Safety

colormap("turbo");
colorbar;
c = colorbar;
c.Label.String = 'Stress [MPa]';
title('Von Mises Stress - Compressive Case',Interpreter='latex');
set(gca, 'XTick', [], 'YTick', []) 
