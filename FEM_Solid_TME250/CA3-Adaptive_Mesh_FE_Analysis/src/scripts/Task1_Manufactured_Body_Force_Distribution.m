clear all
close all
clc

%% Material Properties
L = 40e-3; % m
H = 20e-3; % m
E = 50e9; % Pa
nu = 0.3;

mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

%% Analytical Solution

syms x y real

u_x = x*y*(L-x)^3*(H-y)/(L^3*H);
u_y = x*y*(L-x)*(H-y)^3/(H^3*L);

eps_xx = diff(u_x, x);
eps_yy = diff(u_y, y);
eps_xy = (1/2) * (diff(u_x, y) + diff(u_y, x));

sigma_xx = (2*mu + lambda) * eps_xx + lambda * eps_yy;
sigma_yy = lambda * eps_xx + (2*mu + lambda) * eps_yy;
sigma_xy = 2 * mu * eps_xy;

energy_density = sigma_xx * eps_xx + sigma_yy * eps_yy + 2 * sigma_xy * eps_xy;

energy_norm_squared = int(int(energy_density, y, 0, H), x, 0, L);
energy_norm = sqrt(double(energy_norm_squared));

fprintf('||u||_a = %.10f\n', energy_norm);

save('exact_energy_norm.mat', 'energy_norm'); % To save the exact energy norm

%% Function for Body Forces
fx = -diff(sigma_xx,x)-diff(sigma_xy,y);
fy = -diff(sigma_xy,x)-diff(sigma_yy,y);
matlabFunction(fx,fy,'File','Body_force','Vars',{x,y});