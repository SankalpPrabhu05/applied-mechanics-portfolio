clear all
close all
clc

nx = 50;
ny = 25;

L = 0.04;
H = 0.02;
x = linspace(0,0.04,nx);
y = linspace(0,0.02,ny);
[X,Y] = meshgrid(x,y);

fx = zeros(size(X));
fy = zeros(size(Y));

for i = 1:numel(X)
    [fx(i), fy(i)] = Body_force(X(i), Y(i));
end

Xmm = X*1e3;
Ymm = Y*1e3;

Fmag = sqrt(fx.^2 + fy.^2);

scale = max(Fmag(:));
fxn = fx ./ scale;
fyn = fy ./ scale;


figure
contourf(Xmm, Ymm, Fmag, 20, 'LineStyle', 'none')
colormap('turbo')
hold on

scatter(Xmm(:), Ymm(:), 20, Fmag(:), 'filled')  % magnitude
quiver(Xmm, Ymm, fxn, fyn, 'k')  

axis equal
grid on
xlabel('x [mm]')
ylabel('y [mm]')
title('Body Force Distribution')
