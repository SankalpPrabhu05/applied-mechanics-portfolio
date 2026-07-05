% Define Input Variables 

% ********************* Unit scheme: mm ton s MPa ********************* %  
clear
close all
clc

E = 20;
nu = 0.45;

%%
F11 = linspace(0.5,1.5,200);
S11_Yeoh = zeros(size(F11));
S11_NH = zeros(size(F11));
P11_Yeoh = zeros(size(F11));
%Pv = zeros(4,size(F11,2));

S = zeros(3,3,size(F11,2));
for i=1:size(F11,2)
    Fv = [F11(i); 1; 0; 0];
    [pv, ~, sigma] = Yeoh_func_S33(Fv, E, nu);
    [~, s11nh] = NH_func(Fv, E, nu);
    S11_Yeoh(i) = sigma(1,1);
    S11_NH(i) = s11nh(1,1);
    P11_Yeoh(i) = pv(1,1);
    %S(:,:,i) = sigma;
end

figure;
plot(F11, S11_NH, 'LineWidth', 1);
hold on;
plot(F11, S11_Yeoh, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'r');
hold off;
legend('Neo-Hookean', 'Yeoh', 'Location', 'southeast', 'FontSize', 12);
xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('Cauchy Stress $\sigma$ [MPa]', 'Interpreter', 'latex');
ax = gca;
ax.GridLineStyle = '--';
%ax.GridColor = [0.8, 0.8, 0.8]; % Light grey color
grid on;
title('Nearly Incompressible Material Models');