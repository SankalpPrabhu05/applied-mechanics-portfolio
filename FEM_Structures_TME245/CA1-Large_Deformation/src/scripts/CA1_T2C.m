clear
close all

% Create the structure
mtrlpar.Emod = 20;
mtrlpar.nu = 0.45;

ex = [0 3 0 1.5 1.5 0];
ey = [0 0 2 0 1 1];
def_x = [6 7 4.5 6.4 5.6 5.2];
def_y = [0.7 2.3 1.8 1.2 2 1.1];
edx = def_x - ex;
edy = def_y - ey;
ae = zeros(12,1);

for i=1:6
    ae(2*i-1) = edx(i);
    ae(2*i) = edy(i);
end
t = 100;

[fe_int,Ke_int,Fv_m,Pv_m,Sv_m, dPv_dFv_m] = TRIA6_LDef(ae, ex, ey, mtrlpar, t);
display(Fv_m)
display(Pv_m)
display(fe_int)
display(Ke_int)