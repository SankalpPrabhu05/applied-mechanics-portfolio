function [Ke] = Taylor_hood(ex,ey,E,nu,t)

ngp = 3; % Number of Gauss Points

G = E/(2*(1+nu));
K = E/(3*(1-2*nu));

nodes = size(ex,2);
ndofs_u = nodes*2;
ndofs_p = nodes/2;

Xe= zeros(2,6);
for i = 1:6
    Xe(1,i)=ex(i);
    Xe(2,i)=ey(i);
end

xi_v = [(1/6) (1/6) (2/3);
        (1/6) (2/3) (1/6)];

Se = zeros(ndofs_u,ndofs_u);
Ce = zeros(ndofs_u,ndofs_p);
Be = zeros(ndofs_p,ndofs_p);

A = diag([2 2 1]);
H_v = ones(1,ngp)/6;

for j=1:ngp
    xi = xi_v(:,j);
    H = H_v(j); 
    [B_eps,B_div,B_dev,N_p,detFisop]= Be_Taylor_Hood_implementation(xi,Xe(:,1),Xe(:,2),Xe(:,3),Xe(:,4),Xe(:,5),Xe(:,6));
    Se = Se + B_eps'*G*A*B_dev*detFisop*H*t;
    Ce = Ce + B_div'*N_p*detFisop*H*t;
    Be = Be + N_p'*(1/K)*N_p*detFisop*H*t;
end

Ke = [Se, -Ce;
    -Ce', -Be];