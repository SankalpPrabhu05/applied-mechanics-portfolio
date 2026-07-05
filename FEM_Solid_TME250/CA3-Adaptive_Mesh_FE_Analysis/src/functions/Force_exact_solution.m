function [Fe] = Force_exact_solution(ex,ey,mu,lambda,L,H,thickness)

Xe= zeros(2,3);
for i = 1:3
    Xe(1,i)=ex(i);
    Xe(2,i)=ey(i);
end

xi_v = [(1/6) (1/6) (2/3);
        (1/6) (2/3) (1/6)];
ngp = 3;
H_v = ones(1,ngp)/6;
Fe = zeros(6,1);

for j=1:ngp
   
    xi = xi_v(:,j);
    weight = H_v(j);
    [X,detJ]=Gauss(xi,Xe(:,1),Xe(:,2),Xe(:,3));
    N = [1-xi(1)-xi(2);xi(1);xi(2)];
    [fx,fy] = Body_force(X(1),X(2));
    f_vec = [fx;fy];

    for i = 1:3
        
        Fe(2*i-1) = Fe(2*i-1) + weight * detJ * thickness * N(i) * f_vec(1);
        Fe(2*i)   = Fe(2*i)   + weight * detJ * thickness * N(i) * f_vec(2);
    end
end

