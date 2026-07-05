function [funit,gamma_L] = Loadvector_interpolation(coords,edges,bound_id)

nnodes = size(coords,2);
ndof = nnodes*2;
funit = zeros(ndof,1);
gamma_L=0;

for ec = 1:size(edges,2)
    
    id = edges(end,ec);

    if id==bound_id
        n1 = edges(1,ec);
        n2 = edges(2,ec);

        x1 = coords(:,n1);
        x2 = coords(:,n2);
        L = norm(x2-x1);

        gamma_L = gamma_L+L;

        funit(2*n1)=funit(2*n1)+L/2;
        funit(2*n2)=funit(2*n2)+L/2;
    end

end

end