function [g,N]= gap_lin(X_bar,X,R,u)
a = norm((X+u)-X_bar);
b = (X+u)-X_bar;
N = b/a;
g = a-R-N'*u;

end