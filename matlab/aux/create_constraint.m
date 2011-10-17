function [c_km1] = create_constraint(P, c_k, x_k)
% Copyright 2011, Kenny Erleben, DIKU
%--- Computes c_km1(j) = max_i{ (c_k - x_k)_i | P(i,j)>0 } ----------------
[~,N] = size(P);
d = c_k - x_k;
mask = (P>0)';
D =  repmat(d', N, 1);
D(mask==0) = - Inf;
c_km1 = max( mask .*D , [], 2);
end
