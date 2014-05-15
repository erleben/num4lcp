function [ x ] = aux_psor(A, b, x, c, max_iter, lambda)
% Copyright 2011, Kenny Erleben, DIKU
N = length(b); 
iter = 0;
while iter < max_iter
  for i=1:N
    old_xi = x(i);
    ri     = b(i) + A(i,:)*x;
    Aii    = A(i,i);
    x(i) = max( c(i), old_xi - lambda*(ri / Aii) );
  end
  iter = iter + 1;
end
end