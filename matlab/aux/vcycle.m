function [ x ] = vcycle( k, b_k, x_k, c_k, As, Ps )
% Copyright 2011, Kenny Erleben, DIKU

A_k = As{k};

if k==1
  
  max_iter = 3;
  lambda   = 1.4;
  x        = aux_psor(A_k, b_k, x_k, c_k, max_iter, lambda);
  
else
  
  P      = Ps{k};    %P      = create_prolongation(A_k);
  A_km1  = As{k-1};  %A_km1  = P' * A_k * P;

  r_k    = A_k * x_k + b_k;
  b_km1  = P' * r_k;
  c_km1  = create_constraint(P,c_k,x_k);
  x_km1  = zeros(size(b_km1));
  
  max_iter = 3;
  lambda   = 1;
  x_km1    = aux_psor(A_km1, b_km1, x_km1, c_km1, max_iter, lambda);
  
  x_km1    = vcycle( k-1, b_km1, x_km1, c_km1, As, Ps );
  
  x_k      = x_k + P * x_km1;
  
  max_iter = 3;
  lambda   = 1.7;
  x        = aux_psor(A_k, b_k, x_k, c_k, max_iter, lambda);
  
end
end