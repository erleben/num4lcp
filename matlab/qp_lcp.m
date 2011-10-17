function [x err] = qp_lcp(A, b)
% Copyright 2011, Kenny Erleben, DIKU

%--- Observe A must be symmetric for QP reformulation to make sense -----
  x = quadprog(A, b,  - eye(size(A)),  zeros(size(b)) );
  
  y = abs(A*x + b);
  err = x'*y;
end
