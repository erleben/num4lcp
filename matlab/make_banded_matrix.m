function A = make_banded_matrix(B,N)
% A = MAKE_BANDED_MATRIX(B,N): Make random banded diagonal symmetric matrix
%
% INPUT:
%
%   B -- Number of off-diagonal bands
%   N -- Number of variables
%
% OUTPUT:
%
%   A -- A banded symmetric positive definite matrix
%
% Copyright 2011, Kenny Erleben, DIKU
R    = - rand(B,1);
A    = eye(N,N);
for b=1:B
  for k=1:N-b
    A(k , k+b) = R(b);
    A(k+b,k)   = R(b);
  end
end
[~, D] = eig(A);
offset = min(D(:));
if(offset <= 0 )
  A = A + eye(N,N)*(0.5 - offset);
end
A = sparse(A);
end