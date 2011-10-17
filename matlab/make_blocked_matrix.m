function A = make_blocked_matrix(B,K,F)
% A = MAKE_BLOCKED_MATRIX(B,K,F): Make a random blocked matrix.
%
% INPUT:
%
%   B -- Block size
%   K -- Number of blocks
%   F -- Fillin protentage (number between 0 and 1)
%
% OUTPUT:
%
%   A -- A blocked symmetric positive definite matrix
%
% Copyright 2011, Kenny Erleben, DIKU
N    = B*K;              % Number of variables
R    = rand(N,N);
A    = 0.5 * (R + R');  % Symmetric random matrix values between 0 and 1
for i=1:K
  i_from   = (i-1)*B + 1;
  i_center = i_from + floor(B/2);
  i_to     = i_from + B-1;
  for j=1:K
    
    j_from   = (j-1)*B + 1;
    j_center = j_from + floor(B/2);
    j_to     = j_from + B-1;
    
    if A(i_center, j_center) > F,
      A(i_from:i_to,j_from:j_to) = zeros(B,B);
    end
  end
end
[~, D] = eig(A);
offset = min(D(:));
if(offset <= 0 )
  A = A + eye(N,N)*(0.5 - offset);
end
end