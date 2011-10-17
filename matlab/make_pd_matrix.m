function [A] = make_pd_matrix(N)
% [A] = MAKE_PD_MATRIX(N) Generates a random symmetric positive
% definite matrix.
%
% Input:
%
%    N    -- Number of variables.
%
% Output:
%
%    A    -- A symmetric positive definite matrix.
%
% Copyright 2011, Kenny Erleben, DIKU


% Fast but no control over spectrum
%R    = rand(N,N);
%A = R'*R + 0.1*eye(N,N);

% Slower but tries to keep values bounded
%R    = rand(N,N);
%A    = 0.5 * (R + R'); 
%[~, D] = eig(A);
%offset = min(D(:));
%if(offset <= 0 )
%  A = A + eye(N,N)*(0.1 - offset);
%end

% Fast method and creates a liner spectrum of eigenvalues
d = (1:1:N) ./ N;
R    = rand(N,N);
O = orth(R);
A = O*diag(d)*O';
A = sparse( 0.5*(A+A') );
end