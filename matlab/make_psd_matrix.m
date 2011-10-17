function [A] = make_psd_matrix(N, ratio)
% [A] = MAKE_PSD_MATRIX(N, ratio) Generates a random symmetric positive
% semi definite matrix with a desired eigenvalue spectrum.
%
% Input:
%
%    N    -- Number of variables.
%   ratio -- The fraction of zero eigenvalues.
%
% Output:
%
%    A    -- A symmetric positive semi definite matrix.
%
% Copyright 2011, Kenny Erleben, DIKU

% A = make_pd_matrix(N);
% [V D] = eig(A);
% K = floor( N*ratio );
% if K>0
%   D(1:K,1:K) = 0;
% end
% A = V * D * V';
% A = 0.5*(A + A');

% Faster method that does not need eig routine
d = (1:1:N) ./ N;
K = floor( N*ratio );
if K>0
  d(1:K) = 0;
end
R = rand(N,N);
O = orth(R);
A = O*diag(d)*O';
A = sparse( 0.5*(A+A') );
end