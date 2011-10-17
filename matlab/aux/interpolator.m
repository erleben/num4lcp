function [ P ] = interpolator( A, C, F, S)
% Copyright 2011, Kenny Erleben, DIKU

%--- Get dimensions of problem -------------------------------------------
M = length(A);
N = sum(C(:));

%--- Define interpolatinon weights ---------------------------------------
P            = zeros( M, N);
P( C==1, : ) = eye( N, N);
P( F==1, : ) = 0; %--- This we need to fix

Nj  = A ~= 0;
% Fj  = Nj - S;
% Cj  = repmat( C', M, 1) .* S;
% Fwj = repmat( F', M, 1) .* Fj;
% Fsj = repmat( C', M, 1) .* Fj;
w = 1 ./ (Nj * C);
W = repmat(w,1,M) .* Nj;
P( F==1, : ) = W( F==1, C==1 );

P = sparse(P);

end
