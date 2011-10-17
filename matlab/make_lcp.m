function [x, b] = make_lcp(A,F)
% [x, b] = MAKE_LCP(A): Make LCP
%
% INPUT:
%
%   A -- The coefficient matrix in the LCP
%   F -- The fraction of zero-values in the x-solution.
%
% OUTPUT:
%
%   x -- A solution for the LCP problem.
%   b -- The right hand side vector
%
% Copyright 2011, Kenny Erleben, DIKU

%--- Get number of variables ----------------------------------------------
N = length(A);

%--- Generate a random LCP solution ---------------------------------------
x         = rand(N,1);
x(x < F ) = 0;

%--- Generate a right-hand-side vector that is ----------------------------
%--- consistent with random solution           ----------------------------
b = zeros(N,1);
s = A * x;
b(x>0) = - s(x>0);
end