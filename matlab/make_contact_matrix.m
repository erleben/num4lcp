function [A, x, b] = make_contact_matrix(K)
% [A, x, b] = MAKE_CONTACT_MATRIX(K): Generates a contact matrix that has
% the same fill-in properties as that of the Stewart-Trinkle contact
% formulation using a 2D friction cone.
%
% INPUT:
%
%   K -- Number of contacts.
%
% OUTPUT:
%
%   A -- The contact coefficient matrix.
%
% Copyright 2011, Kenny Erleben, DIKU

%--- Generate a random blocked matrix -------------------------------------
B = 6;      % Block size, one normal, four friction, one slack
F = 0.5;    % Average procentage fill-in
N = B*K;    % Number of variables
A = make_blocked_matrix(B,K,F);

%--- Modify random matrix to have pattern of contact matrix ---------------
for bi=1:K
for bj=1:K
  i = (bi-1)*B + 1;
  j = (bj-1)*B + 1;
  if A(i,j) ~= 0
    A(i+5,j)       = 0.5;  % Set mu (friction) block
    A(i+5,j+5)     = 0;    % set zero diagonal block
    A(i,j+5)       = 0;    % set upper zero block
    A(i+5,j+1:j+4) = -1;   % set -e^T block
    A(i+1:i+4,j+5) = 1;    % set  e   block
  end
end
end
end