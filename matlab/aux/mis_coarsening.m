function [C, F, S] = mis_coarsening( A, gamma )
% [C, F, S] = MIS_COARSENING( A, gamma ) -- Use maximal independent set
% clustering to define a coarsening set.
% Input:
%
%   A     - The A-matrix that should be coarsen
%   gamma - A threshold value between 0 and 1 which indicates the size of
%           the resulting coarse set. The lower gamma value the more
%           aggressive coarsing.
%
% Output:
%   C      - Bitmask for coarse level variables. If C(i) == 1 the i should
%            be on the coarse level.
%
%   F      - Bitmask for the fine level variables. If F(i) == 1 then i
%            should be on the fine level.
%   S      - Influence matrix. S(i,:)==1 is the set of all nodes that
%            strongly influence node i. S(:,j)==1 is the set of all nodes
%            that strongly depend on node j
%
% Copyright 2011, Kenny Erleben, DIKU

%--- By construction we have that 
%---
%---  S(i,j) = 1 
%---
%--- If and only if 
%---
%---  abs( A(i,j) ) > gamma max{ abs(A(i,k)) for all k!=i }
%---
%--- Otherwise
%---
%---  S(i,j) = 0
%---
%--- Observe that
%---
%---   S(i,:)==1 is the set of all nodes that strongly influence node i
%---   S(:,j)==1 is the set of all nodes that strongly depend on node j
%---
B     = sparse( abs( A  - diag(diag(A)) ) );
test  = sparse( gamma*repmat( max( B, [], 1 )' ,1,length(A)) );
S     = sparse( B > test );

lambda = sum(S, 1 )';         % lambda(j) number of nodes that strongly depend on node j
U      = ones(length(A),1);   % All unprocessed nodes
C      = zeros(length(A),1);  % Coarse-level nodes
F      = zeros(length(A),1);  % Fine-level nodes

%--- Keep iterating while U is non-empty
while sum(U)>0, 
 
  [~, k]       = max( lambda );         % Pick unprocessed node k with max influende
  C(k)         = 1;                     % Put k in C
  U(k)         = 0;                     % Remove processed nodes 
  N            = sparse( S(:,k).*U );   % Get unprocessed nodes that depend on k
  F( N==1 )    = 1;                     % Put N in F
  U( N==1 )    = 0;                     % Remove processed nodes

  % Now increase chance neighbors of N will become C nodes
  T = sparse( repmat(N,1,length(A)) .* S );   % Filter: only use N rows
  M = sum(T,1) > 0;                 % Unprocessed neighbors of N
  lambda(M==1) = lambda(M==1) + 1;
  lambda(U==0) = -1;                    
end
S = sparse(S);
end