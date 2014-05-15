function [ As, Ps, data ] = setup( K, A, profiling )
% [ As, Ps ] = SETUP( K, A ): Multilevel setup routine.
%
% input:
%          K -- The number of levels that should be created
%          A -- The initial fine level (K-level) matrix.
%  profiling -- boolean variable indicating if profiling is on or off.
%
% output:
%
%   As   -- A K-by-1 cell array. Entry A{k} holds the k-level A-matrix
%   Ps   -- A K-by-1 cell array. Entry P{k} holds the (k-1) to k level
%           prolongation matrix.
%   data -- If profiling is one then this return value contains info about
%           the setup progress.
%
% Copyright 2011, Kenny Erleben, DIKU
if nargin < 2
  error('Too few arguments');
end
if nargin<3
  profiling = false;
end

Ps = cell(K,1);
As = cell(K,1);
As{K} = sparse(A);

data = [];
for k=K:-1:2
  
  tic;
  [P, t_coa, t_int]= create_prolongation(A);
  t_pro = toc;
    
  tic;
  A = sparse( P' * A * P );
  t_mul = toc;
  fillin = nnz(A)/(length(A).^2);

  Ps{k}  = P;
  As{k-1}  = A;

  if profiling
    data = [data; [t_pro, t_coa, t_int, t_mul, fillin] ];
  end
end
end