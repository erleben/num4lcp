function [P, t_coa, t_int] = create_prolongation(A)
%
% input:
%
%   A       -- The coefficient matrix of the fine level.
%
% output:
%
%    P      -- Prolongation opreator, from coarse to fine level
%              interpolation.
%    t_coa  -- Time in seconds that it took to compute coarsening.
%    t_int  -- Time in seconds that it took to compute interpolator.
%
% Copyright 2011, Kenny Erleben, DIKU
  gamma   = 0.25;
  tic
  [C F S] = mis_coarsening(A, gamma);
  t_coa = toc;
  
  tic
  P = interpolator( A, C, F, S);
  t_int = toc;
end
