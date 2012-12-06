function [x, err, iter, flag, convergence, msg, t_setup] = multilevel(A, b, x0, k, max_iter, tol_rel, tol_abs, profile)
% Copyright 2011, Kenny Erleben, DIKU

addpath('aux'); % Make sure that all subroutines can be found.

% Just a list of human readable text strings to convert the flag return
% code into something readable by writing msg(flag) onto the screen.
msg = {'preprocessing';  % flag = 1
  'iterating';      % flag = 2
  'relative';       % flag = 3
  'absolute';       % flag = 4
  'stagnation';     % flag = 5
  'local minima';   % flag = 6
  'nondescent';     % flag = 7
  'maxlimit'        % flag = 8
  };

if nargin < 2
  error('Too few arguments');
end
flag = 1;

%--- Setup meaningfull default values -------------------------------------
if nargin<3
  x0 = zeros(size(b));
end
if nargin<4
  k = 2;
end
if nargin<5
  max_iter = floor( length(b) /2);
end
if nargin<6
  tol_rel = 0.0001;
end
if nargin<7
  tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if nargin<8
  profile = true;
end

%--- Make sure all values are valid ---------------------------------------
k        = max(k,2);
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);

%--- Setup values needed while iterating ------------------------------------
convergence = zeros(max_iter,1); % Used when profiling to measure the convergence rate
err         = Inf;               % Current error measure
x           = x0;                % Current iterate
c           = zeros(size(x));    % Constraint vector
iter        = 1;                 % Iteration counter

t_start     = tic;
[As, Ps]    = setup(k, A, false); % Setup multilevel information
t_setup     = toc(t_start);

flag        = 2;

while iter <= max_iter
  
  old_x = x;
  x = vcycle( k, b, x, c, As, Ps);
  
  % Compute error/metrit value
  old_err = err;
  
  err = qp_merit(A,b,x);
  
  if profile
    convergence(iter) = err;
  end
  
  % Relative stopping criteria
  if (abs(err - old_err) / abs(old_err)) < tol_rel
    flag = 3;
    break;
  end
  
  % Absolute stopping criteria.
  % Even though we know that the quadratic
  % function is bounded form below we do not know the value of the bound.
  % Therefore the quadratic-merit function is a poor choice for absolute
  % convergence test. Instead we use the same test as in the PGS method.
  y   = abs( A*x + b );
  tst = x'*y;
  if tst < tol_abs
    flag = 4;
    break;
  end
  
  % Stagnation testing
  dx = max( abs( x - old_x ) );
  if dx < eps
    flag = 5;
    break;
  end
  
  iter = iter + 1;
end
if iter>=max_iter
  flag = 8;
  iter = iter - 1;
end
msg = msg{flag};
end