function [x err iter flag convergence msg] = pgs(A, b, x0, max_iter, tol_rel, tol_abs, profile)
% Copyright 2011, Kenny Erleben, DIKU

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

N    = length(b); % Number of variables
flag = 1;

%--- Setup meaningfull default values -------------------------------------
if nargin<3
  x0 = zeros(N,1);
end
if nargin<4
  max_iter = floor(N/2);
end
if nargin<5
  tol_rel = 0.0001;
end
if nargin<6
  tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if nargin<7
  profile = true;
end

%--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);

%--- Setup values needed while iterating ------------------------------------
convergence = zeros(max_iter,1); % Used when profiling to measure the convergence rate
err         = Inf;               % Current error measure
x           = x0;                % Current iterate
iter        = 1;                 % Iteration counter
flag        = 2;

while iter <= max_iter
  
  dx = 0;
  for i=1:N
    old_xi = x(i);
    ri     = b(i) + A(i,:)*x;
    Aii    = A(i,i);
    
    x(i) = max( 0, old_xi - (ri / Aii) );
    dx = max(dx, abs(x(i) - old_xi));
  end
  
  old_err = err;
    
  y   = abs( A*x + b );   % Abs is used to fix problem with negative y-values.
  err = x'*y;
  
  if profile
    convergence(iter) = err;
  end
  
  % Relative stopping criteria
  if (abs(err - old_err) / abs(old_err)) < tol_rel  
    flag = 3;
    break;
  end
  
  % Absolute stopping criteria
  if err < tol_abs   
    flag = 4;
    break;
  end
  
  % Stagnation testing
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
