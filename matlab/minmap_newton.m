function [ x err  iter flag convergence msg] = minmap_newton( A, b, x0, max_iter, tol_rel, tol_abs, profile )
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

%--- Make sure we got good working default values -------------------------
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

%--- Here comes a bunch of magic constants --------------------------------

h       = 1e-7;    % Fixed constant used to evaluate the directional detivative
alpha   = 0.5;     % Step reduction parameter for projected Armijo backtracking line search
beta    = 0.001;   % Sufficent decrease parameter for projected Armijo backtracking line search
gamma   = 1e-28;   % Perturbation values used to fix near singular points in derivative
rho     = eps;     % Descent direction test parameter used to test if the Newton direction does a good enough job.

%--- Setup values need while iterating ------------------------------------

convergence = zeros(max_iter,1); % Used when profiling to measure the convergence rate

err     = Inf;         % Current error measure
x       = x0;          % Current iterate
iter    = 1;           % Iteration counter

flag    = 2;

while (iter <= max_iter )
  
  y = A*x + b;
  
  %--- Test all stopping criteria used ------------------------------------
  phi     = minmap(y,x);         % Calculate min map
  old_err = err;
  err     = 0.5*(phi'*phi);      % Natural merit function
  
  if profile
    convergence(iter) = err;
  end
  if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
    flag = 3;
    break;
  end
  if err < tol_abs   % Absolute stopping criteria
    flag = 4;
    break;
  end
  
  %--- Solve the Newton system --------------------------------------------  
  S       = find( y < x);
  J       = eye(N,N);
  J(S,:)  = A(S,:);
  
  % The more accurate a solution for the Newton direction the better  
  %dx = - J \ phi;

  % However, for large system it might be better to settle with an
  % approximate solution. Risk is that quadratic convergence drops to
  % superlinear, or even linear if accuracy is really bad
  restart = min( N, 20);
  [dx ~ ] = gmres( J, (-phi), restart);   % Use shur method -- its more efficient, I'm lazy here
  
  % If Shur method is used then reduced Newton system will have same matrix
  % properties as A. If A is symmmetric PD then PCG can be used on reduced
  % system. This will be much faster. If A is symmetric PSD then one must
  % find a preconditioner for PCG to work on the reduced A.
  
  % Test if the search direction is smaller than numerical precision. That is if it is too close to zero. 
  if max(abs(dx)) < eps
    flag = 5;
    % Rather than just giving up we may just use the gradient direction
    % instead. However, I am lazy here!
    break;
  end
  
  % Test if we have dropped into a local minimia if so we are stuck
  nabla_phi = phi'*J;
  if norm(nabla_phi) < tol_abs
    flag = 6;
    break;
  end
    
  % Test if our search direction is a 'sufficient' descent direction
  if  nabla_phi*dx  > -rho*(dx'*dx)
    flag = 7;
    % Rather than just giving up we may just use the gradient direction
    % instead. However, I am lazy here!
    % dx = nabla_phi'
    break;
  end
  
  %--- Armijo backtracking combined with a projected line-search ---------
  tau     = 1.0;                  % Current step length
  f_0     = err;
  grad_f  = beta*(nabla_phi*dx);
  x_k     = x;
  
  while true
    
    x_k   = max(0,x + dx*tau); 
    y_k   = A*x_k + b;
    phi_k = minmap( y_k, x_k );
    f_k   = 0.5*(phi_k'*phi_k);
    
    % Perform Armijo codition to see if we got a sufficient decrease
    if ( f_k <= f_0 + tau*grad_f)
      break;
    end
    
    % Test if time-step became too small
    if tau*tau < gamma
      break;
    end
    
    tau = alpha * tau;
  end
  
  % Update iterate with result from Armijo backtracking
  x = x_k;
  
  % Increment the number of iterations
  iter=iter+1;
end

if iter>=max_iter
  flag = 8;
  iter = iter - 1;
end

msg = msg{flag};
  
end

function [ phi ] = minmap(y,x)
% Auxiliary function used by the Minimum map Newton method
% Copyright 2011, Kenny Erleben, DIKU
phi  = min(y,x);
end
