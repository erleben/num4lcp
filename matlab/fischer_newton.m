function [ x err iter flag convergence msg] = fischer_newton( A, b, x0, max_iter, tol_rel, tol_abs, solver, profile )
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
  solver = 'random';
end
if nargin<8
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

flag = 2;

while (iter <= max_iter )
  
  y = A*x + b;
  
  %--- Test all stopping criteria used ------------------------------------
  phi     = fischer(y,x);         % Calculate fischer function
  old_err = err;
  err     = 0.5*(phi'*phi);       % Natural merit function
  
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
  %--- First we find the Jacobian matrix. We wish to avoid singular points
  %--- in Jacobian of the Fischer function
  dx = zeros(N,1);                     % Allocate space for computing the Newton direction
  
  S       = abs(phi)<gamma & abs(x)<gamma;  % Bitmask for singular indices
  I       = find(S==0);                     % Bitmask for non-singular indices
  restart = min(size(A,1),10);              % The number of iterations done before GMRES should restart
  
  switch lower(solver)
    
    case 'random'    % works on full system

      q       = rand(N,1) ./ sqrt(2)  - 1;
      p       = rand(N,1) ./ sqrt(2)  - 1;
      J       = sparse( zeros(N,N) );   
      q(I)    = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
      p(I)    = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
      J(I,I)  = diag(p(I))*eye(length(I)) + diag(q(I))*A(I,I);
      
      [dx ~]  = gmres( J, (-phi), restart);
      
    case 'perturbation'     % works on full system
      
      px          = x;
      dir         = sign( x );
      dir(dir==0) = 1;
      px(S==1)    = gamma * dir(S==1);
      
      p = (px./((y.^2+px.^2).^0.5))-1;
      q = ( y./((y.^2+px.^2).^0.5))-1;
      J = diag(p)*eye(N) + diag(q)*A;
      
      % Billups solution actually pick the gradient direction as the descent
      % direction if the Newton system can not be solved. Look in his
      % thesis at page 85.
      %
      % %if rcond(J) < eps
      % %  dx = -J'*phi;          % Too bad we pick gradient direction
      % %end
      %
      %Instead one may use the pseudo inverse if badly conditioned
      %
      if rcond(J) < eps*2
        dx = pinv(J) * (-phi); % badly conditioned
      else
        dx = J \ (-phi);       % well conditioned
      end
      %
      % Or one can hope for the best and try GMRES
      %
      %[dx ~]  = gmres( J, (-phi), restart);
      %
      
    case 'zero'   % works on reduced system would be similar to random-case if all radom values were set to 1

      J       = sparse( zeros(N,N) );
      q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
      p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
      J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);
      
      [dx(I) ~] = gmres( J(I,I), (-phi(I)), restart);

    case 'approximation'  % works on reduced system assuming zero Jacobian for all bad values

      J       = sparse( zeros(N,N) );
      q       = (y(I)./((y(I).^2+x(I).^2).^0.5))-1;
      p       = (x(I)./((y(I).^2+x(I).^2).^0.5))-1;
      J(I,I)  = diag(p)*eye(length(I)) + diag(q)*A(I,I);
      
      fun       = @(dx) ( fischer( A(I,I)*(x(I)+ dx*h) + b(I), x(I) + dx*h) - phi(I) ) / h;
      [dx(I) ~] = gmres( fun, (-phi(I)), restart);

    otherwise
      disp('Unknown solver method for Newton subsystem')
  end
  
  % Test if the search direction is smaller than numerical precision. That is if it is too close to zero. 
  if max(abs(dx)) < eps
    flag = 5;
    % Rather than just giving up we may just use the gradient direction
    % instead. However, I am lazy here!
    %  dx = nabla_phi'
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
    %  dx = nabla_phi'
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
    phi_k = fischer( y_k, x_k );
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
  
end

function [ phi ] = fischer(y,x)
% Auxiliary function used by the Fischer-Newton method
% Copyright 2011, Kenny Erleben, DIKU
phi  = (y.^2 + x.^2).^0.5 - y - x;
end
