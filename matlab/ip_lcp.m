function [ x err iter flag convergence msg ] = ip_lcp( A, b, x0, max_iter, tol_rel, tol_abs, solver, profile)
% Copyright 2011, 2012, Kenny Erleben, DIKU
% Copyright 2012,       Michael Andersen, DIKU
%
% TODO List:
%   - Add merit/potential reduction method (Solver, currently not used)
%   - Add a Params struct, containing parameters for subsolvers.
%
%  Given  A,b then the LCP is defined by
%
%         A x + b >= 0
%               x >= 0
%     x^T(A x + b) = 0
%
%  We rewrite this using the slack variable y = Ax + b as
%
%       (Ax - y + b)_i = 0
%              x_i y_i = 0
%            x_i, y_i >= 0
%
%  Idea is to solve this iteratively while keeping x_i^k, y_i^k > 0. To
%  ensure invariant we introduce a postive value gamma such that
%
%       (Ax - y + b)_i = 0
%      x_i y_i - gamma = 0
%             x_i, y_i > 0
%
%  Here gamma is a relaxed complimentarity measure (more on this below).
%
%  Now we introduce the Kojima mapping which we use to reformulate the
%  relaxed LCP formulation above
%
%
%                        |  A x - y + b        |
%        F(x,y:gamma) =  |  M W e  - gamma e   | = 0
%
%  where M and W are diagonal matrices with M_ii = x_i and W_ii = y_i and
%  e is the nD vector of ones. In the k'th iteration we pick gamma to be
%
%          gamma^k = sigma (x^ k)^T (y^k)/n
%
%  where n is the dimenstion of the problem and the relaxation (centering)
% parameter is 0 < sigma < 1.
%
%  In the k'th iteration of our framework we solve F(x,y: gamma^k)=
%  0 approximately. The solutions of the Kojima map for all positive
%  values of gamma define the central path which is a trajectory that
%  leads to the solution of the LCP as gamma tends to zero.
%
%  Using Taylor series expansion we define the Newton equation
%
%  |  A    -I  | | dx^k |   | - A x^k + y^k - b     |
%  |           | |      | = |                       |
%  |  W^k  M^k | | dy^k |   | - M^k W^k e + gamma e |
%
% Having solved for the Newton direction ( dx^k, dy^k ) we do the Newton update
%
%     x^(k+1)  = x^k  + tau^k dx^k
%     y^(k+1)  = y^k  + tau^k dy^k
%
% where tau^k is the step length. The step length is chosen as the
% maximum positive value such that
%
%    x^(k+1)_i > 0 and y^(k+1)_i > 0 for all i    (**)
%
% That is we choose tau^k = min( tau^x, tau^y ) where
%
%  tau^x = max{ 0 < tau  < 1 |  x^k + tau dx^k >= (1-alpha) x^k }
%  tau^y = max{ 0 < tau  < 1 |  y^k + tau dy^k >= (1-alpha) y^k }
%
% Here 0 < alpha < 1 controls how far we back away from the maximum step
% for which the conditions (**) are statisfied.
%
% To accelerate the convergence we make alpha approach one as the iterates
% approach the solution

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

N = length(b);

%--- Make sure we got good working default values -------------------------
if nargin < 3
    % This might not be that good a value. Might not be feasible at all.
    x0 = zeros(N,1);
end
if nargin < 4
    max_iter = floor(N/2);
end
if nargin < 5
    tol_rel = 0.0001;
end
if nargin < 6
    tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
end
if nargin<7
    solver = 'central trajectory';
end
if nargin < 8
    profile = true;
end

y0 = A * x0 + b;
e  = ones(N,1);
I  = eye(N,N);

% Make sure we have a feasible starting iterate x0>0 and y0>0
if x0'*y0 <= 0,
    W  = diag(y0);
    M  = diag(x0);
    J  = [ A, -I; W, M];
    Fk = [ A*x0 - y0 + b;  M*W*e - 1 ];
    
    % Help the stability, due to badly conditioned J
    if rcond(J) < 2*eps
        p = pinv(J)*Fk;
    else
%         p  = - inv(J) * Fk;
        p = -J\Fk;
    end
    dx = p(1:N);
    dy = p(N+1:end);
    
    x0 = max(1, abs( x0 + dx ) );
    y0 = max(1, abs( y0 + dy ) );
end

err = inf;
x = x0;
y = y0;
iter = 1;
convergence = zeros(max_iter,1);

flag = 2;

while iter <= max_iter,
    
    old_err = err;
    err = x'*y;
    
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
    
    if (x'*y <= 0)
        display([ 'ERROR  x''*y  = ', num2str(x'*y) ] );
        return;
    end
    
    alpha = iter/ (max_iter+1);    % Goes to one as we iterate
    sigma =  1 ./ iter;            % Goes to zero as we iterate
    gamma = sigma * (x'*y) ./ N;   % Complimentarity measure
    
    %--- Setup Newton system --------------------------------------------
    W  = diag(y);
    M  = diag(x);
    J  = [ A, -I; W, M];
    Fk = [ A*x - y + b;  M*W*e - gamma*e ];
    
    %--- Solve Newton system --------------------------------------------
%     p  = - inv(J) * Fk;        % 2012-01-08 Kenny: Can we do this linear system solve faster? Maybe use some iterative solver?
    p = -J\Fk;                   % 2012-06-18 Michael: This is faster than 
                                 % the above. We might consider inserting a 
                                 % rcond control to see the conditioning
                                 % and then use pinv instead.
    % 2012-01-09 Kenny: Iterative solutions seems to destroy accuracy very
    % fast -- IP really needs a accurate solution for the Newton system.
    % tol = max( 10e-10, eps*norm(Fk) );
    % max_iter = n/2;
    % restart = max(10, n/10);
    % [p ~] = gmres(J, -Fk, restart, tol, max_iter);
    
    dx = p(1:N);
    dy = p(N+1:end);
    
    %--- Determine equal step lengths ------------------------------------
    xx = ((1-alpha)*x  - x) ./ dx;
    yy = ((1-alpha)*y  - y) ./ dy;
    
    if(~isempty(xx(dx < 0)))
        taux =  min( 1, min ( xx(dx < 0) ) );
    else
        taux = 1;
    end
    if(~isempty(yy(dy < 0)))
        tauy =  min( 1, min ( yy(dy < 0) ) );
    else
        tauy = 1;
    end
    tau  =  min(taux, tauy);
    
    if ( (tau<0) || (tau > 1) )
        display([ 'ERROR  Tau   = ', num2str(tau) ] );
        return;
    end
        
    %--- Do the Newton update --------------------------------------------
    x = x + taux*dx;  % 2012-01-09 Kenny: Split update seems to give faster convergence rate! At least for dense sym PSD LCPs.
    y = y + tauy*dy;
    
    if min(x(:)) <= 0
        display([ 'ERROR  Min x = ', num2str(min(x(:))) ] );
        return;
    end
    if min(y(:)) <= 0
        display([ 'ERROR  Min y = ', num2str(min(y(:))) ] );
        return;
    end
    
    %--- Increment iteration counter -------------------------------------
    iter = iter + 1;
end

if iter >= max_iter
    flag = 8;
    iter = iter - 1;
end

% Only track convergence as far as we have taken iterates
convergence = convergence(1:iter);
% Just return the message, not the cell of messages
msg = msg{flag};

end

% Some notes on how to solve the Newton system:
%
% The Newton system is given as
%
%  | A -I  |  | x | = | v |
%  | W  M  |  | y |   | w |
%
% From the last row we find
%
%     y = inv(M)( w - W x)
%
% Substitution into first row gives the reduced system
%
%   (A + inv(M) W) x = v + inv(M) * w
%
% which we can solve for x.
%
%  Observe if A is sym PSD then A + inv(M)*W is sym PD and the reduced
%  system can be solved using PCG method.
%
%  If A is non-sym then A+inv(M)W is also unsym and we may use GMRES to
%  find a solution for the reduced system.
%
%  In general we do not need A to be non-singular as long as A + inv(W)*M
%  is non-singular we can find an unique solution for the Newton system.
%