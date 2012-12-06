function [ x err iter flag convergence msg ] = penalized_fischer_burmeister( A, b, x0, max_iter,tol_rel, tol_abs, profile )
% Copyright 2012, Michael Andersen, DIKU
% Copyright 2011, Kenny Erleben, DIKU
% Based on fischer_newton.m from num4lcp (http://code.google.com/p/num4lcp/)

% addpath('../../kanzow/Kanzow/');

% Based on the paper: A Penalized Fischer-Burmeister NCP-Function:
% Theoretical Investigation And Numerical Results


%%% Human readable exit flag
    msg = {'preprocessing';  % flag = 1
           'iterating';      % flag = 2
           'relative';       % flag = 3
           'absolute';       % flag = 4
           'stagnation';     % flag = 5
           'local minima';   % flag = 6
           'nondescent';     % flag = 7
           'maxlimit'        % flag = 8
           };

    flag = 1;

    %%% Constants copied from the paper.
    % beta \in (0,1)
    beta = 0.5;
    alpha = 0.5;
    % sigma \in (0,1/2)
    sigma = 0.001;
    % p > 2
    p = 2.1;
    % rho > 0
    rho = 1e-10;
    % lambda \in [0,1]
    l = 0.95;
%     l = 1;
    
    [m,n] = size(A);
    
    if m ~= n
        error(['Matrix A is not square: ' num2str(m) 'x' num2str(n)]);
    end
    
    if nargin < 3
        x0 = ones(n,1);
    end
    if nargin < 4
        max_iter = floor(n/2);
    end
    if nargin < 5
        tol_rel = 0.00001;
    end
    if nargin < 6
        tol_abs = 10*eps;
    end
    if nargin < 7
        profile = false;
    end
    
    %--- Ensure that all values are valid ---------------------------------
    max_iter = max(max_iter,1);
    tol_rel  = max(tol_rel,0);
    tol_abs  = max(tol_abs,0);
    x0       = max(0,x0);
        
    x = x0;
        
    convergence = zeros(max_iter,1);
%     convergence = [];
    err = 1e20;
    iter = 1;
    
    flag = 2;
    while iter <= max_iter

        y = A*x + b;
        
        %--- Test all stopping criteria used ------------------------------
        phi_l = phi_lambda(x,y,l);
        old_err = err;
        err = 0.5*(phi_l'*phi_l);
        
%         fprintf('err: %2.15f\n', err)
        
        if profile
            %%% Required for the rpi-simulator to track err
            convergence(iter) = err;
%             convergence = [convergence err];
        end
        
        if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
            flag = 3;
            break;
        end        
        if err < tol_abs
            flag = 4;
            break;
        end
        
        Vk = produceV(x,y,A,l);
%         Vk2 = produceV2(x,y,A,l);
%         Vk = -Vk2;
%         pause
        if rcond(Vk) < eps*2
            dx = pinv(Vk) * (-phi_l); % badly conditioned
        else
            dx = Vk \ (-phi_l);       % well conditioned
        end

%         dx = Vk \ (-phi_l);

        nabla_psi = phi_l'*Vk;
        if max(abs(dx)) < eps
            flag = 5;
            break;
%             dx = -nabla_psi';
        end
                
        if norm(nabla_psi) < tol_abs
            flag = 6;
            break;
        end
                
        if nabla_psi*dx > -rho*(dx'*dx)^p
            flag = 7;
            break;
%             dx = -nabla_psi';
        end
         
%         size(nabla_psi)
%         size(dx)
        
        %--- Armijo backtracking combined with a projected line-search ----
        tau     = 1.0;                  % Current step length
        f_0     = err;
        grad_f  = sigma*(nabla_psi*dx);
        x_k     = x;
        gamma = 1e-16;
        while true
            
            x_k   = max(0,x + dx*tau);
            y_k   = A*x_k + b;
            f_k = psi_lambda( x_k, y_k, l );
            
            % Perform Armijo codition to see if we got a sufficient decrease
            if ( f_k <= f_0 + tau*grad_f)
                break;
            end
            
            % Test if step-size became too small
            if tau*tau < gamma
                break;
            end
            
            tau = alpha * tau;
        end
        
%         lk = 0;
%         x_k = x;
%         f_0 = err;
%         grad_f = sigma*(nabla_psi*dx);
%         while true
%             % This makes it converge, but it is wrong
%             % x_k = max(0,x_k + beta^lk*dx);
%             % Correct, but step size becomes in the order of 1e-15 :(
%             x_k = max(0,x + beta^lk*dx);
%             y_k = A*x_k + b;
%             f_k = psi_lambda(x_k, y_k, l);
%             if f_k <= f_0 + beta^lk*grad_f
%                 break;
%             end
%             if beta^lk < eps
%                 break;
%             end
%             lk = lk + 1;
%         end
%         beta^lk
        % Update iterate with result from Armijo backtracking
        x = x_k;
%         pause
        % Increment iterations
        iter = iter + 1;
    end
    
    if iter>=max_iter
        flag = 8;
        iter = iter - 1;
    end

    msg = msg{flag};
end

function Vk = produceV(x,y, A, l)
%% Produces Vk
% Write parameters and such.
%   
%
%

SMALL = 1e-10;
n = size(x,1);
z = zeros(n,1);
S1 = abs(x) <= SMALL & abs(y) <= SMALL;
S2 = x > SMALL & y > SMALL;
NS = ~(S1 | S2);
Da = zeros(n,1);
Db = zeros(n,1);
%%% Does not really seem to matter whether we use rand or ones!
z(S1) = ones(length(find(S1)),1);
% z(S1) = 0.5*ones(length(find(S1)),1);
% z(S1) = rand(length(find(S1)),1);

denomS1 = sqrt(z(S1).^2+(A(S1,S1)*z(S1)).^2);
Da(S1) = l*(z(S1)./denomS1-1);
Db(S1) = l*((A(S1,S1)*z(S1))./denomS1-1);

denomS2 = sqrt(x(S2).^2+y(S2).^2);
Da(S2) = l*(x(S2)./denomS2-1)+(1-l)*y(S2);
Db(S2) = l*(y(S2)./denomS2-1)+(1-l)*x(S2);

denomNS = sqrt(x(NS).^2+y(NS).^2);
Da(NS) = l*(x(NS)./denomNS-1);
Db(NS) = l*(y(NS)./denomNS-1);

Vk = diag(Da) + diag(Db)*A;

end

function Vk = produceV2(x, y, A, l)

    grad_ab = gradient_cck(x,y,l);
%     grad_ab = kanzow_gradient(x,y,l);
    Da = grad_ab(:,1);
    Db = grad_ab(:,2);
    Vk = diag(Da) + diag(Db)*A;

end

function phi_l = phi_lambda(a,b,lambda)
%% CCK NCP-function, a convex composition of Fischer-Burmeister and CCK NCP
%
%   Input
%       a -> A column vector size = (n,1)
%       b -> A column vector size = (n,1)
%       l -> A fixed lambda value used to weight the input.
%
%   Output
%       phi_l -> A column vector with the result of the Fischer-Burmeister
%                NCP function, with size = (n,1)

    phi_l = lambda*phi_fb(a,b)+(1-lambda)*(max(0,a).*max(0,b));
    
end

function phi = phi_fb(a, b)
%% Fischer-Burmeister NCP-function
%
%   Input
%       a -> A column vector size = (n,1)
%       b -> A column vector size = (n,1)
%
%   Output
%       phi -> A column vector with the result of the Fischer-Burmeister
%              NCP function, with size = (n,1)

    phi = sqrt(a.^2+b.^2) - (a + b);

end

function psi_l = psi_lambda(a,b,l)
%% Natural merit function of phi_lambda. psi_lambda : R^n -> R
%
%   Input:
%       a -> A column vector of size = (n,1)
%       b -> A column vector of size = (n,1)
%       l -> Real value describing weight of penalization term
%
%   Output:
%       psi_l -> Merit values at x^(k) based on phi_l(x^(k))
 
    phi_l = phi_lambda(a,b,l);
    psi_l = 0.5*(phi_l'*phi_l);

end