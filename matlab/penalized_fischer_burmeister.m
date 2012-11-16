function [ x err iter flag convergence msg ] = penalized_fischer_burmeister( A, b, x0, max_iter,tol_rel, tol_abs, profile )
% Copyright 2012, Michael Andersen, DIKU
% Copyright 2011, Kenny Erleben, DIKU
% Based on fischer_newton.m from num4lcp (http://code.google.com/p/num4lcp/)

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
    % sigma \in (0,1/2)
    sigma = 0.001;
    % p > 2
    p = 2.1;
    % rho > 0
    rho = 1e-10;
    % lambda
    l = 0.95;
    
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
    
    convergence = [];
    err = 1e20;
    iter = 0;
    
    flag = 2;
    while iter < max_iter

        y = A*x + b;
        
        %--- Test all stopping criteria used ------------------------------
        phi_l = phi_lambda(x,y,l);
        old_err = err;
        err = 0.5*(phi_l'*phi_l);
        
        fprintf('err: %2.15f\n', err)
        
        if profile
            convergence = [convergence err];
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

        % Check conditioning of Vk
        if rcond(Vk) < 2*eps
            dx = pinv(Vk) * (-phi_l);
        else
            dx = Vk \ (-phi_l);
        end
    
        if max(abs(dx)) < eps
            flag = 5;
            break;
        end
        
        nabla_psi = Vk'*phi_l;
        if norm(nabla_psi) < tol_abs
            flag = 6;
            break;
        end
                
        if nabla_psi'*dx > -rho*(dx'*dx)^p
            flag = 7;
            break;
        end
        
        %--- Armijo backtracking combined with a projected line-search ----
%         tau     = 1.0;                  % Current step length
%         f_0     = err;
%         grad_f  = sigma*(nabla_psi'*dx);
%         x_k     = x;
%         
%         while true
%             
%             x_k   = max(0,x + dx*tau);
%             y_k   = A*x_k + b;
%             f_k = psi_lambda( x_k, y_k, l );
%             
%             % Perform Armijo codition to see if we got a sufficient decrease
%             if ( f_k <= f_0 + tau*grad_f)
%                 break;
%             end
%             
%             % Test if step-size became too small
%             if tau*tau < gamma
%                 break;
%             end
%             
%             tau = alpha * tau;
%         end
        
        lk = 0;
        x_k = x;
        f_0 = err;
        grad_f = sigma*(nabla_psi);
        while true
            % This makes it converge, but it is wrong
            % x_k = max(0,x_k + beta^lk*dx);
            % Correct, but step size becomes in the order of 1e-15 :(
            x_k = max(0,x + beta^lk*dx);
            y_k = A*x_k + b;
            f_k = psi_lambda(x_k, y_k, l);
            if f_k <= f_0 + beta^lk*grad_f
                break;
            end
            if beta^lk < eps
                break;
            end
            lk = lk + 1;
        end

        % Update iterate with result from Armijo backtracking
        x = x_k;
        
        % Increment iterations
        iter = iter + 1;
    end
    msg = msg{flag};
end

function V = produceV(x, y, A,l)
% Calculates V \in R^{n \times n} such that V is an element of
% C-subdifferential \partial_{C}\Phi_{\lambda}(x) 
%%% Following Algorithm 2.4

n = size(x,1);
%%% Generate S1 = {i | x_i == 0, F_i(x) == 0} | MASK
S1 = (x == 0) & (y == 0);
% S1 = (x < gamma) & (y < gamma);
%%% Generate S2 = {i | x_i > 0, F_i(x) > 0} | MASK
S2 = (x > 0) & (y > 0);
% S2 = (x >= gamma) & (y >= gamma);
%%% Generate i \notin S1 U S2 | MASK
NS1US2 = ~(S1 | S2);

z = zeros(n, 1);
%%% z_i = 1 if i \in S1, z_i = 0 if i \notin S1
z(S1) = ones(size(find(S1),1),1);

V = zeros(n,n);

I_S1 = find(S1);
V(I_S1,I_S1) = diag(l*(1 - z(I_S1)./(z(I_S1).^2+(A(I_S1,I_S1)*z(I_S1)).^2).^0.5))*eye(length(I_S1)) + ...
        diag(l*(1 - (A(I_S1,I_S1)*z(I_S1))./(z(I_S1).^2+(A(I_S1,I_S1)*z(I_S1)).^2).^0.5))*A(I_S1,I_S1);

I_S2 = find(S2);
V(I_S2,I_S2) = diag((l*(1 - x(I_S2)./(x(I_S2).^2+y(I_S2).^2).^0.5+(1-l)*y(I_S2))))*eye(length(I_S2)) + ...
        diag((l*(1 - y(I_S2)./(x(I_S2).^2+y(I_S2).^2).^0.5+(1-l)*x(I_S2))))*A(I_S2,I_S2);

I_N = find(NS1US2);
V(I_N,I_N) = diag(l*(1 - x(I_N)./(x(I_N).^2+y(I_N).^2).^0.5))*eye(length(I_N)) + ...
        diag(l*(1 - y(I_N)./(x(I_N).^2+y(I_N).^2).^0.5))*A(I_N,I_N);
end

% function V = produceV(x, y, A,l)
% % Calculates V \in R^{n \times n} such that V is an element of
% % C-subdifferential \partial_{C}\Phi_{\lambda}(x) 
% %%% Following Algorithm 2.4
% 
% n = size(x,1);
% %%% Generate S1 = {i | x_i == 0, F_i(x) == 0} | MASK
% S1 = (x == 0) & (y == 0);
% % S1 = (x < gamma) & (y < gamma);
% %%% Generate S2 = {i | x_i > 0, F_i(x) > 0} | MASK
% S2 = (x > 0) & (y > 0);
% % S2 = (x >= gamma) & (y >= gamma);
% %%% Generate i \notin S1 U S2 | MASK
% NS1US2 = ~(S1 | S2);
% 
% z = zeros(n, 1);
% %%% z_i = 1 if i \in S1, z_i = 0 if i \notin S1
% z(S1) = ones(size(find(S1),1),1);
% 
% e = ones(n,1);
% I = eye(n,n);
% V = zeros(n,n);
% I_S1 = find(S1);
% for i = 1 : length(I_S1)
%     V(I_S1(i),:) = l*(1 - z(I_S1(i))/norm([z(I_S1(i)) A(I_S1(i),:)*z]))*I(I_S1(i),:) + ...
%         l*(1 - (A(I_S1(i),:)*z)/norm([z(I_S1(i)) A(I_S1(i),:)*z]))*A(I_S1(i),:);
% end
% I_S2 = find(S2);
% for j = 1 : length(I_S2)
%     V(I_S2(j),:) = (l*(1 - x(I_S2(j))/norm([x(I_S2(j)) y(I_S2(j))]))+(1-l)*y(I_S2(j)))*I(I_S2(j),:) + ...
%         (l*(1 - y(I_S2(j))/norm([x(I_S2(j)) y(I_S2(j))]))+(1-l)*x(I_S2(j)))*A(I_S2(j),:);
% end
% I_N = find(NS1US2);
% for k = 1 : length(I_N)
%     V(I_N(k),:) = l*(1 - x(I_N(k))/norm([x(I_N(k)) y(I_N(k))]))*I(I_N(k),:) + ...
%         l*(1 - y(I_N(k))/norm([x(I_N(k)) y(I_N(k))]))*A(I_N(k),:);
% end
% end

function [ phi ] = phi_lambda(a, b, l)

phi = l*fischer(a, b) + (1 - l).*phi_plus(a, b);

end

function [ phi ] = fischer(y,x)
% Auxiliary function used by the Fischer-Newton method
% Copyright 2011, Kenny Erleben, DIKU
phi  = (y.^2 + x.^2).^0.5 - y - x;
end

function [ phi ] = phi_plus(a, b)
phi = max(a,0).*max(b,0);
end

function [psi] = psi_lambda(x, y, l)
psi = 0.5*(phi_lambda(x,y,l)'*phi_lambda(x,y,l));
end