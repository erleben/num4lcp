function [ theta ] = qp_merit(A,b,x)
% Copyright 2011, Kenny Erleben, DIKU
theta = ((x'*A *x) ./ 2) + x'*b;
%y = abs(A *x + b);
%theta = x' * y;
%H = min(x,y);
%theta = H'*H;
end