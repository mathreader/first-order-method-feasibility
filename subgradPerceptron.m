function [x, k, stored_steps] = subgradPerceptron(A, b, e, eps, ...
    max_iter, step_size_flag)
%SUBGRADPERCEPTRON Subgradient version of the perceptron algorithm
%   Goal: Perform the subgradient method with updates similar to that of 
%           the perceptron.
%   Inputs: A - m*n matrix representing coefficients of the inequalities
%                   representing the half spaces
%           b - m*1 vector representing constant value of the inequalities
%           e - n*m matrix, each column of e represents the initial point e_i.
%           x0 - n*1 vector representing starting point of the 
%                   subgradient method
%           eps - tolerance of error in the subgradient method
%           optim_f - optimal f, to be used with Polyak's method
%           max_iter - maximum number of iterations we would like to run
%           step_size_flag - 0 if we're using 1/k as step size.
%                            1 if we're using eps/|g|^2 as step size.
%                            2 if we're using (f(x_k)-f')/|g|^2 
% 
%   Outputs: k - number of iteration used
%            x - final vector obtained

% Fetch the dimension of matrices
m = size(A, 1);
n = size(A, 2);

% Construct the last interior point 
insert_e = zeros(n+1, 1);
insert_e(n+1) = 1;

% Construct the new matrix
new_A = [A -b; insert_e'];
%new_A

% Construct the final interior point 
new_e = [e zeros(n, 1); ones(1, m+1)];
%new_e

% RHS of the inequality is simply 0
new_b = zeros(m+1, 1);
%new_b

% Initialization is simply 0
new_x0 = zeros(n+1, 1);
%new_x0

[x, k, stored_steps] = subgradMethodAlt(new_A, new_b, new_e, new_x0, eps, ...
    0, max_iter, step_size_flag);

%new_A*x-new_b

y = x(1:n) / x(n+1);
x = y;

end

