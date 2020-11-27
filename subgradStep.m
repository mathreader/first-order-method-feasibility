function [x_new] = subgradStep(A, b, e, x, eps, optim_f, k, temp_numer, ...
    step_size_flag)
%SUBGRADSTEP Perform one single subgradient step
%   Goal: Perform one single subgradient step given input
%   Inputs: A - m*n matrix representing coefficients of the inequalities
%                   representing the half spaces
%           b - m*1 vector representing constant value of the inequalities
%           e - n*m matrix, each column of e represents the initial point e_i.
%           x - n*1 vector representing the current point evaluated by the
%                   subgradient method
%           eps - tolerance of error in the subgradient method
%           optim_f - optimal f, to be used with Polyak's method
%           k - current iteration number, used in certain subgradient
%               methods
%           temp_numer - Constant b-diag(A*e), used in the computations
%           step_size_flag - 0 if we're using 1/k as step size.
%                            1 if we're using eps/|g|^2 as step size.
%                            2 if we're using (f(x_k)-f')/|g|^2 
%   Outputs: x_new - the value of x after the subgradient step.


% See Notes 3.1(4), find the maximum index of gamma
gamma = (A * x - b) ./ temp_numer + 1;
[gamma_max, j] = max(gamma);

% See Notes 3.1(5), find the subgradient. 
% pi is replaced by phi to avoid clashing names with the constant pi
phi = e(:, j) + (b(j) - A(j, :)*e(:, j)) * ...
    (x - e(:, j)) / (A(j, :) * (x - e(:, j)));

% See Notes 3.1(6), this gives the update step
g = A(j, :)' / (A(j, :) * (phi - e(:, j)));

if (step_size_flag == 0)
    x_new = x - g / k;
elseif (step_size_flag == 1) 
    x_new = x - eps * g / (norm(g)^2);
else
    x_new = x - (gamma_max - optim_f) * g / (norm(g)^2);
end


end

