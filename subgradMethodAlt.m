function [x, k, stored_steps] = subgradMethodAlt(A, b, e, x0, eps, optim_f, ...
    max_iter, step_size_flag)
%SUBGRADMETHODALT Ordinary subgradient method
%   Goal: Perform the subgradient method 
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


% Before we start the subgradient method, check if some e_i already 
% satisfy the conditions. If yes, print that value.
flag_complete = 0;

for i=1:1:n
    if (min(A * e(:, i) - b) >= 0)
        fprintf('One of initial points already satisfy the conditions.\n')
        x = e(:, i);
        k = 0;
        stored_steps = x;
        flag_complete = 1;
        break
    end
end


if (flag_complete == 0)
    % Init index and initial value of x
    k = 1;
    x = x0;

    % Initialize matrix to store the coordinates of each step
    stored_steps = zeros(n, max_iter+1);
    stored_steps(:, 1) = x0;
    
    % Precompute the constant b-diag(A*e)
    temp_numer = b-diag(A*e);
    
    fprintf('Starting subgradient method.\n')
    while (k <= max_iter)
        if (mod(k, 10000) == 0)
            fprintf('Iteration %d\n', k);
        end
        % Check convergence of current solution
        if (min(A * x - b) > 0)
            fprintf('Subgradient method has converged to a solution.\n')
            break
        end

        % If current solution doesn't fit, proceed to subgrad step:
        x = subgradStep(A, b, e, x, eps, optim_f, k, temp_numer, step_size_flag);

        % Store the vector x in the current iteration
        stored_steps(:, k+1) = x;

        % Update index
        k = k + 1;
    end
end

end

