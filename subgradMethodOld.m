function [x, k, stored_steps] = subgradMethodOld(A, b, e, x0, eps, optim_f, ...
    max_iter, step_size_flag)
%SUBGRADMETHODOLD Ordinary subgradient method
%   Goal: Perform the subgradient method 
%   Inputs: A - m*n matrix representing coefficients of the inequalities
%                   representing the half spaces
%           b - m*1 vector representing constant value of the inequalities
%           e - n*m matrix, each column of e represents the initial point e_i.
%           x0 - n*1 vector representing starting point of the 
%                   subgradient method
%           eps - tolerance of errorin the subgradient method
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
    
    % Precompute the constant C
    temp_numer = b-diag(A*e);
    
    c = b ./ temp_numer - 1;
    
    fprintf('Starting subgradient method.\n')
    while (k <= max_iter)
        if (mod(k, 10000) == 0)
            fprintf('Iteration %d\n', k);
        end
        % Check convergence of current solution
        if (min(A * x - b) >= 0)
            fprintf('Subgradient method has converged to a solution.\n')
            break
        end

        % If current solution doesn't fit, proceed to subgrad step:

        % See Notes 3.1(4), find the maximum index of gamma
        gamma = (A * x) ./ temp_numer - c;
        [gamma_max, j] = max(gamma);
        %for i=1:1:m
        %    gamma(i) = (A(i, :) * (x - e(:, i))) / (b(i) - A(i, :) * e(:, i));
        %end

        % See Notes 3.1(5), find the subgradient. 
        % pi is replaced by phi to avoid clashing names with the constant pi
        phi = e(:, j) + (b(j) - A(j, :)*e(:, j)) * ...
            (x - e(:, j)) / (A(j, :) * (x - e(:, j)));

        % See Notes 3.1(6), this gives the update step
        g = A(j, :)' / (A(j, :) * (phi - e(:, j)));
        
        if (step_size_flag == 0)
            x = x - g / k;
        elseif (step_size_flag == 1) 
            x = x - eps * g / (norm(g)^2);
        else
            x = x - (gamma_max - optim_f) * g / (norm(g)^2);
        end

        % Store the vector x in the current iteration
        stored_steps(:, k+1) = x;

        % Update index
        k = k + 1;
    end
end

end

