function [x_final, k, k_inner, sol_type_flag, eps_current] = subgradTwo(A,...
    b, e, x0, eps, eps_shrink, optim_f, max_iter)
%SUBGRADTWO Implement two subgradient methods
%   Goal: Perform the subgradient method via two layers of algorithms
%   Inputs: A - m*n matrix representing coefficients of the inequalities
%                   representing the half spaces
%           b - m*1 vector representing constant value of the inequalities
%           e - n*m matrix, each column of e represents the initial point e_i.
%           x0 - n*1 vector representing starting point of the 
%                   subgradient method
%           eps - tolerance of error in the subgradient method
%           eps_shrink - shrinking of error 
%           optim_f - optimal f, to be used with Polyak's method
%           max_iter - maximum number of iterations we would like to run

% Fetch the dimension of matrices
m = size(A, 1);
n = size(A, 2);

% Init the number of iterations
k = 1;

% Before we start the subgradient method, check if some e_i already 
% satisfy the conditions. If yes, print that value.
flag_complete = 0;

for i=1:1:n
    if (min(A * e(:, i) - b) >= 0)
        fprintf('One of initial points already satisfy the conditions.\n')
        x = e(:, i);
        % stored_steps = x;
        flag_complete = 1;
        break
    end
end


% Proceed to main loop if we don't satisfy the conditions.
if (flag_complete == 0)
    % Initialize useful parameters
    x_polyak = x0;
    x_inner = x0;
    k_inner = 0;
    eps_current = eps;
    
    % Indicates whether inner subgradient method is started.
    % 0 implies not started, 1 implies started.
    inner_subgrad_start_flag = 0;
    inner_subgrad_restart_flag = 0;
    
    % Indicate type of current solution
    % 0 implies no solution is found
    % 1 implies solution is found with Polyak's method
    % 2 implies solution is found with inner subgradient method
    sol_type_flag = 0;
    
    % Precompute the constant b-diag(A*e) and c
    temp_numer = b-diag(A*e);
    c = b ./ temp_numer - 1;
    
    % Precompute gamma(0)
    zero_gamma = - c;
    [zero_gamma_max, ~] = max(zero_gamma);
    
    while (k <= max_iter)
        % Part 1: Check of convergence and set start flag of inner 
        %   subgradient method
        
        % Check convergence of current solution from Polyak's method
        if (min(A * x_polyak - b) >= 0)
            fprintf('Polyak subgradient method has converged to a solution.\n')
            sol_type_flag = 1;
            break
        end
        
        % If current solution is not in the feasible region, 
        % decide whether the current solution is sufficiently close to 
        % start the subgradient method with fixed epsilon.
        
        % Compute the distance between current solution and 0
        temp_gamma = (A * x_polyak - b) ./ temp_numer + 1;
        [temp_gamma_max, ~] = max(temp_gamma);
        
        temp_dist = temp_gamma_max - zero_gamma_max;
        
        % If the current solution is smaller than eps_current, start the
        % inner subgradient method / restart the eps
        if (temp_dist < eps_current * eps_shrink)
            inner_subgrad_start_flag = 1;
            inner_subgrad_restart_flag = 1;
            eps_current = eps_current * eps_shrink;
            %eps_current
        end
        
        % Perform the inner subgradient method, if prompted
        if (inner_subgrad_start_flag == 1)
            % Restart the subgradient method, if needed
            if (inner_subgrad_restart_flag == 1)
                x_inner = x_polyak;
                k_inner = 0;
                inner_subgrad_restart_flag = 0;
            end
            
            % Check convergence of current solution from inner subgrad method
            if (min(A * x_inner - b) >= 0)
                fprintf('Inner subgradient method has converged to a solution.\n')
                sol_type_flag = 2;
                break
            end
            
            % Perform the inner subgradient method
            x_inner = subgradStep(A, b, e, x_inner, eps_current, optim_f, ...
                k_inner, temp_numer, 1);
            
            k_inner = k_inner + 1;
        end
        
        % Then apply one step of Polyak's method to current solution
        x_polyak = subgradStep(A, b, e, x_polyak, eps, optim_f, k, ...
            temp_numer, 2);
        k = k + 1;
    end  
end


if (sol_type_flag == 1)
    x_final = x_polyak;
elseif (sol_type_flag == 2)
    x_final = x_inner;
else
    x_final = x;
end

end

