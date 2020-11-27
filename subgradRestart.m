function [x_final, k_polyak, restart_iter_store, restart_count, sol_type_flag] = ...
    subgradRestart(A, b, e, x0, eps_init, ...
    eps_shrink, restart_num, optim_f, max_iter_polyak, max_iter_restart)
%SUBGRADRESTART Restart method for subgradient method
%   Goal: Perform the subgradient method via the restart method
%   Inputs: A - m*n matrix representing coefficients of the inequalities
%                   representing the half spaces
%           b - m*1 vector representing constant value of the inequalities
%           e - n*m matrix, each column of e represents the initial point e_i.
%           x0 - n*1 vector representing starting point of the 
%                   subgradient method
%           eps_init - tolerance of error in the subgradient method
%           eps_shrink - shrinking of error 
%           restart_num - number of restart methods we're running
%           optim_f - optimal f, to be used with Polyak's method
%           max_iter_polyak - maximum number of iterations for Polyak's
%           method.
%           max_iter_restart - maximum number of iterations we would like to run


% 0. Initialization

% Fetch the dimension of matrices
m = size(A, 1);
n = size(A, 2);


% 1. Trivial checks

% Before we start the subgradient method, check if some e_i already 
% satisfy the conditions. If yes, print that value.
flag_complete = 0;

for i=1:1:n
    if (min(A * e(:, i) - b) >= 0)
        fprintf('One of initial points already satisfy the conditions.\n')
        x = e(:, i);
        sol_type_flag = 0;
        flag_complete = 1;
        break
    end
end


% 2. Main method

% Phase I: Run Polyak's method until we get to a point x_k satisfying
% gamma(x_k)<=3/2

x_polyak = x0;
k_polyak = 1;

% Precompute the constant b-diag(A*e) and c
temp_numer = b-diag(A*e);
c = b ./ temp_numer - 1;

temp_gamma = (A * x_polyak) ./ temp_numer - c;
[temp_gamma_max, ~] = max(temp_gamma);

if (flag_complete == 0)
    % Init the starting point and number of iterations for Polyak's method

    while (k_polyak <= max_iter_polyak)
        % Check convergence of current solution from Polyak's method
        if (min(A * x_polyak - b) >= 0)
            fprintf('Polyak subgradient method has converged to a solution.\n')
            sol_type_flag = 1;
            break
        end

        % Check whether we've reached close enough for Polyak's method
        % Compute gamma(x_polyak)
        temp_gamma = (A * x_polyak) ./ temp_numer - c;
        [temp_gamma_max, ~] = max(temp_gamma);

        if (temp_gamma_max <= 3 / 2)
            fprintf('Polyak subgradient method has converged to a close point.\n')
            break
        end

        % If we haven't converged to a close solution, apply one step of 
        % Polyak's method to current solution
        x_polyak = subgradStep(A, b, e, x_polyak, eps, optim_f, k_polyak, ...
                temp_numer, 2);

        % Update the iteration
        k_polyak = k_polyak + 1;
    end
end

% Phase II. Perform the restart method

if (k_polyak >= max_iter_polyak)
    fprintf('Polyak method failed to find a close point.')
elseif (flag_complete == 0)
    % Proceed to main loop of restart method.
    
    % (1) Initialize the subgradient methods
    % Initial points are all the polyak values
    restart_x_store = repmat(x_polyak, 1, restart_num);
    
    % Initial iteration numbers are set to 1.
    restart_iter_store = ones(restart_num, 1);
    
    % Initial iteration numbers are set to 1.
    restart_count = zeros(restart_num, 1);
    
    % Epsilon values of restart methods
    restart_eps = (eps_init * eps_shrink .^ (0:(restart_num-1)))';
    
    % Goal value for each restart method
    restart_goal = repmat(temp_gamma_max, restart_num, 1) - restart_eps;
    
    % Flag Indicates whether the subgradient method needs to compare its 
    % value with a value passed from the previous subgradient method
    % 0 implies not needed started, 1 implies started.
    restart_subgrad_flag = zeros(restart_num, 1);
    
    % Matrix to store the temporary passed values.
    restart_val_temp_store = zeros(n, restart_num);
    
    
    % (2) Main loop running the subgradient methods
    while (restart_iter_store(1) <= max_iter_restart && flag_complete == 0)
        for i=1:1:restart_num
            % Fetch current x value
            x_current = restart_x_store(:, i);

            % Check convergence of current subgradient method
            if (min(A * x_current - b) >= 0)
                fprintf(['Subgradient method with epsilon %d has converged ', ...
                    'to a solution.\n'], restart_eps(i))
                fprintf('This correspond to machine number %d.\n', i)
                sol_type_flag = 1+i;
                flag_complete = 1;
                break
            end

            % Fetch the current function value
            temp_gamma = (A * x_current) ./ temp_numer - c;
            [temp_gamma_max, ~] = max(temp_gamma);

            x_best = x_current;
            gamma_best = temp_gamma_max;

            % Check if there is another value passed from previous iteration
            % satisfying the assumption
            if (restart_subgrad_flag(i) == 1)
                x_restart = restart_val_temp_store(:, i);
                temp_gamma_restart = (A * x_restart) ./ temp_numer - c;
                [temp_gamma_restart_max, ~] = max(temp_gamma_restart);

                % If the restart value is better, replace the best by the
                % restart value
                if (temp_gamma_max > temp_gamma_restart_max)
                    x_best = x_restart;
                    gamma_best = temp_gamma_restart_max;
                end
            end


            % Check whether we should pass to the next subgradient method
            if (gamma_best < restart_goal(i) && i ~= restart_num)
                % Restart the current method
                restart_x_store(:, i) = x_best;
                restart_goal(i) = restart_goal(i) - restart_eps(i);

                % Pass value to the next vector
                restart_val_temp_store(:, i+1) = x_best;
                restart_subgrad_flag(i+1) = 1;

                % Clear inbox
                restart_subgrad_flag(i) = 0;

                restart_count(i) = restart_count(i) + 1;
                restart_iter_store(i) = 1;
            else 
                % Clear inbox
                restart_subgrad_flag(i) = 0;

                % Make one iteration
                x_current = subgradStep(A, b, e, x_current, restart_eps(i), ...
                    optim_f, restart_iter_store(i), temp_numer, 1);
                restart_x_store(:, i) = x_current;
                restart_iter_store(i) = restart_iter_store(i) + 1;
            end   
        end
    end
end


if (sol_type_flag == 1)
    x_final = x_polyak;
elseif (sol_type_flag >= 2)
    x_final = restart_x_store(:, sol_type_flag - 1);
else
    x_final = x;
end


end

