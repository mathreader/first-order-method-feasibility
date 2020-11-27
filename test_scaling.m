
% This file tests the subgradient method with same matrix size, while 
% different number of variables are scaled.

% Fix a certain dimension of matrix
n = 400;
m = 2000;

scale_num = linspace(0, 2000, 9);

% Repetition of experiments
rep_num = 5;

store_polyak = zeros(9, 1);
store_restart = zeros(9, 1);

for k=1:1:rep_num
    
    % First generate the random matrix 
    A = rand(m, n) * 2 - 1;
    b = -rand(m, 1);
    
    fprintf('Start itertaion %d.\n', k)
    
    for scale_ind=1:1:9
        % 1. Scale the matrix with certain rows
        b(1:scale_num(scale_ind)) = b(1:scale_num(scale_ind))/10;
        
        % 2. Generation of initial points e_i
        e = zeros(n, m);
        for i=1:1:m
            e(:, i) = (b(i) + norm(A(i, :))) * A(i, :)' / (norm(A(i, :))^2);
        end
        
        % 3. Initialization of x_0
        % x0 is picked to be some e_i. Currently, we just set it to e_1
        x0 = e(:, 1);
        
        % 4. Initialization of step size factor eps
        % We start with epsilon <= 1-gamma(0)
        temp = zeros(n, 1);
        
        gamma_zero = zeros(m, 1);
        for i=1:1:m
            gamma_zero(i) = - (A(i, :) *  e(:, i)) / (b(i) - A(i, :) * e(:, i));
        end
        max_gamma_zero = max(gamma_zero);
        eps = 1 - max_gamma_zero;
        
        % 5. Convergence configurations
        % max_iter is the maximum number of iterations we would like to run
        max_iter = 1000000;

        % 6. Additional parameters for the restart method
        eps_start = 1/2;
        eps_shrink = 1/2;
        restart_num = 20;
        max_iter_polyak = 1000;
        max_iter_restart = 100000;
        
        % Main method for running Polyak's method
        fprintf('Running Polyak of iteration %d.\n', scale_ind)

        [~, l, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
            max_iter, 2);

        store_polyak(scale_ind) = store_polyak(scale_ind) + l;

        % Main method for running restart method
        fprintf('Running restart method of iteration %d.\n', scale_ind)

        [~,k_polyak,restart_iter_store,~,~] = ...
            subgradRestart(A, b, e, x0, eps_start, eps_shrink, ...
            restart_num, max_gamma_zero, max_iter_polyak, max_iter_restart);

        store_restart(scale_ind) = store_restart(scale_ind) + k_polyak + ...
            restart_iter_store(1);
        
        b(1:scale_num(scale_ind)) = b(1:scale_num(scale_ind))*10;
    end
    
end

store_polyak = store_polyak / rep_num
store_restart = store_restart / rep_num
store_polyak ./ store_restart


