
% This file tests the subgradient method with different matrix sizes

% 1-a. Inputs
% n is the dimension of space we're working with 
%n_range = linspace(200,800,4);

n_range = [100, 200, 400, 800];

% m is the number of inequalities we have
%m_range = linspace(1000,4000,4);

m_range = [1000, 2000, 4000, 8000];

% Repetition of experiments
rep_num = 5;

store_polyak = zeros(4, 4);
store_restart = zeros(4, 4);

for n_ind=1:1:4
    for m_ind=1:1:4
        n = n_range(n_ind);
        m = m_range(m_ind);
        
        fprintf('Running with %d variables and %d inequalities.\n', n, m)
        avg_polyak = 0;
        avg_restart = 0;
        for k=1:1:rep_num
            % 1. For each input dimension, generate
            % Currently, each entry of A come from a uniform distribution in [-1, 1]
            % and each entry of b come from a uniform distribution in [0, 1]
            A = rand(m, n) * 2 - 1;
            b = -rand(m, 1);
            
            % 2. Generation of initial points e_i
            % See notes 2.1(2). Each column of e represent e_i in the notes.
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
            fprintf('Running Polyak of iteration %d.\n', k)

            [~, l, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
                max_iter, 2);
            
            avg_polyak = avg_polyak + l;
            
            % Main method for running restart method
            fprintf('Running restart method of iteration %d.\n', k)
            
            [~,k_polyak,restart_iter_store,~,~] = ...
                subgradRestart(A, b, e, x0, eps_start, eps_shrink, ...
                restart_num, max_gamma_zero, max_iter_polyak, max_iter_restart);
            
            avg_restart = avg_restart + k_polyak + restart_iter_store(1);
        end
        store_polyak(n_ind, m_ind) = avg_polyak / rep_num;
        store_restart(n_ind, m_ind) = avg_restart / rep_num;
    end
end

store_polyak
store_restart
store_polyak ./ store_restart

