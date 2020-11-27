
% This file tests the case when we put some sparsity into the matrix

% Fix a certain dimension of matrix
n = 400;
m = 2000;

rep_num = 5;
prob_nonzero = linspace(0.2, 1, 5);

store_polyak = zeros(5, 1);
store_restart = zeros(5, 1);


for k=1:1:rep_num
    % First generate the random matrix 
    A_orig = rand(m, n) * 2 - 1;
    b= -rand(m, 1);
    
    for prob_ind=1:1:5
        A = zeros(m, n);
        rand_mat = rand(m, n);
        
        % Insert sparsity by setting values to 0
        for i=1:1:m
            for j=1:1:n
                if rand_mat(i, j) < prob_nonzero(prob_ind)
                    A(i, j) = A_orig(i, j);
                else
                    A(i, j) = 0;
                end
            end
        end
        
        % Make sure no row has all 0's 
        %TODO
        
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
        fprintf('Running Polyak of iteration %d.\n', prob_ind)

        [~, l, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
            max_iter, 2);

        store_polyak(prob_ind) = store_polyak(prob_ind) + l;

        % Main method for running restart method
        fprintf('Running restart method of iteration %d.\n', prob_ind)

        [~,k_polyak,restart_iter_store,~,~] = ...
            subgradRestart(A, b, e, x0, eps_start, eps_shrink, ...
            restart_num, max_gamma_zero, max_iter_polyak, max_iter_restart);

        store_restart(prob_ind) = store_restart(prob_ind) + k_polyak + ...
            restart_iter_store(1); 
    end
end

store_polyak = store_polyak / rep_num
store_restart = store_restart / rep_num
store_polyak ./ store_restart

