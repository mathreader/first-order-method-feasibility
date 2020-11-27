
% Test 1. Testing the subgradient method with randomized matrices

% 1-a. Inputs
% n is the dimension of space we're working with 
n = 200;

% m is the number of inequalities we have
m = 800;


% 1-b. Generation of random half-spaces

% Currently, each entry of A come from a uniform distribution in [-1, 1]
% and each entry of b come from a uniform distribution in [0, 1]
A = rand(m, n) * 2 - 1;
b = -rand(m, 1);

% % Here is some code if we want the matrix to be shrinked by some factor
% w = ones(m, 1);
% w(1:100) = 100;
% A = repmat(w, 1, n) .* A;
% b(1:100) = b(1:100)/10;


% 1-c. Generation of initial points e_i
% Each column of e represent e_i in the notes.
e = zeros(n, m);

for i=1:1:m
    e(:, i) = (b(i) + norm(A(i, :))) * A(i, :)' / (norm(A(i, :))^2);
end


% 1-d. Initialization of x_0
% x0 is picked to be some e_i. Currently, we just set it to e_1
x0 = e(:, 1);

% % Precompute the constant b-diag(A*e) and c
% temp_numer = b-diag(A*e);
% c = b ./ temp_numer - 1;
% 
% temp_gamma = (A * x0) ./ temp_numer - c;
% [temp_gamma_max, ~] = max(temp_gamma);
% temp_gamma_max


% 1-e. Initialization of step size factor eps
% We start with epsilon <= 1-gamma(0)
temp = zeros(n, 1);

gamma_zero = zeros(m, 1);
for i=1:1:m
    gamma_zero(i) = - (A(i, :) *  e(:, i)) / (b(i) - A(i, :) * e(:, i));
end
%gamma_zero
max_gamma_zero = max(gamma_zero);
eps = 1 - max_gamma_zero


% For the subgradient method with two different methods:
eps_start = 1;
eps_shrink = 1/2;


% 1-f. Convergence configurations
% max_iter is the maximum number of iterations we would like to run
max_iter = 100000;

% step_size_flag decides the way we update step size.
% step_size_flag = 0 means we're using 1/k as step size.
% step_size_flag = 1 means we're using eps/|g|^2 as step size 
% step_size_flag = 2 means we're using (f(x_k)-f')/|g|^2 as step size 
step_size_flag = 2;


% Upper bounding the number of iterations needed.
% Upper bound D(\overline{x}) by D(0)

% dist = zeros(n, 1);
% for i=1:1:m
%     dist = norm(e(:, i));
% end
% max_dist = max(dist);
% max_num_iter = (2*max_dist / eps)^2;
% fprintf('Upper bound on iterations is %d \n', max_num_iter);


% fprintf('1. Running the 1/n step size.\n')
% 
% [~, k, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
%     max_iter, 0);
% 
% fprintf('The number of iterations used is %d\n\n', k)
% 
fprintf('2. Running the eps/|g_k|^2 step size.\n')

[~, k, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
    max_iter, 1);

fprintf('The number of iterations used is %d\n\n', k)

fprintf('3. Running the (f(x_k)-f*)/|g_k|^2 step size.\n')

[~, k, ~] = subgradMethodAlt(A, b, e, x0, eps, max_gamma_zero, ...
    max_iter, 2);

fprintf('The number of iterations used is %d\n\n', k)


% fprintf(['4. Running the combination of two subgradient methods. ', ...
%     'Assume that we know the feasible point 0.\n'])
% 
% [x, k1, k2, flag, eps_current] = subgradTwo(A, b, e, x0, eps_start, eps_shrink, ...
%     max_gamma_zero, max_iter);
% 
% if (flag == 1)
%     fprintf('Solution obtained using Polyak method, with iterations %d\n\n', k1)
% elseif (flag == 2)
%     fprintf('Solution obtained using inner subgradient method.\n')
%     fprintf('Number of Polyak method iterations used is %d.\n', k1-k2)
%     fprintf('Number of fixed subgradient method used is %d.\n', k2)
%     fprintf('Fixed subgradient method ran at epslion %d\n\n', eps_current)
% end
% 
% 
% fprintf(['5. Running the combination of two subgradient methods.', ...
%     'Assume that we do not know the feasible point 0.\n'])
% 
% [x, k1, k2, flag, eps_current] = subgradTwo(A, b, e, x0, eps_start, eps_shrink, ...
%     1, max_iter);
% 
% if (flag == 1)
%     fprintf('Solution obtained using Polyak method, with iterations %d\n', k1)
% elseif (flag == 2)
%     fprintf('Solution obtained using inner subgradient method.\n')
%     fprintf('Number of Polyak method iterations used is %d.\n', k1-k2)
%     fprintf('Number of fixed subgradient method used is %d.\n', k2)
%     fprintf('Fixed subgradient method ran at epslion %d\n\n', eps_current)
% end


% fprintf('6. Testing the subgradient method with restart:\n')
% 
% % Additional parameters
% eps_start = 1/2;
% eps_shrink = 1/2;
% restart_num = 20;
% max_iter_polyak = 1000;
% max_iter_restart = 100000;
% 
% [x_final,k_polyak,restart_iter_store,restart_count,flag] = subgradRestart(A, b, ...
%     e, x0, eps_start, eps_shrink, restart_num, max_gamma_zero, ...
%     max_iter_polyak, max_iter_restart);
% 
% if (flag == 1)
%     fprintf('Solution obtained using Polyak method, with iterations %d\n', k1);
% elseif (flag >= 2)
%     fprintf('Solution obtained using restart subgradient method.\n');
%     fprintf('Initialize with Polyak takes %d iterations.\n', k_polyak);
%     fprintf('Total iterations used is %d.\n', restart_iter_store(1));
%     
%     fprintf('Iterations from last restart for each subgradient method: \n')
%     restart_iter_store'
%     fprintf('Number of restarts for each subgradient method:\n')
%     restart_count'
% end


fprintf('7. Testing the subgradient method with Perceptron + fixed epsilon:\n')

[x, k, ~] = subgradPerceptron(A, b, e, 1/2, max_iter, 1);

min(A*x-b)

fprintf('The number of iterations used is %d\n\n', k)


fprintf('8. Testing the subgradient method with Perceptron + Polyak:\n')

[x, k, ~] = subgradPerceptron(A, b, e, 1/2, max_iter, 2);

min(A*x-b)

fprintf('The number of iterations used is %d\n\n', k)






