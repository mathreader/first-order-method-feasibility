
% 1. Results of experiments with different dimensions of matrices

% Matrix for number of iterations used for [200, 400, 600, 800] variables
% with [1000, 2000, 3000, 4000] coordinates

A = load('polyak_avg_iter_200_800_1000_4000_4_4.mat');
B = load('restart_avg_iter_200_800_1000_4000_4_4.mat');

% Rows represent number of variables, column represent number of restarts

A.store_polyak ./ B.store_restart

% Matrix for number of iterations used for [100, 200, 400, 800] variables
% with [1000, 2000, 4000, 8000] coordinates

A = load('polyak_avg_iter_100_800_1000_8000_4_4_exp.mat');
B = load('restart_avg_iter_100_800_1000_8000_4_4_exp.mat');

A.store_polyak ./ B.store_restart


% 2. Results of experiments with different # of scaling in matrices

% We scale the coordinates by factor of 10

A = load('polyak_avg_iter_scale_800_0_2000_9.mat');
B = load('restart_avg_iter_scale_800_0_2000_9.mat');

A.store_polyak ./ B.store_restart


% 3. Results of experiments with different # of scaling in matrices

A = load('polyak_avg_sp_0_2_1_0_5.mat');
B = load('restart_avg_sp_0_2_1_0_5.mat');

A.store_polyak ./ B.store_restart

