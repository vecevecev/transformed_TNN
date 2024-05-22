% main script of transformed TNN
% To use other transforms, replace the function ''min_transform_tnn_Framelet'' by ''min_transform_tnn'''


clear; close all; clc

rng(0)

addpath(genpath('../'))

load foreman.mat

num_frame = 100;
X = X(:, :, 1:num_frame);

[n_1, n_2, n_3] = size(X);


tol = 1e-3;

% sample
p = 0.5;
omega = find(rand(n_1 * n_2 * n_3, 1) < p);
M = zeros(n_1, n_2, n_3);
M(omega) = X(omega);

fprintf('Framelet:\n')
% Generate Framelet Filters (B_spline)
frame = 1;
level = 4;
[D, R] = GenerateFrameletFilter(frame);
N_3 = n_3 * size(D, 1) * level;
% min transform tnn
rho = 1.2;
alpha_0 = 1e-2;
alpha_max = 1e6;
eps = 1e-2;
max_iter = 3e3;
debug = 1;
[X_hat, TX_hat, num_iter] = min_transform_tnn_Framelet(M, omega, rho, D, R, level, alpha_0, alpha_max, eps, max_iter, debug, X);
TX = Fold(FraDecMultiLevel(Unfold(X, size(X), 3), D, level), [n_1, n_2, N_3], 3);

mPSNR = mPSNR(X_hat, X);
mSSIM = mSSIM(X_hat, X);

fprintf('num iter: %d\n', num_iter)
fprintf('tubal rank of TX: %d, estimated TX: %d\n', [tubal_rank(TX, tol), tubal_rank(TX_hat, tol)]);
fprintf('tubal rank of X: %d and estimated X: %d\n', [tubal_rank(X, tol), tubal_rank(X_hat, tol)]);
fprintf('mPSNR: %.2f, mSSIM: %.2f\n', [mPSNR, mSSIM]);


for i_frame = 1: num_frame
    subplot(1,2,1)
    imshow(X_hat(:, :, i_frame))
    title('Framelet')
    subplot(1,2,2)
    imshow(X(:, :, i_frame))
    title('Original')
    pause(0.03)
end