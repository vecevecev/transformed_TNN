clear; clc; close all;

addpath(genpath('../'))

n = 50;
n_1 = n;
n_2 = n;
n_3 = 20;

tol = 1e-3;

% 生成TX
eps_gen_TX = 1e-10;
max_iter_gen_TX = 1e5;

num_monte = 100;

v_1 = zeros(num_monte, 1);
v_2 = zeros(num_monte, 1);
v_3 = zeros(num_monte, 1);
v_4 = zeros(num_monte, 1);
v_5 = zeros(num_monte, 1);
v_6 = zeros(num_monte, 1);
v_7 = zeros(num_monte, 1);
v_8 = zeros(num_monte, 1);
v_9 = zeros(num_monte, 1);
v_10 = zeros(num_monte, 1);
v_11 = zeros(num_monte, 1);
v_12 = zeros(num_monte, 1);


arr_v_1_1 = zeros(num_monte, 1);
arr_v_1_2 = zeros(num_monte, 1);
arr_v_1_3 = zeros(num_monte, 1);
arr_v_1_4 = zeros(num_monte, 1);
arr_v_1_5 = zeros(num_monte, 1);
arr_v_1_6 = zeros(num_monte, 1);

arr_v_2_1 = zeros(num_monte, 1);
arr_v_2_2 = zeros(num_monte, 1);
arr_v_2_3 = zeros(num_monte, 1);
arr_v_2_4 = zeros(num_monte, 1);
arr_v_2_5 = zeros(num_monte, 1);
arr_v_2_6 = zeros(num_monte, 1);

arr_v_3_1 = zeros(num_monte, 1);
arr_v_3_2 = zeros(num_monte, 1);
arr_v_3_3 = zeros(num_monte, 1);
arr_v_3_4 = zeros(num_monte, 1);
arr_v_3_5 = zeros(num_monte, 1);
arr_v_3_6 = zeros(num_monte, 1);

for i_monte = 1 : num_monte
    disp(i_monte)
    T_1 = RandOrthMat_C(n_3);
    T_1 = T_1 / max(svd(T_1));
    T_2 = RandOrthMat_C(n_3);
    T_2 = T_2 / max(svd(T_2));
    T_1_pinv = pinv(T_1);
    T_2_pinv = pinv(T_2);
    
    r = 20;
        
    while true
        X_0 = 1/sqrt(2) * (randn(n_1, n_2, n_3) + 1j * randn(n_1, n_2, n_3));
        [TX_1, TX_2, num_iter_gen_TX] = alternating_projection_two_transform(X_0, T_1, T_2, r, r, eps_gen_TX, max_iter_gen_TX, tol);
        X_1 = mode3_transform(TX_1, T_1_pinv);
        X_2 = mode3_transform(TX_2, T_2_pinv);
        if (tubal_rank(TX_1, tol) == r) && (tubal_rank(TX_2, tol) == r) && (max(abs(X_1(:) - X_2(:))) < eps_gen_TX)
            X = X_1;
            break;
        end
    end
    
    T_3 = [T_1; T_2];
    T_3 = T_3 / max(svd(T_3));
    T_3_pinv = pinv(T_3);
    
    [U_1, S_1, V_1] = pagesvd(mode3_transform(X, T_1), 'econ');
    U_1 = U_1(:, 1:r, :);
    V_1 = V_1(:, 1:r, :);
    
    [U_2, S_2, V_2] = pagesvd(mode3_transform(X, T_2), 'econ');
    U_2 = U_2(:, 1:r, :);
    V_2 = V_2(:, 1:r, :);
    
    [U_3, S_3, V_3] = pagesvd(mode3_transform(X, T_3), 'econ');
    U_3 = U_3(:, 1:r, :);
    V_3 = V_3(:, 1:r, :);
    
    v_1(i_monte) = max(abs(T_1_pinv'), [], 'all') * max(abs(mode3_transform(pagemtimes(U_1, pagetranspose(conj(V_1))), T_1')), [], 'all');
    v_2(i_monte) = max(abs(T_2_pinv'), [], 'all') * max(abs(mode3_transform(pagemtimes(U_2, pagetranspose(conj(V_2))), T_2')), [], 'all');
    v_3(i_monte) = max(abs(T_3_pinv'), [], 'all') * max(abs(mode3_transform(pagemtimes(U_3, pagetranspose(conj(V_3))), T_3')), [], 'all');
    
    v_4(i_monte) = max(abs(T_1_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_1, pagetranspose(conj(V_1))), T_1')).^2, [1, 3])));
    v_5(i_monte) = max(abs(T_2_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_2, pagetranspose(conj(V_2))), T_2')).^2, [1, 3])));
    v_6(i_monte) = max(abs(T_3_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_3, pagetranspose(conj(V_3))), T_3')).^2, [1, 3])));
    
    v_7(i_monte) = max(abs(T_1_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_1, pagetranspose(conj(V_1))), T_1')).^2, [2, 3])));
    v_8(i_monte) = max(abs(T_2_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_2, pagetranspose(conj(V_2))), T_2')).^2, [2, 3])));
    v_9(i_monte) = max(abs(T_3_pinv'), [], 'all') * max(sqrt(sum(abs(mode3_transform(pagemtimes(U_3, pagetranspose(conj(V_3))), T_3')).^2, [2, 3])));
    
    %%
    v_1_1 = 0;
    v_2_1 = 0;
    v_3_1 = 0;
    v_1_2 = 0;
    v_2_2 = 0;
    v_3_2 = 0;
    
    v_1_3 = 0;
    v_2_3 = 0;
    v_3_3 = 0;
    v_1_4 = 0;
    v_2_4 = 0;
    v_3_4 = 0;
    
    
    v_1_5 = 0;
    v_2_5 = 0;
    v_3_5 = 0;
    v_1_6 = 0;
    v_2_6 = 0;
    v_3_6 = 0;
    for i = 1 : n_1
        for j = 1: n_2
            for k = 1 : n_3
                e_ijk = zeros(n_1, n_2, n_3);
                e_ijk(i, j, k) = 1;
                %% P_S
                P_S_T_e_ijk_1 = P_S(mode3_transform(e_ijk, T_1), U_1, V_1);
                P_S_T_tilde_e_ijk_1 = P_S(mode3_transform(e_ijk, T_1_pinv'), U_1, V_1);
                v_1_1 = max([v_1_1, max([norm(P_S_T_e_ijk_1, 'fro')^2, norm(P_S_T_tilde_e_ijk_1, 'fro')^2])]);
                v_1_2 = max([v_1_2, max([norm(mode3_transform(P_S_T_tilde_e_ijk_1, T_1'), 'fro')^2, norm(mode3_transform(P_S_T_e_ijk_1, T_1_pinv), 'fro')^2])]);
    
                P_S_T_e_ijk_2 = P_S(mode3_transform(e_ijk, T_2), U_2, V_2);
                P_S_T_tilde_e_ijk_2 = P_S(mode3_transform(e_ijk, T_2_pinv'), U_2, V_2);
                v_2_1 = max([v_2_1, max([norm(P_S_T_e_ijk_2, 'fro')^2, norm(P_S_T_tilde_e_ijk_2, 'fro')^2])]);
                v_2_2 = max([v_2_2, max([norm(mode3_transform(P_S_T_tilde_e_ijk_2, T_2'), 'fro')^2, norm(mode3_transform(P_S_T_e_ijk_2, T_2_pinv), 'fro')^2])]);
    
                P_S_T_e_ijk_3 = P_S(mode3_transform(e_ijk, T_3), U_3, V_3);
                P_S_T_tilde_e_ijk_3 = P_S(mode3_transform(e_ijk, T_3_pinv'), U_3, V_3);
                v_3_1 = max([v_3_1, max([norm(P_S_T_e_ijk_3, 'fro')^2, norm(P_S_T_tilde_e_ijk_3, 'fro')^2])]);
                v_3_2 = max([v_3_2, max([norm(mode3_transform(P_S_T_tilde_e_ijk_3, T_3'), 'fro')^2, norm(mode3_transform(P_S_T_e_ijk_3, T_3_pinv), 'fro')^2])]);
    
                %% P_U
                P_U_T_e_ijk_1 = P_U(mode3_transform(e_ijk, T_1), U_1, V_1);
                P_U_T_tilde_e_ijk_1 = P_U(mode3_transform(e_ijk, T_1_pinv'), U_1, V_1);
                v_1_3 = max([v_1_3, max([norm(P_U_T_e_ijk_1, 'fro')^2, norm(P_U_T_tilde_e_ijk_1, 'fro')^2])]);
                v_1_4 = max([v_1_4, max([norm(mode3_transform(P_U_T_tilde_e_ijk_1, T_1'), 'fro')^2, norm(mode3_transform(P_U_T_e_ijk_1, T_1_pinv), 'fro')^2])]);
    
                P_U_T_e_ijk_2 = P_U(mode3_transform(e_ijk, T_2), U_2, V_2);
                P_U_T_tilde_e_ijk_2 = P_U(mode3_transform(e_ijk, T_2_pinv'), U_2, V_2);
                v_2_3 = max([v_2_3, max([norm(P_U_T_e_ijk_2, 'fro')^2, norm(P_U_T_tilde_e_ijk_2, 'fro')^2])]);
                v_2_4 = max([v_2_4, max([norm(mode3_transform(P_U_T_tilde_e_ijk_2, T_2'), 'fro')^2, norm(mode3_transform(P_U_T_e_ijk_2, T_2_pinv), 'fro')^2])]);
    
                P_U_T_e_ijk_3 = P_U(mode3_transform(e_ijk, T_3), U_3, V_3);
                P_U_T_tilde_e_ijk_3 = P_U(mode3_transform(e_ijk, T_3_pinv'), U_3, V_3);
                v_3_3 = max([v_3_3, max([norm(P_U_T_e_ijk_3, 'fro')^2, norm(P_U_T_tilde_e_ijk_3, 'fro')^2])]);
                v_3_4 = max([v_3_4, max([norm(mode3_transform(P_U_T_tilde_e_ijk_3, T_3'), 'fro')^2, norm(mode3_transform(P_U_T_e_ijk_3, T_3_pinv), 'fro')^2])]);
    
                %% P_V
                P_V_T_e_ijk_1 = P_V(mode3_transform(e_ijk, T_1), U_1, V_1);
                P_V_T_tilde_e_ijk_1 = P_V(mode3_transform(e_ijk, T_1_pinv'), U_1, V_1);
                v_1_5 = max([v_1_5, max([norm(P_V_T_e_ijk_1, 'fro')^2, norm(P_V_T_tilde_e_ijk_1, 'fro')^2])]);
                v_1_6 = max([v_1_6, max([norm(mode3_transform(P_V_T_tilde_e_ijk_1, T_1'), 'fro')^2, norm(mode3_transform(P_V_T_e_ijk_1, T_1_pinv), 'fro')^2])]);
    
                P_V_T_e_ijk_2 = P_V(mode3_transform(e_ijk, T_2), U_2, V_2);
                P_V_T_tilde_e_ijk_2 = P_V(mode3_transform(e_ijk, T_2_pinv'), U_2, V_2);
                v_2_5 = max([v_2_5, max([norm(P_V_T_e_ijk_2, 'fro')^2, norm(P_V_T_tilde_e_ijk_2, 'fro')^2])]);
                v_2_6 = max([v_2_6, max([norm(mode3_transform(P_V_T_tilde_e_ijk_2, T_2'), 'fro')^2, norm(mode3_transform(P_V_T_e_ijk_2, T_2_pinv), 'fro')^2])]);
    
                P_V_T_e_ijk_3 = P_V(mode3_transform(e_ijk, T_3), U_3, V_3);
                P_V_T_tilde_e_ijk_3 = P_V(mode3_transform(e_ijk, T_3_pinv'), U_3, V_3);
                v_3_5 = max([v_3_5, max([norm(P_V_T_e_ijk_3, 'fro')^2, norm(P_V_T_tilde_e_ijk_3, 'fro')^2])]);
                v_3_6 = max([v_3_6, max([norm(mode3_transform(P_V_T_tilde_e_ijk_3, T_3'), 'fro')^2, norm(mode3_transform(P_V_T_e_ijk_3, T_3_pinv), 'fro')^2])]);
            end
        end
    end

    
    arr_v_1_1(i_monte) = v_1_1;
    arr_v_1_2(i_monte) = v_1_2;
    arr_v_1_3(i_monte) = v_1_3;
    arr_v_1_4(i_monte) = v_1_4;
    arr_v_1_5(i_monte) = v_1_5;
    arr_v_1_6(i_monte) = v_1_6;
    
    arr_v_2_1(i_monte) = v_2_1;
    arr_v_2_2(i_monte) = v_2_2;
    arr_v_2_3(i_monte) = v_2_3;
    arr_v_2_4(i_monte) = v_2_4;
    arr_v_2_5(i_monte) = v_2_5;
    arr_v_2_6(i_monte) = v_2_6;

    arr_v_3_1(i_monte) = v_3_1;
    arr_v_3_2(i_monte) = v_3_2;
    arr_v_3_3(i_monte) = v_3_3;
    arr_v_3_4(i_monte) = v_3_4;
    arr_v_3_5(i_monte) = v_3_5;
    arr_v_3_6(i_monte) = v_3_6;
end

function P_S_X = P_S(X, U, V)
    P_S_X = pagemtimes(pagemtimes(U, pagetranspose(conj(U))), X) + pagemtimes(X, pagemtimes(V, pagetranspose(conj(V)))) - pagemtimes(pagemtimes(U, pagetranspose(conj(U))), pagemtimes(X, pagemtimes(V, pagetranspose(conj(V)))));
end

function P_U_X = P_U(X, U, V)
    P_U_X = pagemtimes(pagemtimes(U, pagetranspose(conj(U))), X);
end

function P_V_X = P_V(X, U, V)
    P_V_X = pagemtimes(X, pagemtimes(V, pagetranspose(conj(V))));
end