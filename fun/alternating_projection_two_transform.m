function [TX_1, TX_2, iter] = alternating_projection_two_transform(X, T_1, T_2, r_1, r_2, eps, max_iter, tol)
    % gen a rank r TX based on X as init
    T_1_pinv = pinv(T_1);
    TX_1 = mode3_transform(X, T_1);
    r_t_1 = tubal_rank(TX_1, tol);

    T_2_pinv = pinv(T_2);
    TX_2 = mode3_transform(X, T_2);
    r_t_2 = tubal_rank(TX_2, tol);

    for iter = 1 : max_iter
        X_t = X;
        
        % X -> T_pinv(svd(T(X)))
        TX_1 = mode3_transform(X, T_1);
        % project TX to rank-r manifold
        [U, S, V] = pagesvd(TX_1, 'econ');
        max_singular_value = max(S, [], 'all');
        for k = 1 : size(T_1, 1)
            if sum(diag(S(:, :, k)) > tol * max_singular_value) > r_1
                TX_1(:, :, k) = U(:, 1:r_1, k) * S(1:r_1, 1:r_1, k) * V(:, 1:r_1, k)';
            end
        end
        X = mode3_transform(TX_1, T_1_pinv);
        X = X / norm(X, 'fro');

        TX_2 = mode3_transform(X, T_2);
        % project TX to rank-r manifold
        [U, S, V] = pagesvd(TX_2, 'econ');
        max_singular_value = max(S, [], 'all');
        for k = 1 : size(T_2, 1)
            if sum(diag(S(:, :, k)) > tol * max_singular_value) > r_2
                TX_2(:, :, k) = U(:, 1:r_2, k) * S(1:r_2, 1:r_2, k) * V(:, 1:r_2, k)';
            end
        end
        X = mode3_transform(TX_2, T_2_pinv);
        X = X / norm(X, 'fro');

        if mod(iter, 100) == 0
            fprintf('iter: %d, ', iter)

            r_t_1 = tubal_rank(TX_1, tol);
            fprintf('rank of TX_1: %d, ', r_t_1)

            r_t_2 = tubal_rank(TX_2, tol);
            fprintf('rank of TX_2: %d, ', r_t_2)

            fprintf('adj max abs error: %.5e\n', max(abs(X(:) - X_t(:))))
        end
        
        if (max(abs(X(:) - X_t(:))) < eps) && (r_t_1 == r_1) && (r_t_2 == r_2)
            break;
        end
    end
