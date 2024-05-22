function [X, iter] = min_transform_tnn(M, T, Omega, rho, alpha_0, alpha_max, eps, max_iter)
    X = M;
    T_pinv = pinv(T);
    alpha = alpha_0;
    T_X = mode3_transform(X, T);
    Z = 0;
    N_3 = size(T, 1);
    Y = T_X;
    for iter = 1 : max_iter
        X_t = X;
        Y_t = Y;
        T_X_t = T_X;
        % compute Y
        [U, S, V] = pagesvd(T_X + Z / alpha, 'econ');
        for k = 1 : N_3
            s = diag(S(:, :, k)) - (1/alpha);
            r = sum(s > 0);
            if r > 0
                Y(:, :, k) = U(:, 1:r, k) * diag(s(1:r)) * V(:, 1:r, k)';
            else
                Y(:, :, k) = 0;
            end
        end

        % compute X
        X = mode3_transform(Y - Z/alpha, T_pinv);
        X(Omega) = M(Omega);
        T_X = mode3_transform(X, T);

        % computing Z
        Z = Z + alpha * (T_X - Y);

        alpha = min(rho * alpha, alpha_max);

        if max([max(abs(X(:) - X_t(:))), max(abs(Y(:) - Y_t(:))), max(abs(T_X(:) - Y(:))), max(abs(T_X(:) - T_X_t(:)))]) < eps
            break;
        end
    end
end