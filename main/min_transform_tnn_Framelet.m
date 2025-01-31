function [X, T_X, iter] = min_transform_tnn_Framelet(M, Omega, rho, D, R, level, alpha_0, alpha_max, eps, max_iter, debug, X_T)
    [n_1, n_2, n_3] = size(M);
    N_3 = n_3 * size(D, 1) * level;
    
    wLevel = -1;
    nfilter = 1;
    mu_level = ones(n_1, n_2, N_3);
    if wLevel <= 0
        for ki=1 : level
            for ii = 1 : size(D, 1)
                mu_level(:,:,n_3*(ki-1)*size(D, 1)+(ii-1)*n_3+1: n_3*(ki-1)*size(D, 1)+ii*n_3) = mu_level(:,:,n_3*(ki-1)*size(D, 1)+(ii-1)*n_3+1 : n_3 * (ki - 1) * size(D, 1)+ii*n_3) * nfilter*norm(D(ii,:));
            end
            nfilter = nfilter * norm(D(ii,:));
        end
    else
        for ki = 1 : level
            for ii = 1 : size(D, 1)-1
                mu_level(:,:,n_3*(ki-1)*size(D, 1)+(ii-1)*n_3+1: n_3*(ki-1)*size(D, 1)+ii*n_3)=mu_level(:,:,n_3*(ki-1)*size(D, 1)+(ii-1)*n_3+1: n_3*(ki-1)*size(D, 1)+ii*n_3)*nfilter;
            end
            nfilter = nfilter * wLevel;
        end
    end

    X = M;
    alpha = alpha_0;
    T_X = Fold(FraDecMultiLevel(Unfold(X, size(X), 3), D, level), [n_1, n_2, N_3], 3);
    Z = 0;
    Y = T_X;
    for iter = 1 : max_iter
        X_t = X;
        Y_t = Y;
        T_X_t = T_X;
        % compute Y
        [U, S, V] = pagesvd(T_X + Z / alpha, 'econ');
        for k = 1 : N_3
            s = diag(S(:, :, k)) - (mu_level(1, 1, k) / alpha);
            r = sum(s > 0);
            if r > 0
                Y(:, :, k) = U(:, 1:r, k) * diag(s(1:r)) * V(:, 1:r, k)';
            else
                Y(:, :, k) = 0;
            end
        end

        % compute X
        X = Fold(FraRecMultiLevel(Unfold(Y - Z/alpha, size(T_X), 3), R, level), size(X),3);
        X(Omega) = M(Omega);
        T_X = Fold(FraDecMultiLevel(Unfold(X, size(X), 3), D, level), [n_1, n_2, N_3], 3);

        % computing Z
        Z = Z + alpha * (T_X - Y);

        alpha = min(rho * alpha, alpha_max);

        if max([max(abs(X(:) - X_t(:))), max(abs(Y(:) - Y_t(:))), max(abs(T_X(:) - Y(:))), max(abs(T_X(:) - T_X_t(:)))]) < eps
            break;
        end

        if debug && (mod(iter, 10) == 0 || iter == 1)
            fprintf('iter: %d, ', iter)
            fprintf('alpha: %.2f, ', alpha)
            fprintf('adj abs diff: %.5f, ', max([max(abs(X(:) - X_t(:))), max(abs(Y(:) - Y_t(:))), max(abs(T_X(:) - Y(:))), max(abs(T_X(:) - T_X_t(:)))]))
            fprintf('mPSNR: %.2f, mSSIM: %.2f\n', [mPSNR(X, X_T), mSSIM(X, X_T)])
        end
    end
end
