function r = tubal_rank(A, tol)
    S = pagesvd(A);
    r = 0;
    max_singular_value = max(S, [], 'all');
    for k = 1 : size(A, 3)
        r = max(sum(S(:, :, k) > tol * max_singular_value), r);
    end
end