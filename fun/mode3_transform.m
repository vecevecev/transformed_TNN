function A_transformed = mode3_transform(A, T)
    assert(size(A, 3) == size(T, 2))
    A_transformed = reshape((T * reshape(A, size(A, 1) * size(A, 2), size(A, 3)).').', size(A, 1), size(A, 2), size(T, 1));
end