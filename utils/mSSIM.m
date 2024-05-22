function val = mSSIM(A, B)
    assert(all(size(A) == size(B)))
    val = [];
    for i = 1 : size(A, 3)
        val = [val, ssim(A(:, :, i), B(:, :, i))];
    end
    val = mean(val);
end