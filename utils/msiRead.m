function msi_img = msiRead(path)
    arr = dir(path);
    cnt_cell = 1;
    arr_cell = {};
    for i = 1 : numel(arr)
        if strfind(arr(i).name, '.png')
            arr_cell{cnt_cell} = arr(i).name;
            cnt_cell = cnt_cell + 1;
        end
    end

    arr_cell = sort(arr_cell);
    msi_img = zeros(size(imread([path, '/', arr_cell{1}]), 1), size(imread([path, '/', arr_cell{1}]), 2), numel(arr_cell));
    for j = 1 : numel(arr_cell)
        img = imread([path, '/', arr_cell{j}]);
        msi_img(:, :, j) = img;
    end
end