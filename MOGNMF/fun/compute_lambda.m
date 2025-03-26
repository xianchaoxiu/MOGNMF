function lambda_e = compute_lambda(Y)
    % Y 是高光谱数据矩阵，尺寸为 [L, N]
    % L 是光谱带的数量（波段数量），N 是像素数量

    [L, N] = size(Y);  % 获取波段数量和像素数量

    % 初始化累加器
    sum_sparse = 0;

    % 遍历每一个波段，计算稀疏性度量
    for l = 1:L
        % 计算第 l 个波段的 L1 范数和 L2 范数
        norm_L1 = norm(Y(l, :), 1);   % L1 范数
        norm_L2 = norm(Y(l, :), 2);   % L2 范数

        % 计算每个波段的稀疏度贡献
        sparse_contrib = sqrt(N) - (norm_L1 / norm_L2);
        
        % 累加稀疏度贡献
        sum_sparse = sum_sparse + sparse_contrib;
    end

    % 计算 lambda_e
    lambda_e = (1 / sqrt(L)) * (sum_sparse / sqrt(N - 1));
end
