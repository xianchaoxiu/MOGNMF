function [A_vca, indices] = vca(X, P)
    % 输入:
    % X: LxN 的高光谱数据矩阵，其中 L 是波段数，N 是像素数
    % P: 端元的数量
    % 输出:
    % A_vca: LxP 的端元矩阵
    % indices: 选择的端元像素的索引

    [L, N] = size(X);  % L 是波段数，N 是像素数

    % 1. 数据中心化（去均值处理）
    rMean = mean(X, 2);  % 计算每个波段的均值
    X_zeroMean = X - repmat(rMean, 1, N);  % 中心化数据

    % 2. 使用PCA进行初步降维
    [U, ~, ~] = svds(X_zeroMean * X_zeroMean' / N, P);  % 使用奇异值分解保留前 P 个成分
    X_projected = U' * X_zeroMean;  % 将数据投影到低维空间

    % 3. 初始化端元矩阵和索引数组
    A_vca = zeros(L, P);  % 端元矩阵
    indices = zeros(1, P);  % 存储端元对应的像素索引

    % 4. 迭代找到 P 个端元
    for i = 1:P
        % 4.1 计算投影后的数据每列的范数
        projection = sum(X_projected.^2, 1);  % 每列的 L2 范数

        % 4.2 找到投影范数最大的像素索引
        [~, maxIdx] = max(projection);
        indices(i) = maxIdx;  % 存储该端元的索引

        % 4.3 将该像素加入端元矩阵
        A_vca(:, i) = X(:, maxIdx);

        % 4.4 更新投影矩阵，移除该端元的影响
        if i < P
            % 计算投影矩阵并从数据中减去端元的贡献
            endmember = X_projected(:, maxIdx);
            projection_update = endmember' * X_projected / norm(endmember)^2;
            X_projected = X_projected - endmember * projection_update;
        end
    end
end
