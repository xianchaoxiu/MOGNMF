function [S_fcls] = fcls(A, X)
    % 输入:
    % A: LxP 的端元矩阵
    % X: LxN 的高光谱数据矩阵，其中 L 是波段数，N 是像素数
    % 输出:
    % S_fcls: PxN 的丰度矩阵

    [L, P] = size(A);
    N = size(X, 2);
    S_fcls = zeros(P, N); % 初始化丰度矩阵

    % 定义 lsqlin 的选项：非负约束和列和为 1 的约束
    options = optimoptions('lsqlin', 'Display', 'off');

    % 逐像素求解 FCLS 问题
    Aeq = ones(1, P);  % 线性等式约束 (sum(s) == 1)
    beq = 1;

    for i = 1:N
        % 求解非负最小二乘问题，带有线性等式约束 (sum(s) == 1)
        S_fcls(:, i) = lsqlin(A, X(:, i), [], [], Aeq, beq, zeros(P, 1), [], [], options);
    end
end
