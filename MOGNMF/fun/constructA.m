function A = constructA(X, k, d)
% k: 邻居数量
% d: 图的阶数
% A: Cell 数组，存储每个视角的多阶图
sigma=1;
for v = 1 : 2
    if v == 1
        % 空间视角：基于空间坐标计算欧几里得距离，并使用高斯核函数生成邻接矩阵
        [W_spatial, ~] = construct_distance_matrices(X, k, sigma);
        W = W_spatial; % 使用空间视角的邻接矩阵
    else
        % 光谱视角：基于光谱欧几里得距离生成邻接矩阵
        [~, W_spectral] = construct_distance_matrices(X, k, sigma);
        W = W_spectral; % 使用光谱视角的邻接矩阵
    end
    
    A{v, 1} = W; % 第一阶邻接矩阵
    
    % 生成高阶图
    for i = 2 : d
        A{v, i} = A{v, i-1} * A{v, 1}; % 递推计算高阶图
    end
end
end
% 根据空间和光谱距离计算邻接矩阵
function [W_spatial, W_spectral] = construct_distance_matrices(X, k, sigma)
[M, N, L] = size(X);
n = M * N;
[grid_x, grid_y] = meshgrid(1:M, 1:N);
coords = [grid_x(:), grid_y(:)]; 

% 计算空间距离并使用高斯核公式
W_spatial = zeros(n, n);
for i = 1:n
    for j = i+1:n % 只计算上三角部分
        dist_spatial = (coords(i,1) - coords(j,1))^2 + (coords(i,2) - coords(j,2))^2;
        W_spatial(i, j) = exp(-dist_spatial / (sigma^2));
        W_spatial(j, i) = W_spatial(i, j); % 直接填充对称元素
    end
end

% 展平高光谱数据
X = permute(X, [2, 1, 3]);
X_reshaped = reshape(X, M * N, L)';

W_spectral = zeros(n, n);
for i = 1:n
    for j = i+1:n % 只计算上三角部分
        dist_spectral = norm(X_reshaped(:,i) - X_reshaped(:,j))^2;
        W_spectral(i, j) = exp(-dist_spectral / (sigma^2));
        W_spectral(j, i) = W_spectral(i, j); % 直接填充对称元素
    end
end
% % 稀疏化邻接矩阵，保留每个像素的 k 个最近邻居
W_spatial = sparsify_matrix(W_spatial, k);
W_spectral = sparsify_matrix(W_spectral, k);
end

function W_sparse = sparsify_matrix(W, k)
n = size(W, 1);
W_sparse = zeros(n, n);

for i = 1:n
    [~, idx] = sort(W(i,:), 'descend'); 
    neighbors = idx(1:k);
    W_sparse(i, neighbors) = W(i, neighbors); 
end

W_sparse = (W_sparse + W_sparse') / 2;

end