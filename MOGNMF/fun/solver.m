function [S, obj, H] = solver(num, V, mu, beta, K, A)
    H = ones(V, K) ./ (V * K);
    SS = 0;
    for v = 1:V
        SS = SS + A{v, 1};
    end

    SS = (SS + SS') / 2;
    iter =20;
    tolerance = 1e-6;  % 设置目标函数值变化的阈值
    prev_obj = inf;    % 初始目标函数值设为无穷大

    for t = 1:iter
        % 更新 S
        B = zeros(num, num);
        for v = 1:V
            for k = 1:K
                B = B + H(v, k) * A{v, k};
            end
        end

        S = max(0, B / (2 * mu + 2));

        % 更新 H
        R = zeros(V * K, V * K);
        P = zeros(V, K);
        for v = 1:V
            R(K * (v - 1) + 1:K * v, K * (v - 1) + 1:K * v) = ones(K, K);
            for k = 1:K
                P(v, k) = sum(sum((S - A{v, k}).^2));
            end
        end
        R = 0 * R + diag(ones(V * K, 1));
        Q = reshape(P', [], 1);
        H_ba = quadprog(2 * beta * R, Q, [], [], ones(1, K * V), 1, zeros(K * V, 1), ones(K * V, 1), [], optimset('Display', 'off'));
        tempH = reshape(H_ba, [K, V]);
        H = tempH';
         H = H / sum(H(:)); 

        obj1 = reshape(P, 1, []) * reshape(H, [], 1);
        obj(t) = obj1 + beta * H_ba' * R * H_ba + mu * trace(S' * S);

        % 判断是否收敛
        if abs(obj(t) - prev_obj) < tolerance
            fprintf('Converged at iteration %d\n', t);
            break;
        end

        prev_obj = obj(t);  % 更新上一次的目标函数值
        fprintf('iter: %d, obj: %.6f\n', t, obj(t));
    end
end
