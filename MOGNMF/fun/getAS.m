function[A_F,S_F] = getAS(X, A_init,S_init,Y,lambda,gamma,mu,P,N,A,S)
W=A_init;
H=S_init;
% 更新 A_old 和 S_old
times=10;
eps = 1e-10;
E= zeros(size(X));
Wh=getW(Y);
Dh=diag(sum(Wh, 2));
Lh = Dh - Wh;
object_old = 0.5*norm((X-W*H), 'fro')^2+0.5*lambda*sum(sum(H.^0.5))+mu*0.5*trace(H*Lh*H')+gamma* norm(E, 'fro');;
% Step 5: 迭代优化
num_iterations = 1000;  % 设定迭代次数
for iteration = 1:num_iterations 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Step 5.1: 更新端元矩阵 W
     X1=X-E;
    t1=W*(H*H');
    t2=X1*H';
    t3=t2./max(t1,eps);
    W=W.*t3;
% Step 5.2: 更新丰度矩阵 H
    delta = 15;  % ASC 权重参数
    P=9;
    Wb = [W; delta*ones(1,P)];
    Vb = [X1; delta*ones(1,N)];
     HHt=H; 
    HHt(HHt==0) = eps;                  
    t1=( Wb'*Wb*H+0.5*lambda*(HHt.^(-0.5))+mu*H*Dh); 
    t2=( Wb'*Vb+mu*H*Wh ); 
    t3=(t2./max(t1,eps));  
    H1=H.*t3;
    %
    t1=( Wb'*Wb*H+mu*H*Dh);
    t2=(Wb'*Vb+mu*H*Wh);
    t3=(t1./max(t2,eps));  
    H2=H.*t3;   
    for i = 1:size(H,1)
    for j = 1:size(H,2)
        if H(i,j) >= 10^-4
            Htemp(i,j) = H1(i,j);  % 当 H(i,j) 大于等于 10^-4 时，取 H1(i,j)
        else
            Htemp(i,j) = H2(i,j);  % 当 H(i,j) 小于 10^-4 时，取 H2(i,j)
        end
    end
     end
  H = Htemp;
  E = updateE(X,W,H,gamma);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
  object_new = 0.5*norm((X-W*H), 'fro')^2+0.5*lambda*sum(sum(H.^0.5))+mu*0.5*trace(H*Lh*H')+gamma* norm(E, 'fro');
    ObjDiff =norm(object_old-object_new);
fprintf('Iter %d: ObjDiff: %.6f\n', iteration, ObjDiff);
    object_old = object_new; 
    if(ObjDiff < 10e-4)
        times=times+1;
        if (times ==10)           
           break;
        end
    end
end
A_F=W;
S_F=H;

end

function E = updateE(X,A,S, gamma)
R = X - A*S; % 计算残差矩阵
[m, n] = size(R);
E = zeros(m, n);
for i = 1:m
    row = R(i, :);
    norm_row = norm(row, 2);

    if norm_row > gamma
        factor = max(0, 1 - gamma / norm_row);
        E(i, :) = factor * row;
    else
        E(i, :) = zeros(1, n); 
    end
end
end
