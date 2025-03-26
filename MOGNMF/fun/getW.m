function S=getW(Y)
[M, N, L] = size(Y); 
n = M * N; % 样本总数 n = M * N
num=n;
V=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=2;
 mu=10.^[-2];
% lambda=2^10;
beta=10.^[-1];
results=[];
K=5;
A = constructA(Y,K,d);
for i1=1:length(mu)
for i2=1:length(beta)
result=struct();
[S,obj,H] = solver(num,V,mu(i1),beta(i2),d,A);
 S(S<1e-5)=0;
end
end
