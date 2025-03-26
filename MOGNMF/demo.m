clear all
clear all
clc
addpath(genpath('fun'));
addpath(genpath('data'));
% ********************************************************************** %
% load the groundtruth endmember and abudaces matrix
load 'F1_A_9'; % A L*M endmember matrix
load 'F1_S_9'; % S row*col*M abundance maps

[row, col, M] = size(S);
N = row*col;H = size(A,1); 
S = reshape(S, N, M)'; % colomnwise abundance matrix

% build 3D noiseless HSI
Y = reshape((A*S)',row,col,H); % For real data, please normalize the radiance to a range of 0 - 1.0.
% add white Gaussian noise+
SNR = 20;
noise_type = 'additive'; eta = 0;
[X, n, Cn] = addNoise (Y,noise_type,SNR, eta, 1);

% nonnegative obervation
 X = max(X,eps); % 2D HSI L*N
Y = reshape(X',row,col,H); % 3D HSI row*col*L 

% ## 2.Initialize A and S *********************************************** %
% use region-based vca to initiate A_init
P=9;
A_init =vca(X,P);
S_init = fcls(A_init,X);
% ## plot the initial results ******************************************* %
% show the initial accuracy(SAD and RMSE)
Sam = SAM(A, A_init); 
rmse(S, S_init, Sam(1,:), Sam(2,:));
%% *************** NMF *************** %%
 K=5;
sigma=1;
%find best solution
lambda = compute_lambda(X);
 mu=0.01;
 gamma=1.5;
[A_f, S_f] = getAS(X, A_init,S_init,Y,lambda,gamma,mu,P,N,A,S);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%  
% 打印最终的SAD和RMSE
Sam = SAM(A, A_f); 
rmse(S, S_f, Sam(1,:), Sam(2,:));