% Figure_1D_Reconstruction.m  
%                             
% This script computes reconstructions for the 1D example.
%
% This code is used to generate Figure 3 in the paper.   
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

clear; close all; clc;
%%% Fix random seed %%%
rng('default')

%% Problem Initialization
%%% Generate forward operator (A), data (d), and true solution (s_true) %%%
n = 2^8; inv.n = n;
kappa = 1;
[A, d, s_true] = heat(n, kappa); inv.A = A; inv.s_true = s_true;

M = size(A,1); inv.M = M;
N = size(A,2); inv.N = N;

%%% Add noise to data %%%
level = 0.02; % noise percentage/level
[eta, sigma] = WhiteNoise(d, level);

d_noise = d + eta; inv.dn = d_noise;

%%% Choose problem domain %%%
xmin = 0;                 % Coordinates of left of interval
xmax = 1;                 % Coordinates of right of interval
nvec = n;                 % Number of points in interval (from 'heat_example_setup')
scale = 1;                % Parameters governing length scales.
nx = n; ny = 1;

%%% Choose prior %%%
prior_type = 'P1'; % noninformative
[~, gradP, logP] = Prior(prior_type);
inv.prior_type = prior_type;

%%% Choose kernel %%%
ker_name = 'Matern_3by2';
[kernel,gradkernel] = ker_fcn(ker_name);

Q = priorCov(xmin, xmax, nvec, scale, ker_name);
inv.Q = Q;

%% Optimizer (genGK) 
optim_procedure = 'interior-point';
theta_0 = rand(3,1);
inv.genGK_iter = 22;
[t_opt, fval, exitflag, output, xvals] = optimizer(@objfun_gengk, theta_0, inv, optim_procedure);

%% Perform Recovery
theta_opt = t_opt;

%%% Get optimal matrices/vectors %%%
[R, dR] = R_mat(theta_opt, M);
R = full(R);
[mu, dmu] = mu_vec(theta_opt, N);
[Qm, dQ] = Q_mat(inv.Q, theta_opt);

%%% Specify right-hand-side vector for pcg %%%
res = d_noise - A*mu;
b = (((res.'/R)*A)*Qm).';

%%% Recover signal %%%
tol = 1e-8;
maxit = 200;
[xx,fl0,rr0,it0,rv0] = pcgRecovery(A,Qm,R,mu,b,tol,maxit);

% Calculate Error Between Images
RE = norm(xx - s_true,2)/norm(s_true,2);
disp(['rel. err = ',num2str(RE)])

%% Reconstructions
close all
fS = 16;

%%%
figure(1)
titleLabel_1 = ['Data (',num2str(level*100),'$\%$ Noise)'];
pl1 = plot(linspace(0,1,n),d,'b-','LineWidth',2);
hold on
pl3 = plot(linspace(0,1,n),d_noise,'k.','MarkerSize',10);
title(titleLabel_1,'interpreter','latex','fontsize',fS)
legend([pl1, pl3],'Data','Data + Noise',...
    'interpreter','latex','fontsize',fS)
grid on
set(gca,'fontsize',fS)

%%%
figure(2)
titleLabel_2 = ['RE $=$ ',num2str(RE*100),'$\%$',' ($k = $',num2str(inv.genGK_iter),')'];
p21 = plot(linspace(0,1,n),s_true,'k','LineWidth',2);
hold on
p23 = plot(linspace(0,1,n),xx,'r--','LineWidth',2);
title(titleLabel_2,'interpreter','latex','fontsize',fS)
legend([p21, p23],'True','Reconstructed',...
    'interpreter','latex','fontsize',fS)
grid on
set(gca,'fontsize',fS)
