%  Figure_1D_BoundAlongTrajectory.m                    
% 
% This script investigates the bounds along the optimization trajectory.
%
% This code is used to generate Figure 4 in the paper.   
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
disp(['theta_0 = (',num2str(theta_0(1)),',',num2str(theta_0(2)),',',num2str(theta_0(3)),')'])
[t_opt, fval, exitflag, output, xvals] = optimizer(@objfun_gengk, theta_0, inv, optim_procedure);


%% 
TR = cell(1,1);
TR{1,1} = xvals;
DIFF_TRUE = cell(length(TR),3);
DIFF_APPROX = cell(length(TR),3);

%%

tInit = cputime;
for i = 1:length(TR)
    d1_true = nan(1,length(TR{i}));
    d2_true = nan(1,length(TR{i}));
    d_true = nan(1,length(TR{i}));
    d1_approx = d_true;
    d2_approx = d_true;
    d_approx = d_true;
    disp(['i = ', num2str(i)])

    tStart = cputime;
    for j = 1:length(TR{i})
        
        disp(['j = ', num2str(j)])
        tht = TR{i}(:,j);
        
        %%%
        [F_exact, ~, val_exact] = objfun_full(tht, inv);
        [F_approx, ~, val_approx] = objfun_gengk(tht, inv);
        
        d1_true(j) = abs(val_exact.F(1) - val_approx.F(1))/abs(val_exact.F(1));
        d2_true(j) = abs(val_exact.F(2) - val_approx.F(2))/abs(val_exact.F(2));
        d_true(j) = abs(F_exact - F_approx)/abs(F_exact);
        
        %%%
        n_mc = 10; % sample size for MC trace estimator
        [err_bound, val] = mc_err_bnd(tht, inv, n_mc);
        
        eb1 = val.LD;
        eb2 = val.QU;
        
        d1_approx(j) = eb1(end)/abs(val_exact.F(1));
        d2_approx(j) = eb2(end)/abs(val_exact.F(2));
        d_approx(j) = err_bound(end)/abs(F_exact);
    end
    tEnd = cputime - tStart;
    disp(['cputime for iteration = ',num2str(tEnd/60), ' minutes'])
    DIFF_TRUE{i,1} = d_true;
    DIFF_TRUE{i,2} = d1_true;
    DIFF_TRUE{i,3} = d2_true;
    
    DIFF_APPROX{i,1} = d_approx;
    DIFF_APPROX{i,2} = d1_approx;
    DIFF_APPROX{i,3} = d2_approx;
end
tFinal = cputime - tInit;
disp(['Total cputime elapsed = ', num2str(tFinal)]);

%% plots
close all

i = 1;

fs = 16;
lw = 2;
ms = 10;

fig1 = figure(1);
semilogy(DIFF_TRUE{i,1},'b-','LineWidth',lw,'MarkerSize',ms)
hold on
semilogy(DIFF_APPROX{i,1},'r--','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('Trajectory Iteration', 'Interpreter', 'Latex')
title('RE', 'Interpreter', 'Latex')
legend('genGK','bound','Location','SouthEast','Interpreter','Latex')
ax = gca;
ax.FontSize = fs;
yticks([10^(-10) 10^0])
%%%%%%%%%%

%%%%%%%%%%
fig2 = figure(2);
semilogy(DIFF_TRUE{i,2},'b-','LineWidth',lw,'MarkerSize',ms)
hold on
semilogy(DIFF_APPROX{i,2},'r--','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('Trajectory Iteration', 'Interpreter', 'Latex')
title('RE logdet', 'Interpreter', 'Latex')
legend('genGK','bound','Location','SouthEast','Interpreter','Latex')
ax = gca;
ax.FontSize = fs;
%%%%%%%%%%

%%%%%%%%%%
fig3 = figure(3);
semilogy(DIFF_TRUE{i,3},'b-','LineWidth',lw,'MarkerSize',ms)
hold on
semilogy(DIFF_APPROX{i,3},'r--','LineWidth',lw,'MarkerSize',ms)
grid on
xlabel('Trajectory Iteration', 'Interpreter', 'Latex')
title('RE quadratic', 'Interpreter', 'Latex')
legend('genGK','bound','Location','SouthEast','Interpreter','Latex')
ax = gca;
ax.FontSize = fs;
%%%%%%%%%%

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
