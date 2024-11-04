% Figure_1D_AccuracyOfObjectiveFunction.m  
%                       
% This script sets up and computes error bounds associated with
% approximating the objective function and gradient for the 1D problem. 
%
% It computes the following:
%   (i)    - The relative error of the (genGK) approximation of the 
%            objective function to its exact formulation; this error is 
%            computed at the OPTIMAL hyperparameter values. Thus, we 
%            optimize the full problem first, and then make evaluations at 
%            the optimal hyperparameter values
%   (ii)   - The quadratic component of the relative error
%   (iii)  - The log-determinant component of the relative error
%   (iv)   - The relative error of the GSVD approximation (as compared
%            to the genGK approximation)
%   (v)    - The the generalized singular values associated to the problem
%
% This code is used to generate Figures 1 and 2 in the paper.
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

clear; close all; clc;
%%% Fix random seed %%%
rng('default')

%% Problem Initialization
%%% Generate forward operator (A), data (d), and true solution (s_true) %%%
n = 2^8;
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

%%% Choose prior %%%
prior_type = 'P1'; % noninformative
[~, gradP, logP] = Prior(prior_type);
inv.prior_type = prior_type;

%%% Choose kernel %%%
ker_name = 'Matern_3by2';
[kernel,gradkernel] = ker_fcn(ker_name);

Q = priorCov(xmin, xmax, nvec, scale, ker_name);
inv.Q = Q;

%% Optimizer (FULL)
optim_procedure = 'interior-point';
theta_0 = rand(3,1);
[t_opt_full, fval1, exitflag1, output1, xvals1] = optimizer(@objfun_full, theta_0, inv, optim_procedure);

%% Compute function evaluations

[F_exact, gradF_exact, val_exact] = objfun_full(t_opt_full, inv);

F1_exact = val_exact.F(1);
F2_exact = val_exact.F(2);

num_rec_terms = 80;

% genGK
F_approx_array = nan(num_rec_terms,1);
F1_approx_array = nan(num_rec_terms,1);
F2_approx_array = nan(num_rec_terms,1);

% GSVD
G_approx_array = nan(num_rec_terms,1);
G1_approx_array = nan(num_rec_terms,1);
G2_approx_array = nan(num_rec_terms,1);

parfor i = 1:num_rec_terms
    disp(['Iteration #',num2str(i)])
    
    [F_approx, gradF_approx, val_approx] = objfun_gengk(t_opt_full, inv, i);
    
    F_approx_array(i,1) = F_approx;
    F1_approx_array(i,1) = val_approx.F(1);
    F2_approx_array(i,1) = val_approx.F(2);
    
    [G_approx, gradG_approx, valG_approx] = objfun_gsvd(t_opt_full, inv, i);
    
    G_approx_array(i,1) = G_approx;
    G1_approx_array(i,1) = valG_approx.F(1);
    G2_approx_array(i,1) = valG_approx.F(2);
end

%%% Absolute and relative errors %%%
F_abs_err = abs(F_exact - F_approx_array);
F_rel_err = F_abs_err/abs(F_exact);

F1_abs_err = abs(F1_exact - F1_approx_array);
F1_rel_err = F1_abs_err/abs(F1_exact);

F2_abs_err = abs(F2_exact - F2_approx_array);
F2_rel_err = F2_abs_err/abs(F2_exact);

G_abs_err = abs(F_exact - G_approx_array);
G_rel_err = G_abs_err/abs(F_exact);

G1_abs_err = abs(F1_exact - G1_approx_array);
G1_rel_err = G1_abs_err/abs(F1_exact);

G2_abs_err = abs(F2_exact - G2_approx_array);
G2_rel_err = G2_abs_err/abs(F2_exact);

%% Compute bound estimates (Monte Carlo estimator)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mc = 10;
inv.genGK_iter = num_rec_terms;
[err_bound, val_err_bound] = mc_err_bnd(t_opt_full, inv, n_mc);

%% Plots
%%% Plot relative errors along with bound
fS = 16;
fig1 = figure(1);
semilogy(F_rel_err,'b-','LineWidth',2)
hold on
semilogy(err_bound/abs(F_exact),'r--','Linewidth',2)
xlabel('Bidiag. Param. $k$','Interpreter','Latex')
title(['RE'],'Interpreter','Latex')
legend('genGK','bound','Interpreter','Latex')
grid on
ax1 = gca;
ax1.FontSize = fS;

%%
fig2 = figure(2);
semilogy(F1_rel_err,'b-','LineWidth',2)
hold on
semilogy(val_err_bound.LD/abs(F_exact),'r--','Linewidth',2)
xlabel('Bidiag. Param. $k$','Interpreter','Latex')
title(['RE logdet'],'Interpreter','Latex')
legend('genGK','bound', 'Interpreter','Latex')
grid on
ax2 = gca;
ax2.FontSize = fS;
ax2.YTick = ax1.YTick;
ax2.YLim = ax1.YLim;

%%
fig3 = figure(3);
semilogy(F2_rel_err,'b-','LineWidth',2)
hold on
semilogy(val_err_bound.QU/abs(F_exact),'r--','Linewidth',2)
xlabel('Bidiag. Param. $k$','Interpreter','Latex')
title(['RE quadratic'],'Interpreter','Latex')
legend('genGK','bound','Interpreter','Latex')
grid on
ax3 = gca;
ax3.FontSize = fS;
ax3.YTick = [10^-10,1];

%% 
fig4 = figure(4);
semilogy(F_rel_err,'b-','LineWidth',2)
hold on
semilogy(G_rel_err,'k-.','LineWidth',2)
semilogy(22,F_rel_err(22),'b*','MarkerSize',12,'LineWidth',2)
xlabel('Bidiag. Param./Index $k$','Interpreter','Latex')
title(['RE'],'Interpreter','Latex')
legend('genGK','GSVD','Interpreter','Latex')
grid on
ax4 = gca;
ax4.FontSize = fS;
ax4.YTick = [10^-10,1];

%% 
[R, ~] = R_mat(t_opt_full, M);
[Qm, ~] = Q_mat(inv.Q, t_opt_full);
R_half = sqrtm(R);
S = svd(R_half\(full(A)*sqrtm(Qm*eye(N))));

fig5 = figure(5);
semilogy(S(1:num_rec_terms),'k-.','LineWidth',2)
xlabel('Index $k$','Interpreter','Latex')
title('Generalized Singular Values','Interpreter','Latex')
grid on
ax5 = gca;
ax5.FontSize = fS;
