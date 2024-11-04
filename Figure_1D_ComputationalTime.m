% Figure_1D_ComputationalTime.m  
%
% This script computes wall closck time to compute the objective and
% gradient for different problem sizes.
%      
% This code is used to generate Figure 5 in the paper.      
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

clear; close all; clc;
%%% Fix random seed %%%
rng('default')

%% Fixed Parameters

kappa = 1;
Qflag = 1; % flag to specify dimension when generating Q
level = 0.02; % noise percentage/level

xmin = 0;                 % Coordinates of left of interval
xmax = 1;                 % Coordinates of right of interval
scale = 1;                % Parameters governing length scales.

%%% Choose prior %%%
prior_type = 'P2'; % Gamma
[~, gradP, logP] = Prior(prior_type);

%%% Choose kernel %%%
ker_name = 'Matern_3by2';
[kernel,gradkernel] = ker_fcn(ker_name);

num_samples = 20;
% I = 7:13; 
I = 7:11; 
F_time = nan(length(I),2,num_samples);
genGK_iter = 22; % bidiagonalization parameter'
inv.genGK_iter = genGK_iter;
theta = rand(3,1);

for j = 1:num_samples
    
    tic;
    %%% display sample number %%%
    disp(['sample #',num2str(j)])
    
    %% Timing
    for i = 1:length(I)
        %%% Generate forward operator (A), data (d), and true solution (s_true) %%%
        n = 2^I(i);
        
        %% Problem initialization
        %%% Generate forward operator (A), data (d), and true solution (s_true) %%%
        [A, d, s_true] = heat(n, kappa);
        
        M = size(A,1); inv.M = M;
        N = size(A,2); inv.N = N;
        
        %%% Add noise to data %%%
        eta = randn(size(d(:)));
        nN = norm(eta(:));
        eta = eta / nN;
        eta = level*norm(d(:))*eta;
        sigma = level*norm(d(:))/nN;
        
        d_noise = d + eta;
        
        %%% problem domain characteristic(s) %%%
        nvec = n;                 % Number of points in interval (from 'heat_example_setup')
        
        %%% Setting Inverse Problem Structure Parameter Values
        inv.A = A;
        inv.dn = d_noise;
        inv.s_true = s_true;
        inv.prior_type = prior_type;
        Q = priorCov(xmin, xmax, nvec, scale, ker_name);
        inv.Q = Q;
        
        %%% display problem size %%%
        disp(['size(A) = (',num2str(M),',',num2str(N),')'])
        
        %%% compute (objective,gradient) pair %%%
        tic;
        [F_approx, gradF_approx, ~] = objfun_gengk(theta, inv); % approximate problem
        tF_approx = toc;
        
        tic;
        [F_true, gradF_true, ~] = objfun_full(theta, inv); % full problem
        tF_true = toc;
        
        %%%
        F_time(i,:,j) = [tF_true,tF_approx];
        
    end
    t_samp = toc;
    disp(['time to complete sample = ',num2str(t_samp/60), ' minutes'])
    disp(' ')
end

F_ave = mean(F_time,3);
F_std = std(F_time,0,3);

disp('')

%% Slopes

% Full
T1 = [ones(size(log(F_ave(:,1)))),log((2.^I).')]; 
y1 = log(F_ave(:,1)); 
coeffs1 = T1\y1;
disp(['Slope of Exact Time Data (loglog space) = ',num2str(coeffs1(2))])
ff1 = @(x) exp(coeffs1(1))*(x.^coeffs1(2));

% genGK
T2 = [ones(size(log(F_ave(:,1)))),log((2.^I).')]; 
y2 = log(F_ave(:,2));
coeffs2 = T2\y2;
disp(['Slope of Approx. Time Data (loglog space) = ',num2str(coeffs2(2))])
ff2 = @(x) exp(coeffs2(1))*(x.^coeffs2(2));

%% Plot eval times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

fS = 16;
mS = 40;
lw = 1;

fig1 = figure(1);
p1 = loglog(2.^I,F_ave(:,1),'b.','LineWidth',lw,'MarkerSize',mS);
hold on
p2 = loglog(2.^I,F_ave(:,2),'r.','LineWidth',lw,'MarkerSize',mS);
xran = linspace(2^I(1), 2^I(end),100);
p3 = loglog(xran,ff1( xran ),'b-','LineWidth',lw,'MarkerSize',mS);
p4 = loglog(xran,ff2( xran ),'r-','LineWidth',lw,'MarkerSize',mS);

legend([p1,p2,p3,p4],'Full','genGK',...
    'Interpreter','Latex','Location','NorthWest')
ax = gca;

xticks((2.^I(1:end)))
xticklabels({'128','256','512','1024','2048','4096','8192'})

xlabel('$n$','Interpreter','Latex')
ylabel('Time (Seconds)','Interpreter','Latex')
title(['Time to Compute $\mathcal{F}$ and $\nabla \mathcal{F}$ (genGK $\# = $ ',num2str(genGK_iter), ')'],'Interpreter','Latex')

ax.FontSize = fS;
grid on
