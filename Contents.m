%% Contents.m
%   This folder contains MATLAB code to accompany the paper:
%
%     "Efficient iterative methods for hyperparameter estimation 
%     in large-scale linear inverse problems" 
%      - Hall-Hooper, Saibaba, Chung, and Miller, ACOM, 2024
%
%   The codes require the following packages:
%       REGU Regularization Tools by Per Christian Hansen
%             http://www.imm.dtu.dk/~pcha/Regutools/
%
%       genHyBR: generalized hybrid iterative methods
%             by Julianne Chung and Arvind K. Saibaba
%              https://github.com/juliannechung/genHyBR
% 
%   First run startup_EB.m for setting paths
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

%% Script files used to generate figures in Section 4.1 of the paper
% Inverse Heat Transfer Example

% Accuracy
Figure_1D_AccuracyOfObjectiveFunction
% creates Figures 1 and 2

% Reconstruction
Figure_1D_Reconstruction
% creates Figure 3

% Along trajectory
Figure_1D_BoundAlongTrajectory
% creates Figure 4

% Experiment 3: computation time
Figure_1D_ComputationalTime
% creates Figure 5