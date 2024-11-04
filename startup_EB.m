%% startup file for Empirical Bayes to estimate hyperparameters
%
% These codes were used to generate the figures in the paper:
%     "Efficient iterative methods for hyperparameter estimation 
%     in large-scale linear inverse problems" 
%      - Hall-Hooper, Saibaba, Chung, and Miller, ACOM, 2024
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

% Regularization Tools
path(path,'./regu-master')

% genHyBR Codes
path(path,'./genHyBR-master')
path(path,'./genHyBR-master/genHyBR')
path(path,'./genHyBR-master/toeplitz')

% % Additional supporting files
path(path,'./EBcodes')

close all
clear all
clc