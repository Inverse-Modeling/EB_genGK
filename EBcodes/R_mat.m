function [R, dR, invR] = R_mat(theta, M)
%
% [R, dR, invR] = R_mat(theta, M)
%
% This function computes a noise covariance matrix.
%
%  Input:
%   theta - P-dimensional column vector of hyperparameters; all
%           components must be nonnegative
%   M     - dimension of R
%
%  Output:
%   R     - Covariance matrix
%   dR    - Gradient of covariance matrix
%   invR  - Inverse of covariance matrix
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

% R matrix
R = theta(1)*speye(M,M);        % covariance associated with noise in data
invR = (1/theta(1))*speye(M,M); % inverse of covariance associated with noise in data
dR{1,1} = speye(M,M);           % gradient of R (w.r.t. theta1)
dR{2,1} = 0*speye(M,M);         % gradient of R (w.r.t. theta2)
dR{3,1} = 0*speye(M,M);         % gradient of R (w.r.t. theta3)

end

