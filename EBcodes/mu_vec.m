function [mu, dmu] = mu_vec(theta, N)
% [mu, dmu] = mu_vec(theta, N)
%
% This function computes a mean vector and its derivative
%
% Input:
%   theta : P-dimensional column vector of hyperparameters; all
%           components must be nonnegative
%   N     : dimension of mu
%
% Output:
%   mu   - Mean vector (here we take it to be zero)
%   dmu  - Gradient of mean vector
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

mu = zeros(N,1);        % mu vector
dmu{1,1} = 0*sparse(N,1); %zeros(N,1);  % gradient of mu (w.r.t. theta1)
dmu{2,1} = 0*sparse(N,1); %zeros(N,1);  % gradient of mu (w.r.t. theta2)
dmu{3,1} = 0*sparse(N,1); %zeros(N,1);  % gradient of mu (w.r.t. theta3)

end

