function [err_bound, val] = mc_err_bnd(theta, inv, n_mc)
%
% [err_bound, val] = mc_err_bnd(theta, inv, n_mc)
%
%  This function computes the error bound at parameters theta
%
% Input:
%   theta : P-dimensional column vector of hyperparameters; all
%           components must be nonnegative
%     inv : structure input with fields
%                  prior_type - prior for the hyperparameters (see Prior.m)
%                  dn - observation vector
%                  A - forward operator
%                  s_true - true solution
%                  Q - prior covariance
%   n_mc : number of Monte Carlo samples
%
% Output:
%   err_bound : error bound
%   val : output structure with fields
%           LD - log determinant term
%           QU - quadratic term
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

[M,N] = size(inv.A);
[Q,~] = Q_mat(inv.Q,theta);
[R,~] = R_mat(theta,M);

A = inv.A;
genGK_iter = inv.genGK_iter;

input = HyBR_lsmrset('RegPar', 1,'x_true', inv.s_true,'Iter', ...
  inv.genGK_iter, 'Reorth', 'on');
[~,Bk, Vk] = getGKdecomp2(A, Q, R, inv.dn, genGK_iter, input);

Z = randn([N,n_mc]);
Theta = zeros(genGK_iter,1);
K = @(x) A'*(R\(A*x));
QZ = Q*Z;
KQZ = K(QZ);

for l = 1:genGK_iter
  T = Bk(1:l+1,1:l)'*Bk(1:l+1,1:l);
  v = Vk(:,1:l);
  Theta(l,1) = trace(Z'*( KQZ - v*(T*(v'*QZ))))/n_mc;
end

dn = inv.dn;
[mu, ~] = mu_vec(theta, N);
res = A*mu - dn;
beta1 = sqrt(dot(res, R\res));
err_bound = (1/2)*(( Theta) + (beta1^2)*Theta./(1 + Theta));
val.LD = (1/2)*( Theta );
val.QU = (1/2)*( (beta1^2)*Theta./(1 + Theta) );
end
