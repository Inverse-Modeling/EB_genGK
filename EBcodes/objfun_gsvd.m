function [F, gradF, val] = objfun_gsvd(theta, inv, k)
%
% [F, gradF, val] = objfun_gsvd(theta, inv, k)
%
%   Evaluates the approximate objective function and gradient based on
%   gsvd, at parameters theta
%
% Input:
%        theta : P-dimensional column vector of hyperparameters; all
%           components must be nonnegative
%          inv : structure input with fields
%                  prior_type - prior for the hyperparameters (see Prior.m)
%                  dn - observation vector
%                  A - forward operator
%                  s_true - true solution
%                  Q - prior covariance
%                  genGK_iter - gsvd rank (if k not defined)
%             k : gsvd rank
%
% Output: F : approx objective function evaluated at theta
%         gradF : approx gradient evaluated at theta
%         val : structure output with fields
%                    F - [F1, F2] components of the objective function
%                    dF - components of the gradient bound
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

if ~exist('k','var')
  genGK_iter = inv.genGK_iter;
else
  genGK_iter = k;
end

[~, gradP, logP] = Prior(inv.prior_type);
dn = inv.dn;
A  = inv.A;
hyp_dim = length(theta);

[M,N] = size(A);
[R, dR] = R_mat(theta, M);
[mu, dmu] = mu_vec(theta, N);
[Qm, dQ] = Q_mat(inv.Q, theta);


%%% Objective Function Computation %%%
if isequal(class(A),'function_handle')
  A = A(eye([M,N]),'notransp');
end

R_half = sqrtm(R);

[U,S,V] = svd(R_half\(full(A)*sqrtm(Qm*eye(N))),'econ');
Ahat = U(:,1:genGK_iter)*(S(1:genGK_iter,1:genGK_iter)*V(:,1:genGK_iter)');
Id = eye(size(Ahat,1));

Z = R_half*( ( Ahat*Ahat.' + Id )*R_half );

F1 = (1/2)*logdetstable(Z);
res = A*mu - dn;
L = Z\res;
F2 = (1/2)*dot(res, L);

F = - logP(theta) + F1 + F2; % objective function
val.F = [F1;F2];

%%% Gradient Computation %%%
gradF = nan(hyp_dim,1);
ALPHA = nan(hyp_dim,1);
GAMMA = nan(hyp_dim,1);
DELTA = nan(hyp_dim,1);

for i = 1:hyp_dim
  gradZ_t_i = dR{i,1} + A*(dQ{i,1}*full(A).');
  k_i = Z\gradZ_t_i;

  alpha_i = (1/2)*trace(k_i);
  gamma_i = - (1/2)*( (res.')*(k_i*L) );
  delta_i = ( (L.')*(A*dmu{i,1}) );

  ALPHA(i,1) = alpha_i;
  GAMMA(i,1) = gamma_i;
  DELTA(i,1) = delta_i;

  gradF(i,1) = -gradP{i,1}(theta) + alpha_i + gamma_i + delta_i;
end

val.dF = [ALPHA,GAMMA,DELTA];

end