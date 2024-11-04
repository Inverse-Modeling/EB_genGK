function [F, gradF, val] = objfun_full(theta, inv)
%
% [F, gradF, val] = objfun_full(theta, inv)
%
%   Evaluates the full objective function and gradient at parameters theta
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
%
% Output: F : objective function evaluated at theta
%         gradF : gradient evaluated at theta
%         val : structure output with fields
%                    F - [F1, F2] components of the objective function
%                    ZinvRes
%                    Z
%                    dF - components of the gradient
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

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
Z = R + A*(Qm*(full(A).'));

F1 = (1/2)*logdetstable(Z);
res = A*mu - dn;
L = Z\res;
F2 = (1/2)*dot(res, L);

F = - logP(theta) + F1 + F2; % objective function
val.F = [F1;F2];
val.ZinvRes = L;
val.Z = Z;

%%% Gradient Computation %%%

gradF = nan(hyp_dim,1);
ALPHA = nan(hyp_dim,1);
GAMMA = nan(hyp_dim,1);
DELTA = nan(hyp_dim,1);

for i = 1:hyp_dim
  gradZ_t_i = dR{i,1} + A*(dQ{i,1}*full(A).');
  k_i = Z\gradZ_t_i;

  alpha_i = (1/2)*trace(k_i); % trace(Z^(-1)*dZ)
  gamma_i = - (1/2)*( (res.')*(k_i*L) ); % (A*mu - dn).'*(Z^(-1)*dZ)*(Z^(-1)*(A*mu - dn))
  delta_i = ( (L.')*(A*dmu{i,1}) ); % (Z^(-1)*(A*mu - dn)).'*A*dmu

  ALPHA(i,1) = alpha_i;
  GAMMA(i,1) = gamma_i;
  DELTA(i,1) = delta_i;

  gradF(i,1) = -gradP{i,1}(theta) + alpha_i + gamma_i + delta_i;
end

val.dF = [ALPHA,GAMMA,DELTA];

end