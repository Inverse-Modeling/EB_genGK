function [F, gradF, val] = objfun_gengk(theta, inv, k)
%
% [F, gradF, val] = objfun_gengk(theta, inv, k)
%
%   Evaluates the approximate objective function and gradient based on
%   genGK, at parameters theta
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
%                  genGK_iter - number of genGK iterations (if k not
%                  defined)
%             k : bidiagonalization parameter
%
% Output: F : approx objective function evaluated at theta
%         gradF : approx gradient evaluated at theta
%         val : structure output with fields
%                    F - [F1, F2] components of the objective function
%                    ZinvRes
%                    dF - components of the gradient
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

if ~exist('k','var')
  genGK_iter = inv.genGK_iter;
else
  genGK_iter = k;
end

A  = inv.A;
M = inv.M;
N = inv.N;
[R, dR] = R_mat(theta, M);
[mu, dmu] = mu_vec(theta, N);
[Qm, dQ] = Q_mat(inv.Q, theta);

[~, gradP, logP] = Prior(inv.prior_type);
dn = inv.dn;

hyp_dim = length(theta);
s_true   = inv.s_true;

if isequal(class(A),'function_handle')
  A = funMat(@ (x) A(x, 'notransp'), @(x) A(x, 'transp'), A([],'size'));
  M = inv.M;
  N = inv.N;
end

input = HyBR_lsmrset('RegPar', 1,'x_true', s_true,'Iter',genGK_iter, 'Reorth', 'on');
[Uk, Bk, Vk] = getGKdecomp2(A, Qm, R, dn(:), genGK_iter, input);

%%% objective function computation %%%
Idkp1 = speye(min(size(Bk)) + 1);
res = A*mu - dn;
beta1 = sqrt(dot(res, R\res));
e1 = Idkp1(:,1);

[uB,sB,vB] = svd(Bk,0);
ss = diag(sB).^2;
Dk = Idkp1 + uB*diag(ss)*uB';


try
  F1_k = (1/2)*( logdetstable(R) + sum(log(1 + ss)) );
catch
  msg = ['Idkp1 + Bk*(Bk.'') is not positive definite; occurs at theta = (',...
    num2str(theta(1)),',',num2str(theta(2)),',',num2str(theta(3)),')'];
  error(msg);
end
F2_k = (1/2)*dot(beta1*e1, Dk\(beta1*e1));

F = -logP(theta) + F1_k + F2_k; % objective function
val.F = [F1_k;F2_k];
val.ZinvRes = Dk\(beta1*e1);


%%% gradient computation %%%
gradF = nan(hyp_dim,1);
ALPHA = nan(hyp_dim,1);
GAMMA = nan(hyp_dim,1);
DELTA = nan(hyp_dim,1);

Ek = vB*diag(ss)*vB'; % Bk.'*Bk;

Tk = speye(min(size(Bk))) + Ek;
Ck = Uk*Bk;
Jk = R\res - ((R\Ck)*(Tk\((R\Ck).'*res)));
for i = 1:hyp_dim
  Gk_i = Vk.'*(dQ{i,1}*Vk);
  Sk_i = Ek*Gk_i;

  alpha_i = (1/2)*(trace(Sk_i) + trace(R\dR{i,1}) ...
    - trace(Tk\(Sk_i*Ek)) ...
    - trace(( Tk\((Ck.')*(R\dR{i,1})) )*( R\Ck )));

  gamma_i = -(1/2)*((Ck.'*Jk)'*( Gk_i*(Ck.'*Jk) ) ) ...
    -(1/2)*((Jk.')*(dR{i,1}*Jk));

  delta_i = ( dot(Jk, A*dmu{i,1} ) );

  ALPHA(i,1) = alpha_i;
  GAMMA(i,1) = gamma_i;
  DELTA(i,1) = delta_i;

  gradF(i,1) = -gradP{i,1}(theta) + alpha_i + gamma_i + delta_i;
end

val.dF = [ALPHA,GAMMA,DELTA];

end
