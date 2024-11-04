function [x_recov,flag,relres,iter,resvec] = pcgRecovery(A,Q,R,mu,rhs,tol,maxit)
%
% [x_recov,flag,relres,iter,resvec] = pcgRecovery(A,Q,R,mu,b,tol,maxit)
%
% For fixed set of hyperparameters, use preconditioned conjugate gradient
% to solve for the reconstructed image.
%
% Input:
%       A - forward operator
%       Q - prior covariance matrix
%       R - noise covariance matrix
%       mu - mean vector
%       rhs - right hand side given by
%               Q * A' * R^{-1} * b
%             where b is the observation vector
%       tol - tolerance
%       maxit - maximum number of iterations
%
% Output:
%       x_recov - reconstructed solution
%       flag, relres, iter, resvec - outputs from pcg
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

% function handle specifying left-hand-side matrix for pcg
fcn = @(x) lhs_pcg(x,A,Q,R);

%%% use pcg to solve for Q\x; then solve for x %%%
[Qinv_x_recov,flag,relres,iter,resvec] = pcg(fcn, rhs,tol,maxit);
x_recov = mu + Q * Qinv_x_recov;

end

function A_pcg = lhs_pcg(x,A,Q,R)
% A_pcg = lhs_pcg(x,mu,A,Q,R)  Left-hand-side Ahat matrix in matrix equation
%                             Ahat*(Q\x) = rhs
%
% A function handle specifying the matrix-vector product
%           (Q * A' * R^{-1} * A * Q + Q) * x
%

y = Q*x;
A_pcg = Q*(A.'*(R\(A*y))) + y;
end

