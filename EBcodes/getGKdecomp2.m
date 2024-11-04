function [U, B, V] = getGKdecomp2(A, Q, R, b, maxiter,options)
% [U, B, V] = getGKdecomp2(A, Q, R, b, maxiter,options)
%
%  Performs maxiter steps of generalized Golub-Kahan bidiagonalization.
%
% Input:
%          A - matrix
%       Q, R - covariance matrices
%          b - right hand side
%    maxiter - number of genGKB iterations
%    options - structure from HyBR (see HyBRset)
%
% Output:
%       U, V - "orthonormal" (w.r.t. R and Q) matrices
%          B - bidiagonal matrix
%
%  Refs:
%   Chung and Saibaba. "Generalized Hybrid Iterative Methods for 
%       Large-Scale Bayesian Inverse Problems", SISC, 2017
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBR_lsmrget(options,'Reorth'), {'on'});

% Initialize for memory alocation
[m,n] = size(A);
if reorth
    U = zeros(m,maxiter+1);
    RinvU = zeros(m,maxiter+1);
else
    U = zeros(m,1);
    RinvU = zeros(m,1);
end
V = zeros(n,maxiter);
QV = zeros(n,maxiter);
B = zeros(maxiter+1,maxiter);

beta = normM(b,@(x)R\x);
U(:,1) = (1 / beta)*b;
RinvU(:,1) = R\U(:,1);

for k = 1:maxiter+1 %Iteration (k=1) is just an initialization

    if reorth % Need reorthogonalization
        if k == 1
            v = A'*(R\U(:,k));
        else
            % v = A'*(R\U(:,k)) - B(k, k-1)*V(:,k-1);
            v = A'*RinvU(:,k) - B(k, k-1)*V(:,k-1);
        end

        % Reorthogonalize V
        for j = 1:k-1
            % temp = (V(:,j)'*(Q*v));
            temp = v'*QV(:,j);
            v = v - temp*V(:,j);
        end
        alpha = normM(v,Q);
        v = v / alpha;
        QV(:,k) = Q*v;
        u = A*QV(:,k) - alpha*U(:,k);

        % Reorthogonalize U
        for j = 1:k
            % temp = (U(:,j)'*(R\u));
            temp = u'*RinvU(:,j);
            u = u - temp*U(:,j);
        end
        beta = normM(u, @(x)R\x);
        u = u / beta;
        U(:,k+1) = u;
        RinvU(:,k+1) = R \ u;
    else % Do not need reorthogonalization, save on storage
        if k == 1
            v = A'*(R\U(:));
        else
            v = A'*(R\U(:)) - B(k, k-1)*V(:,k-1);
        end
        alpha = normM(v,Q);
        v = v / alpha;
        Qv = Q*v;
        u = A*Qv - alpha*U(:);

        beta = normM(u, @(x)R\x);
        u = u / beta;
        U = u(:);
    end

    V(:,k) = v;
    B(k,k) = alpha;
    B(k+1,k) = beta;
end
end

function nrm = normM(v, M)
if isa(M, 'function_handle')
    Mv = M(v);
else
    Mv = M*v;
end
nrm = sqrt(v'*Mv);
end