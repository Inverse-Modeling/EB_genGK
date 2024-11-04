function ld = logdetstable(A, varargin)
%  ld = logdetstable(A, varargin)
%
% A stable way to compute the log-determinant
%
% Input: A - matrix
%        additional parameters to identify if triangular
%           tri - true if A is a triagnular matrix
%                 false if A is not triagnular (need to construct Cholesky
%                 factorization)
% Output: ld - logdet(A)
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

if nargin > 1
  tri = varargin{1};
else
  tri = false;
end

%Compute cholesky factorization only if not triangular
if ~tri
  R = chol(A,'upper');
  ld = 2*sum(log(diag(R)));
else
  ld = sum(log(diag(A)));
end


end