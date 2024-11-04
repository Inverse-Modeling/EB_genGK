function [x_optimal, fval, exitflag, output, xvals] = optimizer(objfunc, theta_0, inv, optim_procedure)
%
% [x_optimal, fval, exitflag, output, xvals] = optimizer(objfunc, theta_0, inv, optim_procedure)
%
% Input:
%     objfunc : objective function handle
%     theta_0 : initial guess of hyperparameters
%     inv     : structure of problem parameters
%     optim_procedure : algorithm for fmincon
%
% Output:
%     x_optimal : optimal parameters
%     fval : corresponding function value
%     exitflag : reason for exiting 
%     output : output parameters from fmincon
%     xvals : values of x where fun has been evaluated
%
% Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)

xvals = [];  % This will contain the values of x where fun has been evaluated
fun = @(theta) objfunc(theta, inv);
Aa = [];
Aeq = [];
bb = [];
beq = [];
hyp_dim = length(theta_0);
lb = zeros(hyp_dim,1);
ub = Inf*ones(hyp_dim,1); % original upper bound
nonlcon = [];

options = optimoptions(@fmincon,'Algorithm',optim_procedure,...
  'MaxIterations', 1e8, ...
  'MaxFunctionEvaluations', 1e8, ...
  'SpecifyObjectiveGradient', true, ...
  'FiniteDifferenceType', 'central', ...
  'Display','iter',...
  'PlotFcn', 'optimplotfval',...
  'CheckGradients', false,...
  'OutputFcn', @outfun);

[x_optimal,fval,exitflag,output] = fmincon(fun,theta_0,Aa,bb,Aeq,beq,lb,ub,nonlcon,options);
  
  function stop = outfun(x,optimValues,state)
    stop = false;
    if isequal(state,'iter')
      xvals = [xvals, x];
    end
  end
end