function [x, resnorm, F, exitflag, output, jacob] = newtonraphson(fun, x0, options)
% NEWTONRAPHSON Solve set of non-linear equations using Newton-Raphson method.
%
% [X, RESNORM, F, EXITFLAG, OUTPUT, JACOB] = NEWTONRAPHSON(FUN, X0, OPTIONS)
% FUN is a function handle that returns a vector of residuals equations, F,
% and takes a vector, x, as its only argument. When the equations are
% solved by x, then F(x) == zeros(size(F(:), 1)).
%
% Optionally FUN may return the Jacobian, Jij = dFi/dxj, as an additional
% output. The Jacobian must have the same number of rows as F and the same
% number of columns as x. The columns of the Jacobians correspond to d/dxj and
% the rows correspond to dFi/d.
%
%   EG:  J23 = dF2/dx3 is the 2nd row ad 3rd column.
%
% If FUN only returns one output, then J is estimated using a center
% difference approximation,
%
%   Jij = dFi/dxj = (Fi(xj + dx) - Fi(xj - dx))/2/dx.
%
% NOTE: If the Jacobian is not square the system is either over or under
% constrained.
%
% X0 is a vector of initial guesses.
%
% OPTIONS is a structure of solver options created using OPTIMSET.
% EG: options = optimset('TolX', 0.001).
%
% The following options can be set:
% * OPTIONS.TOLFUN is the maximum tolerance of the norm of the residuals.
%   [1e-6]
% * OPTIONS.TOLX is the minimum tolerance of the relative maximum stepsize.
%   [1e-6]
% * OPTIONS.MAXITER is the maximum number of iterations before giving up.
%   [100]
% * OPTIONS.DISPLAY sets the level of display: {'off', 'iter'}.
%   ['iter']
%
% X is the solution that solves the set of equations within the given tolerance.
% RESNORM is norm(F) and F is F(X). EXITFLAG is an integer that corresponds to
% the output conditions, OUTPUT is a structure containing the number of
% iterations, the final stepsize and exitflag message and JACOB is the J(X).
%
% See also OPTIMSET, OPTIMGET, FMINSEARCH, FZERO, FMINBND, FSOLVE, LSQNONLIN
%
% References:
% * http://en.wikipedia.org/wiki/Newton's_method
% * http://en.wikipedia.org/wiki/Newton's_method_in_optimization
% * 9.7 Globally Convergent Methods for Nonlinear Systems of Equations 383,
%   Numerical Recipes in C, Second Edition (1992),
%   http://www.nrbook.com/a/bookcpdf.php

% Version 0.5
% * allow sparse matrices, replace cond() with condest()
% * check if Jstar has NaN or Inf, return NaN or Inf for cond() and return
%   exitflag: -1, matrix is singular.
% * fix bug: max iteration detection and exitflag reporting typos
% Version 0.4
% * allow lsq curve fitting type problems, IE non-square matrices
% * exit if J is singular or if dx is NaN or Inf
% Version 0.3
% * Display RCOND each step.
% * Replace nargout checking in funwrapper with ducktypin.
% * Remove Ftyp and F scaling b/c F(typx)->0 & F/Ftyp->Inf!
% * User Numerical Recipies minimum Newton step, backtracking line search
%   with alpha = 1e-4, min_lambda = 0.1 and max_lambda = 0.5.
% * Output messages, exitflag and min relative step.
% Version 0.2
% * Remove `options.FinDiffRelStep` and `options.TypicalX` since not in MATLAB.
% * Set `dx = eps^(1/3)` in `jacobian` function.
% * Remove `options` argument from `funwrapper` & `jacobian` functions
%   since no longer needed.
% * Set typx = x0; typx(x0==0) = 1; % use initial guess as typx, if 0 use 1.
% * Replace `feval` with `evalf` since `feval` is builtin.

%% initialize
% There are no argument checks!
x0 = x0(:); % needs to be a column vector
% set default options
oldopts = optimset( ...
    'TolX', 1e-12, 'TolFun', 1e-6, 'MaxIter', 100, 'Display', 'iter');
if nargin<3
    options = oldopts; % use defaults
else
    options = optimset(oldopts, options); % update default with user options
end
FUN = @(x)funwrapper(fun, x); % wrap FUN so it always returns J
%% get options
TOLX = optimget(options, 'TolX'); % relative max step tolerance
TOLFUN = optimget(options, 'TolFun'); % function tolerance
MAXITER = optimget(options, 'MaxIter'); % max number of iterations
DISPLAY = strcmpi('iter', optimget(options, 'Display')); % display iterations
TYPX = max(abs(x0), 1); % x scaling value, remove zeros
ALPHA = 1e-4; % criteria for decrease
MIN_LAMBDA = 0.1; % min lambda
MAX_LAMBDA = 0.5; % max lambda
%% set scaling values
% TODO: let user set weights
weight = ones(numel(FUN(x0)),1);
J0 = weight*(1./TYPX'); % Jacobian scaling matrix
%% set display
if DISPLAY
    fprintf('\n%10s %10s %10s %10s %10s %12s\n', 'Niter', 'resnorm', 'stepnorm', ...
        'lambda', 'rcond', 'convergence')
    for n = 1:67,fprintf('-'),end,fprintf('\n')
    fmtstr = '%10d %10.4g %10.4g %10.4g %10.4g %12.4g\n';
    printout = @(n, r, s, l, rc, c)fprintf(fmtstr, n, r, s, l, rc, c);
end
%% check initial guess
x = x0; % initial guess
[F, J] = FUN(x); % evaluate initial guess
Jstar = J./J0; % scale Jacobian
if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
    exitflag = -1; % matrix may be singular
else
    exitflag = 1; % normal exit
end
if issparse(Jstar)
    rc = 1/condest(Jstar);
else
    if any(isnan(Jstar(:)))
        rc = NaN;
    elseif any(isinf(Jstar(:)))
        rc = Inf;
    else
        rc = 1/cond(Jstar); % reciprocal condition
    end
end
resnorm = norm(F); % calculate norm of the residuals
dx = zeros(size(x0));convergence = Inf; % dummy values
%% solver
Niter = 0; % start counter
lambda = 1; % backtracking
if DISPLAY,printout(Niter, resnorm, norm(dx), lambda, rc, convergence);end
while (resnorm>TOLFUN || lambda<1) && exitflag>=0 && Niter<=MAXITER
    if lambda==1
        %% Newton-Raphson solver
        Niter = Niter+1; % increment counter
        dx_star = -Jstar\F; % calculate Newton step
        % NOTE: use isnan(f) || isinf(f) instead of STPMAX
        dx = dx_star.*TYPX; % rescale x
        g = F'*Jstar; % gradient of resnorm
        slope = g*dx_star; % slope of gradient
        fold = F'*F; % objective function
        xold = x; % initial value
        lambda_min = TOLX/max(abs(dx)./max(abs(xold), 1));
    end
    if lambda<lambda_min
        exitflag = 2; % x is too close to XOLD
        break
    elseif any(isnan(dx)) || any(isinf(dx))
        exitflag = -1; % matrix may be singular
        break
    end
    x = xold+dx*lambda; % next guess
    [F, J] = FUN(x); % evaluate next residuals
    Jstar = J./J0; % scale next Jacobian
    f = F'*F; % new objective function
    %% check for convergence
    lambda1 = lambda; % save previous lambda
    if f>fold+ALPHA*lambda*slope
        if lambda==1
            lambda = -slope/2/(f-fold-slope); % calculate lambda
        else
            A = 1/(lambda1 - lambda2);
            B = [1/lambda1^2,-1/lambda2^2;-lambda2/lambda1^2,lambda1/lambda2^2];
            C = [f-fold-lambda1*slope;f2-fold-lambda2*slope];
            coeff = num2cell(A*B*C);
            [a,b] = coeff{:};
            if a==0
                lambda = -slope/2/b;
            else
                discriminant = b^2 - 3*a*slope;
                if discriminant<0
                    lambda = MAX_LAMBDA*lambda1;
                elseif b<=0
                    lambda = (-b+sqrt(discriminant))/3/a;
                else
                    lambda = -slope/(b+sqrt(discriminant));
                end
            end
            lambda = min(lambda,MAX_LAMBDA*lambda1); % minimum step length
        end
    elseif isnan(f) || isinf(f)
        % limit undefined evaluation or overflow
        lambda = MAX_LAMBDA*lambda1;
    else
        lambda = 1; % fraction of Newton step
    end
    if lambda<1
        lambda2 = lambda1;f2 = f; % save 2nd most previous value
        lambda = max(lambda,MIN_LAMBDA*lambda1); % minimum step length
        continue
    end
    %% display
    resnorm0 = resnorm; % old resnorm
    resnorm = norm(F); % calculate new resnorm
    convergence = log(resnorm0/resnorm); % calculate convergence rate
    stepnorm = norm(dx); % norm of the step
    if any(isnan(Jstar(:))) || any(isinf(Jstar(:)))
        exitflag = -1; % matrix may be singular
        break
    end
    if issparse(Jstar)
        rc = 1/condest(Jstar);
    else
        rc = 1/cond(Jstar); % reciprocal condition
    end
    if DISPLAY,printout(Niter, resnorm, stepnorm, lambda1, rc, convergence);end
end
%% output
output.iterations = Niter; % final number of iterations
output.stepsize = dx; % final stepsize
output.lambda = lambda; % final lambda
if Niter>=MAXITER
    exitflag = 0;
    output.message = 'Number of iterations exceeded OPTIONS.MAXITER.';
elseif exitflag==2
    output.message = 'May have converged, but X is too close to XOLD.';
elseif exitflag==-1
    output.message = 'Matrix may be singular. Step was NaN or Inf.';
else
    output.message = 'Normal exit.';
end
jacob = J;
end

function [F, J] = funwrapper(fun, x)
% if nargout<2 use finite differences to estimate J
try
    [F, J] = fun(x);
catch
    F = fun(x);
    J = jacobian(fun, x); % evaluate center diff if no Jacobian
end
F = F(:); % needs to be a column vector
end

function J = jacobian(fun, x)
% estimate J
dx = eps^(1/3); % finite difference delta
nx = numel(x); % degrees of freedom
nf = numel(fun(x)); % number of functions
J = zeros(nf,nx); % matrix of zeros
for n = 1:nx
    % create a vector of deltas, change delta_n by dx
    delta = zeros(nx, 1); delta(n) = delta(n)+dx;
    dF = fun(x+delta)-fun(x-delta); % delta F
    J(:, n) = dF(:)/dx/2; % derivatives dF/d_n
end
end
