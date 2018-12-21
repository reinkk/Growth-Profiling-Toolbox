function [b,stats] = mq2wls(A,y,V,caller)

%MQ2WLS Weighted Least Squares
%   X = MQ2WLS(A,y) returns the ordinary least squares solution to the
%   linear system of equations A*b = y, i.e., b is the N-by-1 vector that
%   minimizes the sum of squared errors (y - A*b)'*(y - A*b), where A is
%   M-by-N, and y is M-by-1.  When rank(A) < N, MQ2WLS
%   sets the maximum possible number of elements of X to zero to obtain a
%   "basic solution".
%
%   X = MQ2WLS(A,y,W), where W is a vector length M of real positive weights,
%   returns the weighted least squares solution to the linear system A*b =
%   y, i.e., X minimizes (y - A*b)'*diag(W)*(y - A*b).  W typically
%   contains either counts or inverse variances.
%
%   See also MLDIVIDE, SLASH, LSQNONNEG, QR.
%
%   References:
%      [1] Paige, C.C. (1979) "Computer Solution and Perturbation
%          Analysis of Generalized Linear Leat Squares Problems",
%          Mathematics of Computation 33(145):171-183.
%      [2] Golub, G.H. and Van Loan, C.F. (1996) Matrix Computations,
%          3rd ed., Johns Hopkins University Press.
%      [3] Goodall, C.R. (1993) "Computation using the QR Decomposition",
%          in Computational Statistics, Vol. 9 of Handbook of Statistics,
%          edited by C.R. Rao, North-Holland, pp. 492-494.
%      [4] Strang, G. (1986) Introduction to Applied Mathematics,
%          Wellesley-Cambridge Press, pp. 398-399.
%
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 5.15.4.7 $  $Date: 2007/09/18 02:16:54 $

[nobs,nvar] = size(A); % num observations, num predictor variables
nrhs = size(y,2);      % num right-hand sides

% Assume V is full rank until we find out otherwise.  Decide later whether or not
% to factor V using Cholesky, unless it has been specified.
rankV = nobs;

% Error checking.
if size(y,1) ~= nobs
    error('mq2wls:InputSizeMismatch', 'y must have the same number of rows as A.');
end

% V not given, assume the identity.
if nargin < 3 || isempty(V)
    V = [];

% V given as a weight vector.
elseif isvector(V) && numel(V)==nobs && all(V>=0)

else
    error('mq2wls:InvalidWVec', '%s%d%s', ...
          'V must be a ',nobs,'-by-1 weight vector.');
end

outClass = superiorfloat(A,y,V);

if isempty(V) % No weights given, assume proportional to identity.
    V = ones(size(y));
    Aw = A;
    yw = y;
elseif isvector(V) % Weights given, scale rows of design matrix and response.
    D = sqrt(full(V(:)));
    Aw = bsxfun(@times,D,A);
    yw = bsxfun(@times,D,y);
end

% Factor the design matrix, incorporate covariances or weights into the
% system of equations, and transform the response vector.
[Q,R,perm] = qr(Aw,0);
z = Q'*yw;

% Use the rank-revealing QR to remove dependent columns of Aw.
keepCols = (abs(diag(R)) > abs(R(1)).*max(nobs,nvar).*eps(class(R)));
rankAw = sum(keepCols);

if rankAw < nvar
    warning('mq2wls:RankDefDesignMat', 'A is rank deficient to within machine precision.');
    R = R(keepCols,keepCols);
    z = z(keepCols,:);
    perm = perm(keepCols);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of R that were thrown out.
bb = R \ z;
b = zeros(nvar,nrhs,outClass);
b(perm,1:nrhs) = bb;

% Compute the MSE, need it for the std. errs. and covs.
dfe = nobs - rankAw;
if dfe > 0
    mse = full(sum(yw.*conj(yw),1) - sum(z.*conj(z),1)) ./ dfe;
else % rankAw == nobs, and so Ab == y exactly
    mse = zeros(1,nrhs,outClass);
end

% Compute the covariance matrix of the LS estimates, or just their
% SDs.  Fill in zeros corresponding to exact zero coefficients.
Rinv = R \ eye(rankAw,outClass);
S = zeros(nvar,nvar,outClass);
S(perm,perm) = Rinv*Rinv' .* mse; % guaranteed to be hermitian
stdb = sqrt(diag(S));

%Calculate statistics
yhat = A*b;
ybar = sum(V.*y)/sum(V);
residuals = y - yhat;
n = length(y);
p = length(b);
dfe = n-p;
dft = n-1;
if sum(A(:,end)) == size(A,1) %last column has only ones
    dfr = p - 1;
else
    dfr = p;
end
sse = norm(residuals)^2;    % sum of squared errors
ssr = norm(yhat - ybar)^2;  % regression sum of squares
sst = norm(y - ybar)^2;     % total sum of squares;
mse = sse./dfe;
h = sum(Q.*Q,2);
s_sqr_i = (dfe*mse - abs(residuals).^2./(1-h))./(dfe-1);


%Outputs
if isequal(caller,'robfit')
    stats = rankAw;
elseif isequal(caller,'initfit')
    stats.R = R;
    stats.perm = perm;
    stats.rankA = rankAw;
    stats.leverage = h;
    stats.res = y - A*b;
    stats.normr = norm(stats.res);
else
    stats.numobs = n;
    stats.numparam = p;
    stats.coeffs = b;
    stats.residuals = residuals;    
    stats.stdresid = residuals/sqrt(mse*(1-h));
    stats.studresid = residuals./sqrt(s_sqr_i.*(1-h));
    stats.s = sqrt(mse);
    %Degrees of freedom
    stats.dfr = dfr;
    stats.dfe = dfe;    
    stats.dft = dft;   
    %Sum of squares
    stats.ssr = ssr;
    stats.sse = sse;
    stats.sst = sst;
    %Goodness of fit
    stats.rsquare = 1 - sse/sst;
    stats.r2f = (sst/dft)/(sse/dfe);
    stats.pvalr2 = 1 - fcdf(stats.r2f, dft, dfe);
    stats.adjrsquare = 1 - (1-stats.rsquare)*(n-1)/dfe;    
    %t Statistics
    stats.se = stdb;
    ok_ = stats.se > 0;    %avoid divide by zero
    stats.t = NaN(size(b));
    stats.t(ok_) = b(ok_) ./ stats.se(ok_);
    stats.pvalt = 2 * tcdf(-abs(stats.t), dfe);
    %F statistics
    stats.f = (ssr/dfr)/(sse/dfe);
    stats.pvalf = 1 - fcdf(stats.f, dfr, dfe);
    %Least squares parameters
    stats.S = S;
    stats.Rinv = Rinv;
    stats.Jacobian = A;
end