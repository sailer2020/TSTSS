function [C2,Err] = almLasso_mat_func(Y,affine,alpha,q,thr,maxIter,verbose,Lambda)


if (nargin < 2)
    % default subspaces are linear
    affine = false; 
end
if (nargin < 3)
    % default regularizarion parameters
    alpha = 5;
end
if (nargin < 4)
    % default norm in L1/Lq optimization program
    q = 2;
end
if (nargin < 5)
    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    thr = 1*10^-7; 
end
if (nargin < 6)
    % default maximum number of iterations of ALM
    maxIter = 5000; 
end
if (nargin < 7)
    % reporting iterations and errors
    verbose = true; 
end


if (length(alpha) == 1)
    alpha1 = alpha(1);
    alpha2 = alpha(1);
elseif (length(alpha) == 2)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
end

if (length(thr) == 1)
    thr1 = thr(1);
    thr2 = thr(1);
elseif (length(thr) == 2)
    thr1 = thr(1);
    thr2 = thr(2);
end

[D,N] = size(Y);

mu1p = alpha1 * 1/Lambda;
mu2p = alpha2 * 1;

if (~affine)
    mu1 = mu1p;
    mu2 = mu2p;
    P = Y'*Y;
    A = inv(mu1.*P+mu2.*eye(N));
    C1 = zeros(N,N);
    Lambda2 = zeros(N,N);
    err1 = 10*thr1; 
    i = 1;
    while ( err1 > thr1 && i < maxIter )
        Z = A * (mu1.*P+mu2.*C1-Lambda2);
        C2 = shrinkL1Lq(Z+Lambda2./mu2,1/mu2,q);
        Lambda2 = Lambda2 + mu2 .* (Z - C2);
        err1 = errorCoef(Z,C2);
        C1 = C2;
        i = i + 1;
    end
    Err = err1;

else
    mu1 = mu1p;
    mu2 = mu2p;
    P = Y'*Y;
    A = inv(mu1.*P+mu2.*eye(N)+mu2.*ones(N,N));
    C1 = zeros(N,N);
    Lambda2 = zeros(N,N);
    lambda3 = zeros(1,N);
    err1 = 10*thr1; err2 = 10*thr2;
    i = 1;
    % ALM iterations
    while ( (err1 > thr1 || err2 > thr1) && i < maxIter )
        % updating Z
        Z = A * (mu1.*P+mu2.*(C1-Lambda2./mu2)+mu2.*ones(N,N)+repmat(lambda3,N,1));
        % updating C
        C2 = shrinkL1Lq(Z+Lambda2./mu2,1/mu2,q);  
        % updating Lagrange multipliers
        Lambda2 = Lambda2 + mu2 .* (Z - C2);
        lambda3 = lambda3 + mu2 .* (ones(1,N) - sum(Z,1));
        % computing errors
        err1 = errorCoef(Z,C2);
        err2 = errorCoef(sum(Z,1),ones(1,N));

        C1 = C2;
        i = i + 1;
    end
    Err = [err1;err2];
end
