function C2 = ds3solver_regularized(D,p,options)

warning off

if (nargin < 2)
    p = inf;
end

if (nargin < 3)
    [rho_min, rho_max] = computeRegularizer(D,p);
    options.rho = .1 * rho_max;
    options.mu = 1 * 10^-1;
    options.maxIter = 3000;
    options.errThr = 1 * 10^-7;
    options.verbose = true;
end

rho = options.rho;
mu = options.mu;
maxIter = options.maxIter;
errThr = options.errThr;
verbose = options.verbose;
CFD = options.cfd;
ratio = 0.1;

% initialization
[Nr,Nc] = size(D);
terminate = false;
k = 1;
[~,idx] = min(sum(D,2));
C1 = zeros(size(D));
C1(idx,:) = 1;
Lambda = zeros(Nr,Nc);

% running iterations
while (~terminate)
    Z = shrinkL1Lp_DS3(C1-(Lambda+D)./mu,rho/mu*CFD,p);
    C2  = solver_BCLS_closedForm(Z + Lambda ./ mu);
    
    Lambda = Lambda + mu .* (Z - C2);
    
    err1 = errorCoef_DS3(Z,C2);
    err2 = errorCoef_DS3(C1,C2);
    
    if ( k >= maxIter || (err1 <= errThr && err2 <= errThr) )
        terminate = true;
        if (verbose)
            fprintf('Terminating: \n');
            fprintf('||Z-C||= %1.2e, ||C1-C2||= %1.2e, repNum = %3.0f, iteration = %.0f \n\n',err1,err2,length(findRepresentatives(C2)),k);
        end
    else
        k = k + 1;
        if (verbose)
            if (mod(k,100)==0)
                fprintf('||Z-C||= %1.2e, ||C1-C2||= %1.2e, repNum = %3.0f, iteration = %.0f \n',err1,err2,length(findRepresentatives(C2)),k);
            end
        end
    end
    C1 = C2;
end
