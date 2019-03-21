function [sInd, C] = smrs(Y,alpha,r,verbose,Lambda)

if (nargin < 2)
    alpha = 5;
end
if (nargin < 3)
    r = 0;
end
if (nargin < 4)
    verbose = true;
end

q = 2;
regParam = [alpha alpha];
affine = true;
thr = 1 * 10^-7;
maxIter = 5000;
thrS = 0.99; 
thrP = 0.95;
N = size(Y,2);

Y = Y - repmat(mean(Y,2),1,N);
if (r >= 1)
    [~,S,V] = svd(Y,0);
    r = min(r,size(V,1));
    Y = S(1:r,1:r) * V(:,1:r)';
end

C = almLasso_mat_func(Y,affine,regParam,q,thr,maxIter,verbose,Lambda);

sInd = findRep(C,thrS,q);