function cssInd = findRep(C,thr,q)

if (nargin < 4)
    q = 2;
end
if (nargin < 3)
    thr = 0.9;
end

N = size(C,1);

r = zeros(1,N);
for i = 1:N
    r(i) = norm(C(i,:),q);
end
[nrm,nrmInd] = sort(r,'descend');
nrmSum = 0;
for j = 1:N
    nrmSum = nrmSum + nrm(j);
    if (nrmSum / sum(nrm) > thr)
        break;
    end
end
cssInd = nrmInd(1:j);