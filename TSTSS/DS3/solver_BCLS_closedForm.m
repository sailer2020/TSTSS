function C  = solver_BCLS_closedForm(U)

[m,N] = size(U);
V = sort(U,'descend');
activeSet = 1:N;
theta = zeros(1,N);
i = 1;
while (~isempty(activeSet) && i <= m)
    idx = ( V(i,activeSet) - (sum(V(1:i,activeSet),1) - 1) / i ) <= 0;
    theta(activeSet(idx)) = (sum(V(1:i-1,activeSet(idx)),1) - 1) / (i-1);
    activeSet(idx) = [];
    i = i + 1;
end
if ~isempty(activeSet)
    theta(activeSet) = (sum(V(1:m,activeSet),1) - 1) / m;
end

C = max(U-repmat(theta,m,1),0);