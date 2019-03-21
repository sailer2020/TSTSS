function C2 = shrinkL1Lp_DS3(C1,lambda,p)

C2 = [];
if ~isempty(lambda)
    [D,N] = size(C1);
    if (p == 1)
        C2 = max(abs(C1)-repmat(lambda,1,N),0) .* sign(C1);
    elseif (p == 2)

        r = max(sqrt(sum(C1.^2,2)) - lambda,0);
        C2 = repmat(r./(r+lambda),1,N) .* C1;
    elseif(p == inf)
        C2 = zeros(D,N);
        for j = 1:D
            C2(j,:) = L2_Linf_shrink_DS3(C1(j,:)',lambda(j))';
        end
    end
end