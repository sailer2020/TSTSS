function lambda = initLambda(Y)

[~,N] = size(Y);


T = zeros(N,1);
for i = 1:N
    yi = Y(:,i);
    ymean = mean(Y,2);
    T(i) = norm(yi'*(ymean*ones(1,N)-Y));
end
lambda = max(T);