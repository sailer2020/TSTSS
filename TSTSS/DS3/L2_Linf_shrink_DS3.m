function x = L2_Linf_shrink_DS3(y,t)

x = y;
[dummy,o] = sort(abs(y),'descend');
z = y(o);
mz = abs(z);

% find cut-off index
cs = cumsum(abs(z(1:length(z)-1)))./(1:length(z)-1)'-t./(1:length(z)-1)';
d = (cs>abs(z(2:length(z))));
if sum(d) == 0
   cut_index = length(y);
else
   cut_index = min(find(d==1));
end

% shrink coordinates 1 to cut_index
zbar = mean(abs(z(1:cut_index)));
if cut_index < length(y)
   x(o(1:cut_index)) = sign(z(1:cut_index))*max(zbar-t/cut_index,abs(z(cut_index+1)));
else
   x(o(1:cut_index)) = sign(z(1:cut_index))*max(zbar-t/cut_index,0);
end
