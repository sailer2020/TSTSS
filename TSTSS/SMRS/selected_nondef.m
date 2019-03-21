function selected_nondefective_data = selected_nondef(data_feature,data_label)
nondefective_index = find(data_label==0);
nondef_feature = data_feature(nondefective_index,:);

normalize_nondef_feature = zscore(nondef_feature);
Y = normalize_nondef_feature';
instance_size = size(normalize_nondef_feature,1);

alpha = 1; % regularization parameter
r = 0;

verbose = true;
Lambda = initLambda(Y);
Lambda = floor(Lambda/20)*20;

for Lambda = Lambda:5:1000000
    sInd = [];
    [sInd,C] = smrs(Y,alpha,r,verbose,Lambda);
    select_number = size(sInd,2);
    if select_number >= instance_size*0.8
        break
    end
end

selected_nondefective = nondef_feature(sInd,:);
num_selected_nondefective = size(selected_nondefective,1)
selected_nondefective_data = [selected_nondefective zeros(num_selected_nondefective,1)];