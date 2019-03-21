function selected_defective_data = selected_def(data_feature,data_label)
defective_index = find(data_label==1);
def_feature = data_feature(defective_index,:);

normalize_def_feature = zscore(def_feature);
Y = normalize_def_feature';
instance_size = size(normalize_def_feature,1);
alpha = 1; 
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

selected_defective = def_feature(sInd,:);
num_selected_defective = size(selected_defective,1)
selected_defective_data = [selected_defective ones(num_selected_defective,1)];