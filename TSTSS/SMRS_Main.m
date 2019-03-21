function first_stage_selected = SMRS_Main(data)
data_feature = data(:,1:end-1);
data_label = data(:,end);

selected_defective_data = selected_def(data_feature,data_label);
selected_nondefective_data = selected_nondef(data_feature,data_label);

first_stage_selected = [selected_defective_data;selected_nondefective_data];