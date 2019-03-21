clc, clear all

addpath('.\liblinear\');
addpath('.\dataset\');
addpath('.\Classifiers\');
addpath('.\DS3\');
addpath('.\SMRS\');

X_data = load('*****'); %load the data of the prior version

Y_data = load('****'); %load the data of the current version


first_stage_selected = SMRS_Main(X_data);

[s_Ind20, s_Ind30, s_Ind40, s_Ind50, s_Ind60, s_Ind70, s_Ind80, s_Ind90] = DS3_Main(first_stage_selected,Y_data);

% after obtaining the selected modules, you can use any kind of machine learning classifier to construct CVDP model.

