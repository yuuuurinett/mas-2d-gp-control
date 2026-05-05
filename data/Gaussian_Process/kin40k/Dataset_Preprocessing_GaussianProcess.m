%% Dataset_Preprocessing_GaussianProcess
clc; clear; close all;
%% KIN4K
clear;
DataSetRawStruct = load('KIN40K_train.mat');
X = DataSetRawStruct.x;
Y = DataSetRawStruct.y;
N_pretrain = 10000;
NrSet_pretrain = randperm(size(X,1),min(N_pretrain,size(X,1)));
gp = fitrgp(X(NrSet_pretrain,:),Y(NrSet_pretrain,1), ...
	'KernelFunction','ardsquaredexponential','Standardize',false);
SigmaL = gp.KernelInformation.KernelParameters(1:end-1);
SigmaF = gp.KernelInformation.KernelParameters(end);
SigmaN = gp.Sigma;
save('Dataset\Gaussian Process\KIN40K\KIN40K_Hyperparameter.mat', ...
	'SigmaF','SigmaL','SigmaN');
%% POL
clear;
DataSetRawStruct = load('Dataset\Gaussian Process\POL\POL_train.mat');
X = DataSetRawStruct.x;
Y = DataSetRawStruct.y;
N_pretrain = 10000;
NrSet_pretrain = randperm(size(X,1),min(N_pretrain,size(X,1)));
gp = fitrgp(X(NrSet_pretrain,:),Y(NrSet_pretrain,1), ...
	'KernelFunction','ardsquaredexponential','Standardize',false);
SigmaL = gp.KernelInformation.KernelParameters(1:end-1);
SigmaF = gp.KernelInformation.KernelParameters(end);
SigmaN = gp.Sigma;
save('Dataset\Gaussian Process\POL\POL_Hyperparameter.mat', ...
	'SigmaF','SigmaL','SigmaN');
%% POL
clear;
DataSetRawStruct = load('Dataset\Gaussian Process\PUMADYN32NM\PUMADYN32NM_train.mat');
X = DataSetRawStruct.x;
Y = DataSetRawStruct.y;
N_pretrain = 10000;
NrSet_pretrain = randperm(size(X,1),min(N_pretrain,size(X,1)));
gp = fitrgp(X(NrSet_pretrain,:),Y(NrSet_pretrain,1), ...
	'KernelFunction','ardsquaredexponential','Standardize',false);
SigmaL = gp.KernelInformation.KernelParameters(1:end-1);
SigmaF = gp.KernelInformation.KernelParameters(end);
SigmaN = gp.Sigma;
save('Dataset\Gaussian Process\PUMADYN32NM\PUMADYN32NM_Hyperparameter.mat', ...
	'SigmaF','SigmaL','SigmaN');