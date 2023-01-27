function [trained_sas, stateCollection] = ...
    train_sas(trainInput, trainTeacher , sas, nForgetPoints, varargin)
% TRAIN_SAS Trains the output weights of a SAS for forecasting and
% classification tasks
% In the offline case, it computes the weights using the method
% esn.methodWeightCompute(for example 'gaussian_elimination' or 'pseudoinverse')
%
% inputs:
% trainInput = for forecasting this is the input time series of size
% nInputUnits x nTrainingPoints.
% for classification tasks trainInput is an array of structures and each
% entry contains three fields:
%       'input': matrix nInputUnits x nTrainingPoints. nTrainingPoints can
%       change from trial to trial.
%       'class': natural number that determines the class to which the
%       trial belongs to
%       'numClasses': total number of classes
% trainTeacher = teacher time series of size nOutputUnits x nTrainingPoints
% for forecasting tasks. For classification tasks it can be left empty
% since the class information is contained in structures in trainInput 
% esn = an ESN structure, through which we run our input sequence
% nForgetPoints - the first nForgetPoints will be disregarded
%
% optional arguments:
% 'ifRegressInput': 1 or 0 depending on the use or not of the input in the
%   computation of the output weights (default = 0)
%
% outputs:
% trained_esn = an Esn structure with the option trained = 1 and
% outputWeights and outputIntercept set.
% stateCollection = matrix of size nInternalUnits x (nTrainingPoints-nForgetPoints)
%
% Version 1.0, August 8, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega


trained_sas = sas;
trained_sas.ifRegressInput = 0;
args = varargin;
nargs= length(args);
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', trained_sas.ifRegressInput = args{i+1} ;
        otherwise error('the option does not exist');
    end
end


switch trained_sas.learningMode
    case 'offline'
        if ~isstruct(trainInput) && ~isstruct(trainTeacher)
            stateCollection = simulate_sas(trainInput, trainTeacher, trained_sas, nForgetPoints);
            if isempty(trainInput) && ~isempty(trainTeacher)
                stateCollection = stateCollection(:, 1:end-1);
            end
            [trained_sas.outputWeights, trained_sas.outputIntercept] = ...
                feval(trained_sas.methodWeightCompute, trainInput(:, nForgetPoints + 1:end), trainTeacher(:, nForgetPoints + 1:end),...
                stateCollection, trained_sas) ;
        else
            stateCollection = [];
            trainInputCollection = [];
            trainTeacherCollection = [];
            P = length(trainInput);
            numClasses = trainInput(1).numClasses;
            for i=1:P
                [~, lengthTrial] = size(trainInput(i).input);
                stateCollection = [stateCollection, simulate_sas(trainInput(i).input, [], trained_sas, nForgetPoints)];
                temptrainTeacherCollection = zeros(numClasses, lengthTrial - nForgetPoints);
                temptrainTeacherCollection(trainInput(i).class, :) = ones(1, lengthTrial - nForgetPoints);
                trainTeacherCollection = [trainTeacherCollection, temptrainTeacherCollection];
                trainInputCollection = [trainInputCollection, trainInputCollection(:, nForgetPoints + 1:end)];
            end
            [trained_sas.outputWeights, trained_sas.outputIntercept] = ...
                feval(trained_sas.methodWeightCompute, trainInput, trainTeacherCollection,...
                stateCollection, trained_sas) ;
        end
        
    case 'online'
        %     To be done
end

trained_sas.trained = 1 ;


