function [trained_esn, stateCollection] = ...
    train_esn(trainInput, trainTeacher , esn, nForgetPoints, varargin)
% TRAIN_ESN Trains the output weights of an ESN for forecasting and
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
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 30, 2006, H. Jaeger
% Revision 2, Feb 23, 2007, H. Jaeger
% Revision 3, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: adapted to
% different conventions. Only online and offline options
% Revision 4, May 3, 2017, L. Grigoryeva and J.-P. Ortega: adapted to
% classification problems


trained_esn = esn;
trained_esn.ifRegressInput = 0;
args = varargin;
nargs= length(args);
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', trained_esn.ifRegressInput = args{i+1} ;
        otherwise error('the option does not exist');
    end
end


switch trained_esn.learningMode
    case 'offline'
        if ~isstruct(trainInput) && ~isstruct(trainTeacher)
            stateCollection = simulate_esn(trainInput, trainTeacher, trained_esn, nForgetPoints);
            if isempty(trainInput) && ~isempty(trainTeacher)
                stateCollection = stateCollection(:, 1:end-1);
            end
            [trained_esn.outputWeights, trained_esn.outputIntercept] = ...
                feval(trained_esn.methodWeightCompute, trainInput(:, nForgetPoints + 1:end), trainTeacher(:, nForgetPoints + 1:end),...
                stateCollection, trained_esn) ;
        else
            stateCollection = [];
            trainInputCollection = [];
            trainTeacherCollection = [];
            P = length(trainInput);
            numClasses = trainInput(1).numClasses;
            for i=1:P
                [~, lengthTrial] = size(trainInput(i).input);
                stateCollection = [stateCollection, simulate_esn(trainInput(i).input, [], trained_esn, nForgetPoints)];
                temptrainTeacherCollection = zeros(numClasses, lengthTrial - nForgetPoints);
                temptrainTeacherCollection(trainInput(i).class, :) = ones(1, lengthTrial - nForgetPoints);
                trainTeacherCollection = [trainTeacherCollection, temptrainTeacherCollection];
                trainInputCollection = [trainInputCollection, trainInputCollection(:, nForgetPoints + 1:end)];
            end
            [trained_esn.outputWeights, trained_esn.outputIntercept] = ...
                feval(trained_esn.methodWeightCompute, trainInput, trainTeacherCollection,...
                stateCollection, trained_esn) ;
        end
        
    case 'online'
        %     To be done
end

trained_esn.trained = 1 ;


