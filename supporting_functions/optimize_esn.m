function [esn_optim, NRMSE_optim] = optimize_esn(esn, trainInput, trainTeacher, nForgetPoints, varargin)

% Optimizes the specified parameters of an ESN with respect to the training
% error

%%%%% input arguments:
% esn: initial esn whose parameters have to be optimized. Does not need to
%       be trained
% trainInput = input time series of size nInputUnits x nTrainingPoints
% trainTeacher = teacher time series of size nOutputUnits x nTrainingPoints
% nForgetPoints = the first nForgetPoints will be disregarded

% optional arguments:
% 'ifRegressInput': 1 or 0 depending on the use or not of the input in the
%   computation of the output weights (default = 0)
%
% 'ifSaveOutput': 1 or 0 depending if we want to save the intermediate
% optimization steps in the file 'optimization_process.mat' (default = 0)
%
% 'ifMultiStart': integer determining the number of starting initial values
% in the solution of the optimization problem (default = 1). If bigger than
% one it cancels the ifSaveOutput option.
%
% The rest of arguments specify the hyperparameter that needs to be
% optimized. User needs to provide pairs with the name of the hyperparameter
% and the interval to which it belongs. The possibilities are:
%
% 'inputScaling': a scalar multiplyer in front of the inputMask.
%
% 'teacherScaling': a scalar multiplyer in front of the teacherMask.
%
% 'inputShift': a vector nInternalUnits x 1 with additive bias in the
%               reservoir equation
%
% 'spectralRadius': a positive number.
%
% 'ridgeLambda': small constant used for Ridge regularization

% Version 1.0, February 25, 2017
% Version 1.1, August 12, 2017: added multistart option
% Copyright: L. Grigoryeva and J.-P. Ortega

args = varargin;
nargs= length(args);
current = 1;
x0= [];
xlow = [];
xup = [];
esn.ifRegressInput = 0;
ifSaveOutput = 0;
ifMultiStart = 1;
psoptions = psoptimset('Display','iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 50000, 'UseParallel', 'Always', 'CompleteSearch', 'on', 'MaxFunEvals', 1000);
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', esn.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput',
            if args{i+1}
                ifSaveOutput = args{i+1} ;
                psoptions = psoptimset('OutputFcns', @outfun, 'Display','iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 50000, 'UseParallel', 'Always', 'CompleteSearch', 'on', 'MaxFunEvals', 1000);
                pathhistory = [];
                save('optimization_process.mat', 'pathhistory');
            end
        case 'ifMultiStart',
            if args{i+1} > 1
                ifMultiStart = args{i+1};
                ifSaveOutput = 0;
                options_multistart = optimoptions(@fmincon, 'Algorithm', 'sqp', 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'iter', 'UseParallel', 'Always');
            end
        case 'inputScaling',
            x0(current) = esn.inputScaling;
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        case 'teacherScaling',
            x0(current) = esn.teacherScaling;
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        case 'inputShift',
            x0(current) = esn.inputShift(1);
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        case 'spectralRadius',
            x0(current) = esn.spectralRadius;
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        case 'ridgeLambda',
            x0(current) = esn.ridgeLambda;
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        otherwise error('the option does not exist');
    end
end

funHandle = @(x)compute_train_error(x, esn, trainInput, trainTeacher, nForgetPoints, varargin);

if ifSaveOutput
    [optim_x, NRMSE_optim, ~] = patternsearch(funHandle, x0, [], [], [], [], xlow, xup, [], psoptions);
else
    problem = createOptimProblem('fmincon', 'objective', ...
        funHandle, 'x0', x0, 'lb', xlow, 'ub', xup, 'options', options_multistart);
    ms = MultiStart('UseParallel', true);
    [optim_x, NRMSE_optim, ~, ~,~] = run(ms, problem, ifMultiStart);
end

esn_optim = esn;
current = 1;
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', esn_optim.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput';
        case 'ifMultiStart';
        case 'inputScaling',
            esn_optim.inputScaling = optim_x(current);
            current = current + 1;
        case 'teacherScaling',
            esn_optim.teacherScaling = optim_x(current);
            current = current + 1;
        case 'inputShift',
            esn_optim.inputShift = optim_x(current) + zeros(esn_optim.nInternalUnits, 1);
            current = current + 1;
        case 'spectralRadius',
            esn_optim.spectralRadius = optim_x(current);
            esn_optim.reservoirMatrix = esn_optim.spectralRadius * esn_optim.internalWeights_UnitSR;
            current = current + 1;
        case 'ridgeLambda',
            esn_optim.ridgeLambda = optim_x(current);
            current = current + 1;
        otherwise error('the option does not exist');
    end
end
[esn_optim, ~] = ...
    train_esn(trainInput, trainTeacher, esn_optim, nForgetPoints, 'ifRegressInput', esn_optim.ifRegressInput);
end

function NRMSE = compute_train_error(x, esn, trainInput, trainTeacher, nForgetPoints, varargin)
esn.ifRegressInput = 0;
varargin = varargin{1};
args = varargin;
nargs= length(args);
current = 1;
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', esn.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput';
        case 'ifMultiStart';
        case 'inputScaling',
            esn.inputScaling = x(current);
            current = current + 1;
        case 'teacherScaling',
            esn.teacherScaling = x(current);
            current = current + 1;
        case 'inputShift',
            esn.inputShift = x(current) + zeros(esn.nInternalUnits, 1);
            current = current + 1;
        case 'spectralRadius',
            esn.spectralRadius = x(current);
            esn.reservoirMatrix = esn.spectralRadius * esn.internalWeights_UnitSR;
            current = current + 1;
        case 'ridgeLambda',
            esn.ridgeLambda = x(current);
            current = current + 1;
        otherwise error('the option does not exist');
    end
end

[trained_esn, stateCollection] = ...
    train_esn(trainInput, trainTeacher, esn, nForgetPoints, 'ifRegressInput', esn.ifRegressInput);
[stateCollectionTrain, predictedTrainOutput] = simulate_esn(trainInput, trainTeacher, trained_esn, nForgetPoints);
if isempty(trainInput)
    predictedTrainOutput = predictedTrainOutput(1: end-1);
end

NRMSE = compute_NRMSE(predictedTrainOutput, trainTeacher(:, nForgetPoints + 1:end));
end


function [stop,options,optchanged] = outfun(optimValues,options, flag)
stop = false;
switch flag
    case 'init'
        hold on
    case 'iter'
        % Concatenate current point and objective function
        load('optimization_process.mat')
        iteration = optimValues.iteration;
        pathhistory(iteration).fval = optimValues.fval;
        pathhistory(iteration).x = optimValues.x;
        save('optimization_process.mat', 'pathhistory');
        %                 pathhistory.gradient = [pathhistory.gradient optimValues.gradient];
    case 'done'
        hold off
    otherwise
end
options = options;
optchanged = false;
end



