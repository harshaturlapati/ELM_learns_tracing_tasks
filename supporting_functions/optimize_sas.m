function [sas_optim, NRMSE_optim] = optimize_sas(sas, trainInput, trainTeacher, nForgetPoints, varargin)

% Optimizes the specified hyperparameters of a SAS with respect to the training
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
% 'spectralRadiusP': vector of size (orderP + 1) x 1
%
% 'spectralRadiusR': vector of size (orderR + 1) x 1
%
% 'ridgeLambda': small constant used for Ridge regularization

% Version 1.0, August 13, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega

args = varargin;
nargs= length(args);
current = 1;
x0= [];
xlow = [];
xup = [];
sas.ifRegressInput = 0;
ifSaveOutput = 0;
ifMultiStart = 1;
psoptions = psoptimset('Display','iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 50000, 'UseParallel', 'Always', 'CompleteSearch', 'on', 'MaxFunEvals', 1000);
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', sas.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput',
            if args{i+1}
                ifSaveOutput = args{i+1} ;
                psoptions = psoptimset('OutputFcns', @outfun, 'Display','iter', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 50000, 'UseParallel', 'Always', 'CompleteSearch', 'on', 'MaxFunEvals', 5000);
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
            linput = length(sas.inputScaling);
            x0(current:current + linput - 1) = sas.inputScaling;
            xlow = [xlow args{i+1}(1:linput)];
            xup = [xup args{i+1}(linput + 1:end)];
            current = current + linput;
        case 'teacherScaling',
            lteacher = length(sas.teacherScaling);
            x0(current:current + lteacher - 1) = sas.teacherScaling;
            xlow = [xlow args{i+1}(1:lteacher)];
            xup = [xup args{i+1}(lteacher + 1:end)];
            current = current + lteacher;
        case 'inputShift',
            x0(current) = sas.inputShift(1);
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        case 'spectralRadiusP',
            x0(current:current + sas.orderP) = sas.spectralRadiusP';
            xlow = [xlow args{i+1}(1:sas.orderP + 1)];
            xup = [xup args{i+1}(sas.orderP + 2:end)];
            current = current + sas.orderP + 1;
        case 'spectralRadiusR',
            x0(current:current + sas.orderR) = sas.spectralRadiusR';
            xlow = [xlow args{i+1}(1:sas.orderR + 1)];
            xup = [xup args{i+1}(sas.orderR + 2:end)];
            current = current + sas.orderR + 1;
        case 'ridgeLambda',
            x0(current) = sas.ridgeLambda;
            xlow = [xlow args{i+1}(1)];
            xup = [xup args{i+1}(2)];
            current = current + 1;
        otherwise error('the option does not exist');
    end
end

funHandle = @(x)compute_train_error(x, sas, trainInput, trainTeacher, nForgetPoints, varargin);

if ifSaveOutput
    [optim_x, NRMSE_optim, ~] = patternsearch(funHandle, x0, [], [], [], [], xlow, xup, [], psoptions);
else
    problem = createOptimProblem('fmincon', 'objective', ...
        funHandle, 'x0', x0, 'lb', xlow, 'ub', xup, 'options', options_multistart);
    ms = MultiStart('UseParallel', true);
    [optim_x, NRMSE_optim, ~, ~,~] = run(ms, problem, ifMultiStart);
end

sas_optim = sas;
current = 1;
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', sas_optim.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput';
        case 'ifMultiStart';
        case 'inputScaling',
            sas_optim.inputScaling = optim_x(current:current + length(sas.inputScaling) - 1);
            current = current + length(sas.inputScaling);
        case 'teacherScaling',
            sas_optim.teacherScaling = optim_x(current:current + length(sas.teacherScaling) - 1);
            current = current + length(sas.teacherScaling);
        case 'inputShift',
            sas_optim.inputShift = optim_x(current) + zeros(sas_optim.nInternalUnits, 1);
            current = current + 1;
        case 'spectralRadiusP',
            sas_optim.spectralRadiusP = optim_x(current:current + sas.orderP)';
                for i = 1:sas.orderP + 1
                    sas_optim.Pmatrices(:, :, i) = sas_optim.Pmatrices_UnitSR(:, :, i) * sas_optim.spectralRadiusP(i);
                end
            current = current + sas_optim.orderP + 1;
        case 'spectralRadiusR',
            sas_optim.spectralRadiusR = optim_x(current:current + sas.orderR)';
                for i = 1:sas.orderR + 1
                    sas_optim.Rmatrices(:, :, i) = sas_optim.Rmatrices_UnitSR(:, :, i) * sas_optim.spectralRadiusR(i);
                end
            current = current + sas_optim.orderR + 1;
        case 'ridgeLambda',
            sas_optim.ridgeLambda = optim_x(current);
            current = current + 1;
        otherwise error('the option does not exist');
    end
end
[sas_optim, ~] = ...
    train_sas(trainInput, trainTeacher, sas_optim, nForgetPoints, 'ifRegressInput', sas_optim.ifRegressInput);
end

function NRMSE = compute_train_error(x, sas, trainInput, trainTeacher, nForgetPoints, varargin)
sas.ifRegressInput = 0;
varargin = varargin{1};
args = varargin;
nargs= length(args);
current = 1;
for i=1:2:nargs
    switch args{i},
        case 'ifRegressInput', sas.ifRegressInput = args{i+1} ;
        case 'ifSaveOutput';
        case 'ifMultiStart';
        case 'inputScaling',
            sas.inputScaling = x(current:current + length(sas.inputScaling) - 1);
            current = current + length(sas.inputScaling);
        case 'teacherScaling',
            sas.teacherScaling = x(current:current + length(sas.teacherScaling) - 1);
            current = current + length(sas.teacherScaling);
        case 'inputShift',
            sas.inputShift = x(current) + zeros(sas.nInternalUnits, 1);
            current = current + 1;
        case 'spectralRadiusP',
            sas.spectralRadiusP = x(current:current + sas.orderP)';
                for i = 1:sas.orderP + 1
                    sas.Pmatrices(:, :, i) = sas.Pmatrices_UnitSR(:, :, i) * sas.spectralRadiusP(i);
                end
            current = current + sas.orderP + 1;
        case 'spectralRadiusR',
            sas.spectralRadiusR = x(current:current + sas.orderR)';
                for i = 1:sas.orderR + 1
                    sas.Rmatrices(:, :, i) = sas.Rmatrices_UnitSR(:, :, i) * sas.spectralRadiusR(i);
                end
            current = current + sas.orderR + 1;
        case 'ridgeLambda',
            sas.ridgeLambda = x(current);
            current = current + 1;
        otherwise error('the option does not exist');
    end
end

[trained_sas, stateCollection] = ...
    train_sas(trainInput, trainTeacher, sas, nForgetPoints, 'ifRegressInput', sas.ifRegressInput);
[stateCollectionTrain, predictedTrainOutput] = simulate_sas(trainInput, trainTeacher, trained_sas, nForgetPoints);
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



