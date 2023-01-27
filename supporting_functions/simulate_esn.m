function [stateCollection, outputSequence] = simulate_esn(trainInput, trainTeacher, esn,...
    nForgetPoints, varargin)

% simulate_esn_GO  runs the input through the ESN and writes the
% obtained reservoir states into stateCollection.
% The first nForgetPoints will be deleted, as the first few states could be
% not reliable due to initial transients
%
% input args:
% trainInput = input time series of size nInputUnits x nTrainingPoints
% trainTeacher = teacher time series of size nOutputUnits x nTrainingPoints
% esn = an ESN structure, through which we run our input sequence
% nForgetPoints: an integer, may be negative, positive or zero.
%    If positive: the first nForgetPoints will be disregarded (washing out
%    initial reservoir transient)
%    If negative: the network will be initially driven from zero state with
%    the first input repeated |nForgetPoints| times; size(inputSequence,1)
%    many states will be sorted into state matrix
%    If zero: no washout accounted for, all states except the zero starting
%    state will be sorted into state matrix
%
% Note: trainInput and/or trainTeacher may be the empty list [].
% If the inputSequence is empty, we are dealing with a purely
% generative task; states are then computed by teacher-forcing
% trainTeacher. If both trainInput and trainTeacher are [] then the
% ESN is running autonomously; the ESN needs to be trained and in vargin initial state,
% and initial teacher need to be declared
%
% optional input arguments:
% 'startingState': nInternalUnits x 1 vector.
% 'startingTeacher': nOutputUnits x 1 vector.
% 'lengthSimul': integer for the length of autonomous running
%
% output:
% stateCollectMat = matrix of size nInternalUnits x T. T explained below
% outputSequence = only for trained ESNs. matrix of size nOutputUnits x T. T explained below
% The length T of the simulation depends on the way the ESN is simulated:
% If ESN is: 
%     -input and teacher forced then T = (nTrainingPoints-nForgetPoints)
%     -input forced then T = (nTrainingPoints-nForgetPoints)
%     -teacher forced then T = (nTrainingPoints-nForgetPoints+1)
%     -autonomous running then T = (lengthSimul-nForgetPoints)
%
% Version 1.0, February 13, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega



% Initialize starting values
startingState = zeros(esn.nInternalUnits, 1);
startingTeacher = zeros(esn.nOutputUnits, 1);
lengthSimul = 0;
args = varargin;
nargs= length(args);
for i=1:2:nargs
    switch args{i}
        case 'startingState', startingState = args{i+1} ;
        case 'startingTeacher', startingTeacher = args{i+1} ;
        case 'lengthSimul', lengthSimul = args{i+1} ;
        otherwise error('the option does not exist');
    end
end




if isempty(trainTeacher) && ~isempty(trainInput)
    [~, nTrainingPoints] = size(trainInput);
    stateCollection = zeros(esn.nInternalUnits, nTrainingPoints);
        totalstate = esn.inputShift + esn.reservoirMatrix * startingState + esn.inputScaling * esn.inputMask * trainInput(:, 1);
    stateCollection(:, 1) = feval(esn.type, totalstate, esn);
    for i = 2:nTrainingPoints
        totalstate = esn.inputShift + esn.reservoirMatrix * stateCollection(:, i - 1) + ...
            esn.inputScaling * esn.inputMask * trainInput(:, i);
        stateCollection(:, i) = feval(esn.type, totalstate, esn);
    end
elseif isempty(trainInput) && ~isempty(trainTeacher)
    [~, nTrainingPoints] = size(trainTeacher);
    nTrainingPoints = nTrainingPoints + 1;
    stateCollection = zeros(esn.nInternalUnits, nTrainingPoints);
    totalstate = esn.inputShift + esn.reservoirMatrix * startingState + esn.teacherScaling * esn.teacherMask * startingTeacher;
    stateCollection(:, 1) = feval(esn.type, totalstate, esn);
    for i = 2:nTrainingPoints
        totalstate = esn.inputShift + esn.reservoirMatrix * stateCollection(:, i - 1) + ...
            esn.teacherScaling * esn.teacherMask * trainTeacher(:, i - 1);
        stateCollection(:, i) = feval(esn.type, totalstate, esn);
    end
elseif ~isempty(trainInput) && ~isempty(trainTeacher)
    [~, nTrainingPoints] = size(trainInput);

    stateCollection = zeros(esn.nInternalUnits, nTrainingPoints);
    totalstate = esn.inputShift + esn.reservoirMatrix * startingState + esn.inputScaling * esn.inputMask * trainInput(:, 1) +...
        esn.teacherScaling * esn.teacherMask * startingTeacher;
    stateCollection(:, 1) = feval(esn.type, totalstate, esn);
    for i = 2:nTrainingPoints
        totalstate = esn.inputShift + esn.reservoirMatrix * stateCollection(:, i - 1) + ...
            esn.inputScaling * esn.inputMask * trainInput(:, i) + esn.teacherScaling * esn.teacherMask * trainTeacher(:, i - 1);
        stateCollection(:, i) = feval(esn.type, totalstate, esn);
    end
else
    if esn.trained == 0
        error('the ESN is not trained and cannot run autonomously');
    else
        stateCollection = zeros(esn.nInternalUnits, lengthSimul);
        totalstate = esn.inputShift + esn.reservoirMatrix * startingState + esn.teacherScaling * esn.teacherMask * startingTeacher;
        stateCollection(:, 1) = feval(esn.type, totalstate, esn);
        teacher = feval(esn.outputActivationFunction, esn.outputIntercept + esn.outputWeights' * stateCollection(:, 1));
        for i = 2:lengthSimul
            totalstate = esn.inputShift + esn.reservoirMatrix * stateCollection(:, i - 1) + esn.teacherScaling * esn.teacherMask * teacher;
            stateCollection(:, i) = feval(esn.type, totalstate, esn);
            teacher = feval(esn.outputActivationFunction, esn.outputIntercept + esn.outputWeights' * stateCollection(:, i));
        end
    end
end

if esn.trained == 0
    outputSequence = [];
else
    if esn.ifRegressInput == 0
        x = stateCollection;
    else
        x =[stateCollection; trainInput];
    end
    [~, T] = size(x);
    outputSequence = zeros(esn.nOutputUnits, T);
    for i=1:T
        outputSequence(:, i) = feval(esn.outputActivationFunction, esn.outputIntercept + esn.outputWeights' * x(:, i));
    end
end

stateCollection = stateCollection(:, nForgetPoints + 1:end);
outputSequence = outputSequence(:, nForgetPoints + 1:end);

