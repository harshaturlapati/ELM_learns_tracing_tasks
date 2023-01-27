function [stateCollection, outputSequence] = simulate_sas(trainInput, trainTeacher, sas,...
    nForgetPoints, varargin)

% simulate_sas  runs the input through the SAS and writes the
% obtained reservoir states into stateCollection.
% The first nForgetPoints will be deleted, as the first few states could be
% not reliable due to initial transients
%
% input args:
% trainInput = input time series of size nInputUnits x nTrainingPoints
% trainTeacher = teacher time series of size nOutputUnits x nTrainingPoints
% sas = a SAS structure, through which we run our input sequence
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
% SAS is running autonomously; the SAS needs to be trained and in vargin,
% 'startingState' and 'startingTeacher' need to be declared
%
% optional input arguments:
% 'startingState': nInternalUnits x 1 vector.
% 'startingTeacher': nOutputUnits x 1 vector.
% 'lengthSimul': integer for the length of autonomous running
%
% output:
% stateCollectMat = matrix of size nInternalUnits x T. T explained below
% outputSequence = only for trained SASs. matrix of size nOutputUnits x T. T explained below
% The length T of the simulation depends on the way the ESN is simulated:
% If SAS is:
%     -input and teacher forced then T = (nTrainingPoints-nForgetPoints)
%     -input forced then T = (nTrainingPoints-nForgetPoints)
%     -teacher forced then T = (nTrainingPoints-nForgetPoints+1)
%     -autonomous running then T = (lengthSimul-nForgetPoints)
%
% Version 1.0, August 8, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega



% Initialize starting values
startingState = zeros(sas.nInternalUnits, 1);
startingTeacher = zeros(sas.nOutputUnits, 1);
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
    stateCollection = zeros(sas.nInternalUnits, nTrainingPoints);
    totalstate = sas.inputShift + Pimage(sas, trainInput(:, 1)) * startingState + Qimage(sas, trainInput(:, 1), sas.inputScaling);
    stateCollection(:, 1) = feval(sas.type, totalstate, sas);
    for i = 2:nTrainingPoints
        totalstate = sas.inputShift + Pimage(sas, trainInput(:, i)) * stateCollection(:, i - 1) + ...
            Qimage(sas, trainInput(:, i), sas.inputScaling);
        stateCollection(:, i) = feval(sas.type, totalstate, sas);
    end
elseif isempty(trainInput) && ~isempty(trainTeacher)
    [~, nTrainingPoints] = size(trainTeacher);
    nTrainingPoints = nTrainingPoints + 1;
    stateCollection = zeros(sas.nInternalUnits, nTrainingPoints);
    totalstate = sas.inputShift + Rimage(sas, startingTeacher) * startingState + Simage(sas, startingTeacher, sas.teacherScaling);
    stateCollection(:, 1) = feval(sas.type, totalstate, sas);
    for i = 2:nTrainingPoints
        totalstate = sas.inputShift + Rimage(sas, trainTeacher(:, i - 1)) * stateCollection(:, i - 1) + ...
            Simage(sas, trainTeacher(:, i - 1), sas.teacherScaling);
        stateCollection(:, i) = feval(sas.type, totalstate, sas);
    end
elseif ~isempty(trainInput) && ~isempty(trainTeacher)
    [~, nTrainingPoints] = size(trainInput);
    stateCollection = zeros(sas.nInternalUnits, nTrainingPoints);
    totalstate = sas.inputShift + Pimage(sas, trainInput(:, 1)) * startingState + Qimage(sas, trainInput(:, 1), sas.inputScaling) +...
        Rimage(sas, startingTeacher) * startingState + Simage(sas, startingTeacher, sas.teacherScaling);
    stateCollection(:, 1) = feval(sas.type, totalstate, sas);
    for i = 2:nTrainingPoints
        totalstate = sas.inputShift + Pimage(sas, trainInput(:, i)) * stateCollection(:, i - 1) + ...
            Qimage(sas, trainInput(:, i), sas.inputScaling) + ...
            Rimage(sas, trainTeacher(:, i - 1)) * stateCollection(:, i - 1) + Simage(sas, trainTeacher(:, i - 1), sas.teacherScaling);
        stateCollection(:, i) = feval(sas.type, totalstate, sas);
    end
else
    if sas.trained == 0
        error('the ESN is not trained and cannot run autonomously');
    else
        stateCollection = zeros(sas.nInternalUnits, lengthSimul);
        totalstate = sas.inputShift + Rimage(sas, startingTeacher) * startingState + Simage(sas, startingTeacher, sas.teacherScaling);
        stateCollection(:, 1) = feval(sas.type, totalstate, sas);
        teacher = feval(sas.outputActivationFunction, sas.outputIntercept + sas.outputWeights' * stateCollection(:, 1));
        for i = 2:lengthSimul
            totalstate = sas.inputShift + Rimage(sas, teacher) * stateCollection(:, i - 1) + Simage(sas, teacher, sas.teacherScaling);
            stateCollection(:, i) = feval(sas.type, totalstate, sas);
            teacher = feval(sas.outputActivationFunction, sas.outputIntercept + sas.outputWeights' * stateCollection(:, i));
        end
    end
end

if sas.trained == 0
    outputSequence = [];
else
    if sas.ifRegressInput == 0
        x = stateCollection;
    else
        x =[stateCollection; trainInput];
    end
    [~, T] = size(x);
    outputSequence = zeros(sas.nOutputUnits, T);
    for i=1:T
        outputSequence(:, i) = feval(sas.outputActivationFunction, sas.outputIntercept + sas.outputWeights' * x(:, i));
    end
end

stateCollection = stateCollection(:, nForgetPoints + 1:end);
outputSequence = outputSequence(:, nForgetPoints + 1:end);
end
function Pz = Pimage(sas, input)
order = sas.orderP;
if ~isempty(order)
    temp = sas.Pmatrices(:, :, 1);
    for i = 1:order
        temp = temp + sas.Pmatrices(:, :, i + 1) * input^i;
    end
    Pz = temp;
else
    Pz = 0;
end
end

function Qz = Qimage(sas, input, inputScaling)
order = sas.orderQ;
if ~isempty(order)
    if length(inputScaling) == 1
        temp = 0;
        for i = 1:order
            temp = temp + sas.Qmatrices(:, i) * input^i;
        end
        Qz = inputScaling * temp;
    else
        temp = 0;
        for i = 1:order
            temp = temp + inputScaling(i) * sas.Qmatrices(:, i) * input^i;
        end
        Qz = temp;
    end
else
    Qz = 0;
end
end

function Rd = Rimage(sas, teacher)
order = sas.orderR;
if ~isempty(order)
    temp = sas.Rmatrices(:, :, 1);
    for i = 1:order
        temp = temp + sas.Rmatrices(:, :, i + 1) * teacher^i;
    end
    Rd = temp;
else
    Rd = 0;
end
end

function Sd = Simage(sas, teacher, teacherScaling)
order = sas.orderS;
if ~isempty(order)
    if length(teacherScaling) == 1
        temp = 0;
        for i = 1:order
            temp = temp + sas.Smatrices(:, i) * teacher^i;
        end
        Sd = teacherScaling * temp;
    else
        temp = 0;
        for i = 1:order
            temp = temp + teacherScaling(i) * sas.Smatrices(:, i) * teacher^i;
        end
        Sd = temp;
    end
else
    Sd = 0;
end
end

