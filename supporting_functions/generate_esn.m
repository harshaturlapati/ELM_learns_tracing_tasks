function esn = generate_esn(nInputUnits, nInternalUnits, nOutputUnits, varargin)
% Creates an ESN set up for use in multiple-channel output association tasks.
% The number of input, internal, and output
% units have to be set. Any other option is set using the format
% 'name_of_options1',value1,'name_of_option2',value2, etc.
%
%%%%% input arguments:
% nInputUnits: the dimension of the input
% nInternalUnits: size of the Esn
% nOutputUnits: the dimension of the output
%
%%%%% optional arguments:
% 'inputMask': a nInternalUnits x nInputUnits matrix. Mask in front of z.
%
% 'inputScaling': a scalar multiplyer in front of the inputMask.
%
% 'teacherMask': a nInternalUnits x nOutputUnits matrix. Mask in front of d.
%
% 'teacherScaling': a scalar multiplyer in front of the teacherMask.
%
% 'inputShift': a vector nInternalUnits x 1 with additive bias in the
%               reservoir equation
%
% 'connectivity': sparsity degree between 0 and 1 of the reservoir matrix.
%                 0 totally sparse, 1 full matrix
%
% 'spectralRadius': a positive number less than 1.
%
% 'reservoirMatrix': nInternalUnits x nInternalUnits matrix with the
%                    internal weights
% 'ridgeLambda': small constant used for Ridge regularization
%
% 'learningMode': a string ('offline' or 'online')
%     1. Case 'offline': trainInput and trainOutput each represent a
%        time series in an array of size  sequenceDimension x sequenceLength
%     2. Case 'online': trainInput and trainOutput are a single time
%        series, output weights are adapted online
%
% 'reservoirActivationFunction': a string ("tanh", "identity", "sigmoid01") ,
%
% 'outputActivationFunction': a string("tanh", "identity", "sigmoid01") ,
%
% 'inverseOutputActivationFunction': the inverse to
%    outputActivationFunction, one of 'atanh', 'identity', 'sigmoid01_inv'.
%    When choosing the activation function, make sure the inverse
%    activation function is corectly set.
%
% 'methodWeightCompute': a string ('pseudoinverse' or 'gaussian_elimination'). It
%    specifies which method to use to compute the output weights given the
%    state collection matrix and the teacher
%
% 'type': a string ('plain_esn', 'leaky_esn' or 'twi_esn')
% 'trained': a flag indicating whether the network has been trained already
% 'timeConstants': option used in networks with type == "leaky_esn", "leaky1_esn" and "twi_esn".
%                      Is given as column vector of size esn.nInternalUnitsm, where each entry
%                      signifies a time constant for a reservoir neuron.
% 'leakage': option used in networks with type == "leaky_esn" or "twi_esn"
% 'RLS_lambda': option used in online training(learningMode == "online")
% 'RLS_delta': option used in online training(learningMode == "online")
%
% for more information on the Echo State network approach take a look at
% the following tutorial :
% http://www.faculty.iu-bremen.de/hjaeger/pubs/ESNTutorialRev.pdf

% Version 1.0, April 30, 2006
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 6, 2006, H. Jaeger
% Revision 2, Feb 23, 2007, H. Jaeger
% Revision 3, June 27, 2007, H. Jaeger
% Revision 4, July 1, 2007, H. Jaeger:
%    - changed esn.timeConstant to esn.timeConstants
%    - deleted esn.retainment
%    - deleted esn.internalWeights (to enforce that this is set outside
%    this generation script)
% Revision 5, July 29, 2007, H. Jaeger: bugfix (for cases of zero input or
%    output length)
% Revision 6, Jan 28, 2009, H. Jaeger: bugfix (deleted defunct error
%                             catching routine for leaky / twiesn
%                             reservoirs)
% Revision 7, Feb 14, 2017, L. Grigoryeva and J.-P. Ortega: adapted to new
% conventions and new options

%%%% set the number of units
esn.nInternalUnits = nInternalUnits;
esn.nInputUnits = nInputUnits;
esn.nOutputUnits = nOutputUnits;

%init default parameters
esn.inputMask = 2.0 * rand(nInternalUnits, nInputUnits) - 1.0;
esn.inputScaling = 1;
esn.teacherMask = 2.0 * rand(nInternalUnits, nOutputUnits) - 1.0;
esn.teacherScaling = 1;
esn.inputShift = zeros(nInternalUnits, 1);

esn.connectivity = min([10/nInternalUnits 1]);
esn.spectralRadius = 1 ;
esn.internalWeights_UnitSR = generate_internal_weights(nInternalUnits, ...
    esn.connectivity);
esn.reservoirMatrix = esn.spectralRadius * esn.internalWeights_UnitSR;
esn.ridgeLambda = 0;

esn.learningMode = 'offline_singleTimeSeries' ;
esn.reservoirActivationFunction = 'tanh';
esn.outputActivationFunction = 'identity' ;
esn.inverseOutputActivationFunction = 'identity' ;
esn.methodWeightCompute = 'gaussian_elimination' ;

esn.type = 'plain_esn' ;
esn.trained = 0 ;


args = varargin;
nargs= length(args);
for i=1:2:nargs
    switch args{i},
        case 'inputMask', esn.inputMask = args{i+1};
        case 'inputScaling', esn.inputScaling = args{i+1};
        case 'teacherMask', esn.teacherMask = args{i+1};
        case 'teacherScaling', esn.teacherScaling = args{i+1};
        case 'inputShift', esn.inputShift = args{i+1};
        case 'connectivity', esn.connectivity = args{i+1} ;
            esn.internalWeights_UnitSR = generate_internal_weights(nInternalUnits, ...
                esn.connectivity);
            esn.reservoirMatrix = esn.spectralRadius * esn.internalWeights_UnitSR;
        case 'spectralRadius', esn.spectralRadius = args{i+1} ;
            esn.reservoirMatrix = esn.spectralRadius * esn.internalWeights_UnitSR;
        case 'reservoirMatrix', esn.reservoirMatrix = args{i+1} ;
        case 'ridgeLambda', esn.ridgeLambda = args{i+1} ;
        case 'learningMode', esn.learningMode = args{i+1} ;
        case 'reservoirActivationFunction',esn.reservoirActivationFunction=args{i+1};
        case 'outputActivationFunction',esn.outputActivationFunction=  ...
                args{i+1};
        case 'inverseOutputActivationFunction', esn.inverseOutputActivationFunction=args{i+1};
        case 'methodWeightCompute', esn.methodWeightCompute = args{i+1} ;
        case 'type' , esn.type = args{i+1} ;         
        otherwise error('the option does not exist');
    end
end

%%%% error checking
% check that inputScaling has correct format
if size(esn.inputMask,1) ~= nInternalUnits
    error('the size of the inputMask does not match the size of the reservoir');
end
if size(esn.inputMask,2) ~= nInputUnits
    error('the size of the inputMask does not dimensionality of the input signal');
end
if size(esn.teacherMask,1) ~= nInternalUnits
    error('the size of the teacherMask does not match the size of the reservoir');
end
if size(esn.teacherMask,2) ~= nOutputUnits
    error('the size of the teacherMask does not dimensionality of the teacher signal');
end
if esn.connectivity > 1  || esn.connectivity < 0
    error('the connectivity has to be between 0 and 1');
end
if ~strcmp(esn.learningMode,'offline') &&...
        ~strcmp(esn.learningMode,'online')
    error('learningMode should be either "offline" or "online" ') ;
end
if ~((strcmp(esn.outputActivationFunction,'identity') && ...
        strcmp(esn.inverseOutputActivationFunction,'identity')) || ...
        (strcmp(esn.outputActivationFunction,'tanh') && ...
        strcmp(esn.inverseOutputActivationFunction,'atanh')) || ...
        (strcmp(esn.outputActivationFunction,'sigmoid01') && ...
        strcmp(esn.inverseOutputActivationFunction,'sigmoid01_inv')))  ...
        error('outputActivationFunction and inverseOutputActivationFunction do not match');
end


