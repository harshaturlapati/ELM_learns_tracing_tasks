function sas = generate_sas(nInputUnits, nInternalUnits, nOutputUnits, orderP, orderQ, orderR, orderS, varargin)
% Creates a non-homogeneous state-affine system SAS for use in multiple-channel
% output association tasks.
% The number of input, internal, output, and the order of the P and Q
% polynomials have to be set.
% Any other option is set using the format
% 'name_of_options1',value1,'name_of_option2',value2, etc.
%
%%%%% input arguments:
% nInputUnits: the dimension of the input
% nInternalUnits: size of the SAS
% nOutputUnits: the dimension of the output
% orderP: order of P polynomial. P is set to 0 if orderP = []
% orderQ: order of Q polynomial. Q is set to 0 if orderQ = []
% orderR: order of R polynomial. R is set to 0 if orderR = []
% orderS: order of S polynomial. S is set to 0 if orderS = []
%
%%%%% optional arguments:
% 'Pmatrices': a nInternalUnits x nInternalUnits x (orderP + 1) array.
%
% 'Qmatrices': a nInternalUnits x orderQ array.
%
% 'Rmatrices': a nInternalUnits x nInternalUnits x (orderR + 1) array.
%
% 'Smatrices': a nInternalUnits x orderS array.
%
% 'inputScaling': it can be either a scalar multiplyer in front of the Q polynomial 
%                 or an orderQ x 1 vector that contains a multiplier for each monomial. Default = 1
%
% 'teacherScaling': a scalar multiplyer in front of the S polynomial
%                   or an orderS x 1 vector that contains a multiplier for each monomial. Default = 1
%
% 'inputShift': a vector nInternalUnits x 1 with additive bias in the
%               reservoir equation. Default is zero
%
% 'spectralRadiusP': vector of size (orderP + 1) x 1 that determines the
%                   spectral radii of the matrices in the P polynomial.
%
% 'spectralRadiusR': vector of size (orderR + 1) x 1 that determines the
%                   spectral radii of the matrices in the R polynomial.
%
% 'connectivity': sparsity degree between 0 and 1 of the matrices in P and R polynomials.
%                 0 totally sparse, 1 full matrix. If this parameter is
%                 specified together with 'Pmatrices' or 'Rmatrices' then they are
%                 newly random generated
%
% 'ridgeLambda': small constant used for Ridge regularization. Default is 0
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
% 'type': a string ('plain_sas', 'leaky_sas' or 'twi_sas')
% 'trained': a flag indicating whether the SAS has been trained already
% 'timeConstants': option used in networks with type == "leaky_esn", "leaky1_esn" and "twi_esn".
%                      Is given as column vector of size esn.nInternalUnitsm, where each entry
%                      signifies a time constant for a reservoir neuron.
%
% Version 1.0, August 7, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega

%%%% set the number of units
sas.nInternalUnits = nInternalUnits;
sas.nInputUnits = nInputUnits;
sas.nOutputUnits = nOutputUnits;
sas.orderP = orderP;
sas.orderQ = orderQ;
sas.orderR = orderR;
sas.orderS = orderS;

%init default parameters
sas.inputScaling = 1;
sas.teacherScaling = 1;
sas.inputShift = zeros(nInternalUnits, 1);

if ~isempty(orderP)
    sas.spectralRadiusP = (1/(orderP + 1)) * ones(orderP + 1, 1);
else
    sas.spectralRadiusP = [];
end
if ~isempty(orderR)
    sas.spectralRadiusR = (1/(orderR + 1)) * ones(orderR + 1, 1);
else
    sas.spectralRadiusR = [];
end
sas.connectivity = min([10/nInternalUnits 1]);

sas.Pmatrices_UnitSR = zeros(nInternalUnits, nInternalUnits, orderP + 1);
sas.Pmatrices = zeros(nInternalUnits, nInternalUnits, orderP + 1);
if ~isempty(orderP)
    for i = 1:orderP + 1
        sas.Pmatrices_UnitSR(:, :, i) = generate_internal_weights(nInternalUnits, sas.connectivity);
        sas.Pmatrices(:, :, i) = sas.Pmatrices_UnitSR(:, :, i) * sas.spectralRadiusP(i);
    end
    sas.Pmatrices_UnitSR = ndSparse(sas.Pmatrices_UnitSR);
    sas.Pmatrices = ndSparse(sas.Pmatrices);
end
if ~isempty(orderQ)
    sas.Qmatrices = 2.0 * rand(nInternalUnits, orderQ) - 1.0;
else
    sas.Qmatrices = zeros(nInternalUnits, orderQ);
end
sas.Rmatrices_UnitSR = zeros(nInternalUnits, nInternalUnits, orderR + 1);
sas.Rmatrices = zeros(nInternalUnits, nInternalUnits, orderR + 1);
if ~isempty(orderR)
    for i = 1:orderR + 1
        sas.Rmatrices_UnitSR(:, :, i) = generate_internal_weights(nInternalUnits, sas.connectivity);
        sas.Rmatrices(:, :, i) = sas.Rmatrices_UnitSR(:, :, i) * sas.spectralRadiusR(i);
    end
    sas.Rmatrices_UnitSR = ndSparse(sas.Rmatrices_UnitSR);
    sas.Rmatrices = ndSparse(sas.Rmatrices);
end
if ~isempty(orderS)
    sas.Smatrices = 2.0 * rand(nInternalUnits, orderS) - 1.0;
else
    sas.Smatrices = zeros(nInternalUnits, orderS);
end
sas.ridgeLambda = 0;

sas.learningMode = 'offline_singleTimeSeries' ;
sas.reservoirActivationFunction = 'identity';
sas.outputActivationFunction = 'identity' ;
sas.inverseOutputActivationFunction = 'identity' ;
sas.methodWeightCompute = 'gaussian_elimination' ;

sas.type = 'plain_esn' ;
sas.trained = 0 ;


args = varargin;
nargs= length(args);
for i=1:2:nargs
    switch args{i},
        case 'Pmatrices', sas.Pmatrices = args{i+1};
        case 'Qmatrices', sas.Qmatrices = args{i+1};
        case 'Rmatrices', sas.Rmatrices = args{i+1};
        case 'Smatrices', sas.Smatrices = args{i+1};
        case 'inputScaling', sas.inputScaling = args{i+1};
        case 'teacherScaling', sas.teacherScaling = args{i+1};
        case 'inputShift', sas.inputShift = args{i+1};
        case 'spectralRadiusP', sas.spectralRadiusP = args{i+1} ;
        case 'spectralRadiusR', sas.spectralRadiusR = args{i+1} ;
        case 'connectivity', sas.connectivity = args{i+1} ;
            sas.Pmatrices_UnitSR = zeros(nInternalUnits, nInternalUnits, orderP + 1);
            sas.Pmatrices = zeros(nInternalUnits, nInternalUnits, orderP + 1);
            if ~isempty(orderP)
                for i = 1:orderP + 1
                    sas.Pmatrices_UnitSR(:, :, i) = generate_internal_weights(nInternalUnits, sas.connectivity);
                    sas.Pmatrices(:, :, i) = sas.Pmatrices_UnitSR(:, :, i) * sas.spectralRadiusP(i);
                end
                sas.Pmatrices_UnitSR = ndSparse(sas.Pmatrices_UnitSR);
                sas.Pmatrices = ndSparse(sas.Pmatrices);
            end
            sas.Rmatrices_UnitSR = zeros(nInternalUnits, nInternalUnits, orderR + 1);
            sas.Rmatrices = zeros(nInternalUnits, nInternalUnits, orderR + 1);
            if ~isempty(orderR)
                for i = 1:orderR + 1
                    sas.Rmatrices_UnitSR(:, :, i) = generate_internal_weights(nInternalUnits, sas.connectivity);
                    sas.Rmatrices(:, :, i) = sas.Rmatrices_UnitSR(:, :, i) * sas.spectralRadiusR(i);
                end
                sas.Rmatrices_UnitSR = ndSparse(sas.Rmatrices_UnitSR);
                sas.Rmatrices = ndSparse(sas.Rmatrices);
            end
        case 'ridgeLambda', sas.ridgeLambda = args{i+1} ;
        case 'learningMode', sas.learningMode = args{i+1} ;
        case 'reservoirActivationFunction',sas.reservoirActivationFunction=args{i+1};
        case 'outputActivationFunction',sas.outputActivationFunction=  ...
                args{i+1};
        case 'inverseOutputActivationFunction', sas.inverseOutputActivationFunction=args{i+1};
        case 'methodWeightCompute', sas.methodWeightCompute = args{i+1} ;
        case 'type' , sas.type = args{i+1} ;
        otherwise error('the option does not exist');
    end
end

%%%% error checking
% check that inputScaling has correct format
if size(sas.Pmatrices,1) ~= nInternalUnits
    error('the size of the Pmatrix does not match the size of the reservoir');
end
if size(sas.Qmatrices,1) ~= nInternalUnits
    error('the size of the Qmatrix does not match the size of the reservoir');
end
if size(sas.Rmatrices,1) ~= nInternalUnits
    error('the size of the Rmatrix does not match the size of the reservoir');
end
if size(sas.Smatrices,1) ~= nInternalUnits
    error('the size of the Smatrix does not match the size of the reservoir');
end
if sas.connectivity > 1  || sas.connectivity < 0
    error('the connectivity has to be between 0 and 1');
end
if ~strcmp(sas.learningMode,'offline') &&...
        ~strcmp(sas.learningMode,'online')
    error('learningMode should be either "offline" or "online" ') ;
end
if ~((strcmp(sas.outputActivationFunction,'identity') && ...
        strcmp(sas.inverseOutputActivationFunction,'identity')) || ...
        (strcmp(sas.outputActivationFunction,'tanh') && ...
        strcmp(sas.inverseOutputActivationFunction,'atanh')) || ...
        (strcmp(sas.outputActivationFunction,'sigmoid01') && ...
        strcmp(sas.inverseOutputActivationFunction,'sigmoid01_inv')))  ...
        error('outputActivationFunction and inverseOutputActivationFunction do not match');
end


