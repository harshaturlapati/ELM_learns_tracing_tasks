function [outputWeights, outputIntercept] = pseudoinverse(trainInput, trainTeacher, stateCollection, esn)

% PSEUDOINVERSE computes the outputWeights using the standard pseudoinverse
% operation
%
% inputs:
% trainInput = input time series of size nInputUnits x (nTrainingPoints-nForgetPoints) 
% trainTeacher = teacher time series of size nOutputUnits x (nTrainingPoints-nForgetPoints) 
% stateCollectMat = matrix of size nInternalUnits x (nTrainingPoints-nForgetPoints)
%
% output:
% outputWeights = ifRegressInput = 0 then it is matrix of size nInternalUnits x nOutputUnits containing the learnt weights
%   ifRegressInput = 1 then it is matrix of size (nInternalUnits + nInputUnits) x nOutputUnits containing the learnt weights

% outputIntercept = vector of size nOutputUnits x 1 containing the learnt
%   intercept


% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: different
% conventions. Generalized to multivariate inputs and outputs

if esn.ifRegressInput == 0;
    x = stateCollection;
else
    x =[stateCollection; trainInput];
end
y = feval(esn.inverseOutputActivationFunction, trainTeacher);

[N, T] = size(x);
% get the mean of the reservoir output
muX = mean(x, 2);
muY = mean(y, 2);
% compute the GAMMA(0) of the reservoir output
Gamma0 = (1/T) * (x * x') - muX * muX';
% get the Cov(x,y) of the reservoir output and the teaching signal
Cov_XY = (1/T) * x * y' - muX * muY';
% get the reservoir output matrix and intercept
outputWeights = pinv((Gamma0 + esn.ridgeLambda * eye(N))) * Cov_XY;
outputIntercept = muY- outputWeights(:)' * muX;


    
    
    
    