function [NRMSE, MSE] = compute_NRMSE(estimatedOutput, correctOutput)
% Computes the NRMSE between estimated and correct ESN outputs.
% 
% input arguments:
%
% estimatedOutput: array of size nOutputUnits x sampleLength, containing network
% output data
% correctOutput: array of size nOutputUnits x sampleLength, containing the
% original teacher data. 
%
% output:
% err: a row vector of NRMSE's, each corresponding to one of the output
% dimensions.
%
% If length(correctOutput) > length(estimatedOutput), the first
% elements from correctOutput are deleted. This accounts for cases where
% some (nForgetPoints many) initial transient data points were cancelled
% from estimatedOutput, as occurs in calls to test_esn.
%
% Version 1.0, June 6, 2006, H. Jaeger (as compute_error)
% Revision 1, August 17, 2007, H. Jaeger (renamed to compute_NRMSE,
%                    changed length to size)
% Copyright: Fraunhofer IAIS 2006 / Patents pending
% Revision 1, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: adapted to other conventions

  
nEstimatePoints = size(estimatedOutput, 2) ; 


correctVariance = var(correctOutput, 0, 2) ; 
meanerror = sum((estimatedOutput - correctOutput).^2, 2)/nEstimatePoints ; 
MSE = meanerror;
NRMSE = sqrt(meanerror./correctVariance); 
