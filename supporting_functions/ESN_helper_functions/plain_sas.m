function internalState = plain_sas(totalstate , sas , varargin)

% PLAIN_SAS computes the new internal states of the SAS by using the simple
% sas equations
%
% input arguments:
% totalstate: the previous totalstate, vector of size 
%     (sas.nInternalUnits + sas.nInputUnits + sas.nOutputUnits) x 1
% sas: the SAS structure
%
% output: 
% internalState: the updated internal state, size esn.nInternalUnits x 1
%
% Version 1.0, August 7, 2017
% Copyright: L. Grigoryeva and J.-P. Ortega


internalState =  feval(sas.reservoirActivationFunction, ...
    totalstate)  ;   

