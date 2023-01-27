function plot_outputWeights(outweights, figNr)
% plot the output weights (of a trained ESN). For each output channel, a
% new color is used.
%
% input arguments:
% outweights: a weight matrix of size  (nInputUnits +
% nInternalUnits) x nOutputUnits
% 
% figNr: either [] or an integer. If [], a new figure is created, otherwise
% the plot is displayed in a figure window with number figNr
%
% Created June 6, 2006, H. Jaeger
% Revision 1, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: adapted to new
% conventions and new options


if isempty(figNr)
  figure ; clf;
else
  figure(figNr); clf;
end

plot(outweights');

title('output weights'); 
