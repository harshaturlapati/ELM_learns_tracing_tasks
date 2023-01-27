function plot_states(stateMatrix, states, nPoints, figNr, titletext)

% PLOT_STATES plots the internal states of the esn
%
% inputs:
% stateMatrix = matrix of size (nInputUnits + nInternalUnits ) x (nTrainingPoints)
% states = vector of size 1 x n , containing the indices of the internal
% units we want to plot
% nPoints = natural number containing the number of points to plot
% figNr: either [] or an integer. If [], a new figure is created, otherwise
% the plot is displayed in a figure window with number figNr
% titletext: a string which is displayed as title over first panel
%
% example  : plot_states(stateMatrix,[1 2 3 4],200) plots the first 200
% points from the traces of the first 4 units

%
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, June 6, 2006, H. Jaeger
% Revision 2, July 2, 2007, H. Jaeger
% Revision 3, Aug 17, 2007, H. Jaeger
% Revision 1, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: adapted to new
% conventions and new options


if isempty(figNr)
  figure ; clf;
else
  figure(figNr); clf;
end


nStates = length(states) ;

xMax = ceil(sqrt(nStates)) ;
yMax = ceil(nStates /xMax);

for iPlot = 1 : nStates
  subplot(xMax,xMax,iPlot) ;
  plot(stateMatrix(states(1,iPlot),1:nPoints));
  if iPlot == 1
    title(titletext);
  end
end

