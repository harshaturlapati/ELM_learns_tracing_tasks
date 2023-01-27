function plot_sequence(teacherSequence, predictedSequence, nPoints, caption, figNr)
% PLOT_SEQUENCE plots each dimension from the both the teacherSequence and
% the sequence predicted by the network. The first nPoints values from each
% sequence are plotted. 
%
% inputs: 
% teacherSequence: matrix of size nDataPoints x nDimensions
% predictedSequence: matrix of size nDataPoints x nDimensions
% nPoints: a natural number. 
%           The first nPoints from the input time series are plotted
% caption: a string containing the caption of the figure

%
% Created April 30, 2006, D. Popovici
% Copyright: Fraunhofer IAIS 2006 / Patent pending
% Revision 1, Feb 23, 2007, H. Jaeger
% Revision 2, Feb 13, 2017, L. Grigoryeva and J.-P. Ortega: adapted to new
% conventions and new options

if nargin == 3
    caption = '' ; 
end


if isempty(figNr)
  figure ; clf;
else
  figure(figNr); clf;
end

if isempty(nPoints)
  nPoints = size(teacherSequence, 2);
end

nDimensions = length(teacherSequence(:,1)) ; 

for iPlot = 1 : nDimensions
    subplot(nDimensions,1,iPlot) ; 
    %%%% set the caption of the figure    
    title(caption);
    hold on ; 
    plot(teacherSequence(iPlot, 1:nPoints),'r') ; 
    plot(predictedSequence(iPlot, 1:nPoints)) ; 
end
