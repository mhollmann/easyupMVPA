% Show the histogram of the any data.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
% Parameters:
%
% Returns:
%
% Comments:
%
function showHistogram(dataMatrixIn, maskMatrix, infoString)
  
  if(isempty(findobj('Tag', infoString)))    
    h = figure('Tag' , infoString, 'Name' ,infoString);
  else
    figure(findobj('Tag', infoString));
  end

  histData = dataMatrixIn(:);
  histData = histData(maskMatrix(:)>0);
  xStep = -50:1:50;
  hist(histData, xStep);
  
  
  
end