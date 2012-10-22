% Prints the content of the given data-splitter on the screen.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%
%   printDataSplitter(dataSplitter)
%
%   Simple print of dataSplitter-struct.
%   For information about the structure of dataSplitter see function "getDataSplitter".
%
% Parameters:
%   dataSplitter  - dataSplitter to print
%
% Returns:
%
%
% Comments:
%
function printDataSplitter(dataSplitter)

  disp('******* DataSplitter ******')
  disp(['Type: ', dataSplitter.splitterType]);
  disp('SplitMatrix (0=not used, 1=test, 2=training) : ');

  for i=1:size(dataSplitter.splitMatrix,1)
    disp(num2str(dataSplitter.splitMatrix(i,:)));
  end

  disp('***************************');

end