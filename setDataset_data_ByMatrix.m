% Set the data field of a datasets by a given matrix.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = setDataset_data_ByMatrix(dataset, matrix2D)
%
%   This methods sets the data (2D or 4D) for the given dataset and checks the outcome. 
%
% Parameters:
%   dataset     - the datset to set the data2D field
%   dataMatrix  - a 2D or 4D matrix containing the whole data (nmbFeatures x nmbSamples)
%
% Returns:
%   dataset   - the datset with included data2D
%
% Comments:
%
function [dataset] = setDataset_data_ByMatrix(dataset, dataMatrix)

  if(~isfield(dataset, 'type'))
    error('Given dataset has not defined the field "type" !');
  end
  if( ~exist('dataset','var') || ~exist('dataMatrix','var') ) 
    error('Usage of setDataset_data_ByMatrix: [dataset] = setDataset_data_ByMatrix(dataset, dataMatrix)');
  end
  if(isempty(dataMatrix))
    error('Given matrix2D is empty!');
  end  
  
  dataset.data = dataMatrix;
  
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_data_ByMatrix: In the current state the dataset is not suitable for further processing, please see messages before!');
  end

end