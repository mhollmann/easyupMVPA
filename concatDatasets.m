% Concatenates two datasets (appending in sample dimension).
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
% This method concats two datasets. It must be ensured that the data is in the same space!
% If one dataset has 10 samples and the second one has 15 samples the resulting dataset has 
% 25 samples...
%
% If any dataset provides a mask, it will be used. If both do, the mask of the first one is used.
% The concatenation order will be like order of input.
%
% Parameters:
%   dataset1    - the first datset
%   dataset2    - the second datset
%
% Returns:
%   dataset     - the datset that is the concatenation of ds1 and ds2
%
% Comments:
%
function [dataset] = concatDatasets(dataset1, dataset2)
  
  if(dataset1.is2D &&  dataset2.is2D)
    dataset = dataset1;
    if(isempty(dataset1.mask) && ~isempty(dataset2.mask))  
      dataset.mask = dataset2.mask;
    end
    dataset.data   = cat(2, dataset.data, dataset2.data);
  elseif(dataset1.is4D &&  dataset2.is4D)
    
    dataset = dataset1;
    if(isempty(dataset1.mask) && ~isempty(dataset2.mask))  
      dataset.mask = dataset2.mask;
    end
    dataset.data   = cat(4, dataset.data, dataset2.data);
  else
    error('CONCAT DATASETS: Please check the dataset: dataset-type is not defined correctly OR input datasets are not of the same type!');
  end
    
  dataset.chunks   = cat(2, dataset.chunks, dataset2.chunks);
  dataset.classIDs = cat(2, dataset.classIDs, dataset2.classIDs);
   
end

