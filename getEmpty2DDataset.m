% Returns an empty 2D dataset (i.e. for EEG-data).
%
% Author: Maurice Hollmann
% Date  : 08/07
%
% Description:
%
%   [dataset] = getEmpty2DDataset()
%
%   Returns an empty 2D dataset (i.e. for EEG-data) , that can be "filled" with content.
%   This is aequivalent to getEmpty4DDataset() in 4D case (i.e. fMRI time series). 
%
%   Returned sructure:
%
%   dataset.type                   = 'dataset2D';
%   dataset.is4D = false;
%   dataset.is2D = true;
%   dataset.data2D                 = [];
%   dataset.chunks                 = [];
%   dataset.classIDs               = [];
%   dataset.mask                   = [];
%   dataset.featureSelectionMap    = [];
%   dataset.processingHistory      = {};
%
% Parameters:
%
% Returns:
%
%   dataset          -> empty dataset
%
% Comments:
%
function [dataset] = getEmpty2DDataset()

  dataset = {};
  
  dataset.type = 'dataset2D';
  dataset.is4D = false;
  dataset.is2D = true;

  dataset.data                   = [];
  dataset.chunks                 = [];
  dataset.classIDs               = [];
  dataset.mask                   = [];
  dataset.featureSelectionMap    = [];
  dataset.processingHistory      = {};
  
end