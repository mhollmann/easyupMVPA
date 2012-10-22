% Returns an empty 4D dataset (i.e. for fMRI-data).
%
% Author: Maurice Hollmann
% Date  : 08/07
%
% Description:
%
%   [dataset] = getEmpty4DDataset()
%
%   Returns an empty 4D dataset (i.e. for fMRI-data) , that can be "filled" with content:
%   This is aequivalent to getEmpty2DDataset() in 2D case (i.e. eeg data). 
%
%   Returned sructure:
%
%   dataset.type                  = 'dataset4D';
%   dataset.is4D = true;
%   dataset.is2D = false;
%   dataset.dataFilelist          = [];
%   dataset.data                  = [];
%   dataset.chunks                = [];
%   dataset.classIDs              = [];
%   dataset.mask                  = [];
%   dataset.featureSelectionMap   = [];
%   dataset.data_3DNiftiHdr       = [];
%   dataset.processingHistory     = {};
%
% Parameters:
%
% Returns:
%
%   dataset          -> empty dataset
%
% Comments:
%
function [dataset] = getEmpty4DDataset()

  dataset = {};
  
  dataset.type = 'dataset4D';
  dataset.is4D = true;
  dataset.is2D = false;
  
  dataset.dataFilelist          = [];
  dataset.data                  = [];
  dataset.chunks                = [];
  dataset.classIDs              = [];
  dataset.mask                  = [];
  dataset.featureSelectionMap   = [];
  dataset.data_3DNiftiHdr       = [];
  dataset.processingHistory     = {};
  
end