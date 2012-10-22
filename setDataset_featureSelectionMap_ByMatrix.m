% Set the featureSelectionMap (1D or 3D) field of a dataset by a given matrix.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = setDataset_featureSelectionMap_ByMatrix(dataset, mapMatrix)
%
%   This high-level function sets the featureSelectionMap of a dataset. This map can be used to select features in 3D or 1D space.
%   The parameter mapMatrix is a 3D matrix (in dataset4D) or 1D matrix (in dataset2D).
%   If a mask is defined the resulting selection-map in the dataset will contain just ones at coordinates that fall into the mask.
%   The mapMatrix can contain non-zero elements all these will be included.
%
%
% Parameters:
%   dataset     - the datset to set the classIDs for
%   mapMatrix   - a 3D or 1D matrix , containing non-zero and zero elements - non-zero = include , zero = exclude
%
% Returns:
%   dataset   - the datset with included featureSelectionMap3D (this mask just contains zeros and ones)
%
% Comments:
%
function [dataset] = setDataset_featureSelectionMap_ByMatrix(dataset, mapMatrix)

  if( ~exist('dataset','var') || ~exist('mapMatrix','var') ) 
    error('Usage of setDataset_featureSelectionMap_ByMatrix: [dataset] = setDataset_featureSelectionMap_ByMatrix(dataset, mapMatrix)');
  end
  
  %4D type
  if(dataset.is4D)
  
    %check dimensions
    sizeData = size(dataset.data);
    sizeMask = size(mapMatrix);
    if(sizeData(1)~=sizeMask(1) || sizeData(2)~=sizeMask(2) || sizeData(3)~=sizeMask(3))
        error('The dimensions of feature selection map do not fit the space of dataset.data!');
    end

    if(isfield(dataset, 'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
      tmpMap = dataset.featureSelectionMap;
    else  
      tmpMap = [];
    end

    mapMatrix(mapMatrix>0)=1;
    mapMatrix(mapMatrix<0)=1;
    dataset.featureSelectionMap = uint8(mapMatrix);

    %combine an existing mask3D and given map (just voxels in both will survive)
    if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
     dataset.featureSelectionMap(dataset.mask==0) = 0;
    end
  
  %2D type
  elseif(dataset.is2D)
    
    %check dimensions
    sizeData = size(dataset.data);
    sizeMask   = size(mapMatrix);
    if(sizeData(1)~=sizeMask(1))
        error('The dimensions of mapMatrix do not fit the space of dataset.data!');
    end

    if(isfield(dataset, 'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
      tmpMap = dataset.featureSelectionMap;
    else  
      tmpMap = [];
    end

    mapMatrix(mapMatrix>0)=1;
    mapMatrix(mapMatrix<0)=1;
    dataset.featureSelectionMap = uint8(mapMatrix);

    %combine an existing mask and given map (just voxels in both will survive)
    if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
     dataset.featureSelectionMap(dataset.mask==0) = 0;
    end
    
  else
    error('Please check the dataset: field "type" is not defined!');
  end    
    
  %check dataset
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_featureSelectionMap_ByMatrix: In the current state the dataset is not suitable for further processing, please see messages before!');
  end
  
end