% Set the mask (1D or 3D) field of a dataset by a given matrix.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% [dataset] = setDataset_mask_ByMatrix(dataset, maskMatrix)
%
% Description:
% This methods sets a mask  (1D or 3D) for the dataset (both types dataset2D and 4D supported) for excluding i.e. non-brain voxels. 
% If dataset is 4D-type the mask has to be in 3D space. 
% If dataset is 2D-type the mask has to be a Vector. 
% All non-zero elements in mask image are included in mask.
%
% Parameters:
%   dataset   - the datset to set the mask for
%   matrix    - 1D or 3D matrix with zeros and ones (zeros exclude features)
%
% Returns:
%   dataset   - the datset with included mask (this mask just contains zeros and ones)
%
% Comments:
%
function [dataset] = setDataset_mask_ByMatrix(dataset, maskMatrix)
  
  %4D dataset
  if(dataset.is4D)
    
    if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
      tmpMask = dataset.mask;
    else  
      tmpMask = [];
    end

    %check dimensions
    sizeData = size(dataset.data);
    sizeMask = size(maskMatrix);

    if(sizeData(1)~=sizeMask(1) || sizeData(2)~=sizeMask(2) || sizeData(3)~=sizeMask(3))
      error('The dimensions of maskImageFile do not fit the space of dataset.data!');
    end
    maskMatrix(maskMatrix>0)=1;
    maskMatrix(maskMatrix<0)=1;
    dataset.mask = maskMatrix;

    %check dataset
    if(~checkDataset(dataset))
      dataset.mask = tmpMask;
      error('Could not set field mask, dataset.mask is unchanged!');
    end
    
  %2D dataset
  elseif(dataset.is2D)
    
    if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
      tmpMask = dataset.mask;
    else  
      tmpMask = [];
    end

    %check dimensions
    sizeData = size(dataset.data);
    if(sizeData(1)~= length(maskMatrix))
      error('The length of maskMatrix do not fit the nmb of samples in dataset.data2D!');
    end
    maskMatrix(maskMatrix>0)=1;
    maskMatrix(maskMatrix<0)=1;
    
    if(size(maskMatrix,2) == 1 && size(maskMatrix,1) > size(maskMatrix,2))
      dataset.mask = uint8(maskMatrix);
    else
      dataset.mask = uint8(maskMatrix');
    end

    %check dataset
    if(~checkDataset(dataset))
      dataset.mask = tmpMask;
      error('Could not set field mask, dataset.mask is unchanged!');
    end
     
  end
  
  
end