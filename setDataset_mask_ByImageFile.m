% Set the mask field of 4D datasets by a given nifti-file.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
% This methods sets a mask for the dataset (in example for excluding non-brain voxels). 
% The mask has to be in 3D space. All non-zero elements in mask image are included in mask.
%
% Parameters:
%   dataset        - the datset to set the mask for
%   maskImageFile  - string: the image file (3D image as .nii or .hdr) image is containing non-zero and zero elements - non-zero = include , zero = exclude
%
% Returns:
%   dataset   - the datset with included mask (this mask just contains zeros and ones)
%
% Comments:
%
function [dataset] = setDataset_mask_ByImageFile(dataset, maskImageFile)

  if(~exist('dataset','var') || ~exist('maskImageFile','var'))
      error('Usage of setDataset_mask_ByImageFile: [dataset] = setDataset_mask_ByImageFile(dataset, maskImageFile [string pointing to the file (.nii or .hdr)]');
  end
  
  if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
    tmpMask = dataset.mask;
  else  
    tmpMask = [];
  end
  
  %Load all the data in the filelist and concatenate in dataset.data4D
  disp('INFO: setDataset_mask3D_ByImageFile: Loading mask image file ...');
  disp(maskImageFile);
    
  %check existance of the file
  if(~exist(maskImageFile, 'file'))
    error(['Could not read mask image file: ', maskImageFile]);
  end
    
  %Read the data (nii or hdr)
  %the nii loading function is slightly changed for better memory performance
  dataNii = load_nii(maskImageFile);
  
  %check dimensions
  sizeData = size(dataset.data);
  sizeMask = size(dataNii.img);
  
  if(sizeData(1)~=sizeMask(1) || sizeData(2)~=sizeMask(2) || sizeData(3)~=sizeMask(3))
      error('The dimensions of maskImageFile do not fit the space of dataset.data!');
  end
  dataNii.img(dataNii.img>0)=1;
  dataNii.img(dataNii.img<0)=1;
  dataset.mask = dataNii.img;
  
  
  %check dataset
  if(~checkDataset(dataset))
    dataset.mask = tmpMask;
    error('Could not set field mask, dataset.mask is unchanged!');
  end
    
end