% Set the class-IDs of a dataset by a given vector.
%
% Author: Maurice Hollmann
% Date  : 0/910
%
% Description:
%
%   [dataset] = setDataset_classIDs_ByVector(dataset, classIDsVector)
%
%   This methods sets the classIDs of the given dataset (i.e. [0 0 0 1 1 1]). 
%
% Parameters:
%   dataset         - the datset to set the classIDs for
%   classIDsVector  - A vector holding class IDs
%
% Returns:
%   dataset   - the datset with included classIDs
%
% Comments:
%
function [dataset] = setDataset_classIDs_ByVector(dataset, classIDsVector)
  
  if( ~exist('dataset','var') || ~exist('classIDsVector','var') ) 
    error('Usage of setDataset_classIDs_ByVector: [dataset] = setDataset_chunks_ByAttribFile(dataset, classIDsVector [1xN vector])');
  end

  if(size(classIDsVector,2) == 1 && size(classIDsVector,1) > size(classIDsVector,2))
    dataset.classIDs = uint8(classIDsVector');
  else
    dataset.classIDs = uint8(classIDsVector);
  end
  
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_classIDs_ByVector: In the current state the dataset is not suitable for further processing, please see messages before!');
  end

end