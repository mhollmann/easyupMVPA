% Set the chunks of a dataset by a given vector.
%
% Author: Maurice Hollmann
% Date  : 0/910
%
% Description:
%
%   [dataset] = setDataset_chunks_ByVector(dataset, chunkVector)
%
%   This methods sets the chunks of the given dataset. Chunks can be used to
%   subdivide a dataset. One may for example set a separate chunk for all transition 
%   scans for easy removing after preprocessing. 
%
% Parameters:
%   dataset      - the datset to set the data4D for
%   chunkVector  - a vector holing chunk information
%
% Returns:
%   dataset   - the datset with included chunks
%
% Comments:
%
function [dataset] = setDataset_chunks_ByVector(dataset, chunkVector)
  
  if( ~exist('dataset','var') || ~exist('chunkVector','var') ) 
    error('Usage of setDataset_chunks_ByVector: [dataset] = setDataset_chunks_ByAttribFile(dataset, chunkVector [1xN vector])');
  end

  if(size(chunkVector,2) == 1 && size(chunkVector,1) > size(chunkVector,2))
    dataset.chunks = uint16(chunkVector');
  else
    dataset.chunks = uint16(chunkVector);
  end
  
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_chunks_ByVector: In the current state the dataset is not suitable for further processing, please see messages before!');
  end

end