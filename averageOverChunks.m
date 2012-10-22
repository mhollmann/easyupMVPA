% Averages data in dataset over the chunks with equal ID.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%  This method does an averaging over chunks. That means all samples belonging to the same chunk are
%  averaged in sample dimension (e.g. over time if a sample is a fMRI scan). The number of samples in
%  the result dataset is the number of unique chunk is in input-chunks. 
%
%  For averaging it just makes sense if elements with the SAME chunk ID have the SAME class label. This
%  function does not explicitely control for that!
%
% Parameters:
%   dataset     - the dataset for averaging (types dataset2D or dataset4D)
%
% Returns:
%   dataset    - the dataset with samples averaged over chunks
%
% Comments:
%
function [dataset] = averageOverChunks(dataset)
  
  if ~exist('dataset','var')
    error('Usage of averageOverChunks: [dataset] = averageOverChunks(dataset)');
  end
  
  minChunkID      = min(dataset.chunks);
  maxChunkID      = max(dataset.chunks);
  uniqueChunkIDs  = unique(dataset.chunks);
  
  if(dataset.is2D)
    
      size2d      = size(dataset.data);
      newData2D   = zeros(size2d(1), length(uniqueChunkIDs));
      newChunks   = zeros(1,length(uniqueChunkIDs), 'uint16');
      newClassIDs = zeros(1,length(uniqueChunkIDs), 'uint16');
      d2DIndex    = 1;

      %loop over every chunk in dataset
      %check if chunk id does exist
      for i=minChunkID:maxChunkID

        chunkRefs = dataset.chunks == i;

        %is id inside chunks?
        if(sum(chunkRefs) > 0)
          %average over this chunk id
          avg = mean(dataset.data(:,chunkRefs),2);
          newData2D(:,d2DIndex) = avg;
          %set id for new chunks
          tmpChunks = dataset.chunks(chunkRefs);
          newChunks(d2DIndex)   = tmpChunks(1);
          %set id for new class id
          tmpClassIDs = dataset.classIDs(chunkRefs);
          newClassIDs(d2DIndex) = tmpClassIDs(1);
          d2DIndex = d2DIndex +1;
        end

      end
      dataset.data = newData2D;
  
  elseif(dataset.is4D)
    
      size4d      = size(dataset.data);
      newData4D   = zeros(size4d(1), size4d(2), size4d(3), length(uniqueChunkIDs));
      newChunks   = zeros(1,length(uniqueChunkIDs), 'uint16');
      newClassIDs = zeros(1,length(uniqueChunkIDs), 'uint16');
      d4DIndex    = 1;

      %loop over every chunk in dataset
      %check if chunk id does exist
      for i=minChunkID:maxChunkID

        chunkRefs = dataset.chunks == i;

        %is id inside chunks?
        if(sum(chunkRefs) > 0)
          %average over this chunk id
          avg = mean(dataset.data(:,:,:,chunkRefs),4);
          newData4D(:,:,:,d4DIndex) = avg;
          %set id for new chunks
          tmpChunks = dataset.chunks(chunkRefs);
          newChunks(d4DIndex)   = tmpChunks(1);
          %set id for new class id
          tmpClassIDs = dataset.classIDs(chunkRefs);
          newClassIDs(d4DIndex) = tmpClassIDs(1);
          d4DIndex = d4DIndex +1;
        end

      end
      dataset.data = newData4D;
    
  else 
    
     error('AVERAGE OVER CHUNKS: Please check the dataset: dataset-field "type" is not defined OR input dataset type is invalid!');
  end
 
  %set chunks and classIDs
  dataset.chunks   = newChunks;
  dataset.classIDs = newClassIDs;

  checkDataset(dataset);
  
end % end function