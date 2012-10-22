% Checks a given dataset for correctness (i.e. size of defined matrices).
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   checkDataset(dataset)
%
%   Checks a dataset: Dimensions of subfields etc.
%
% Parameters:
%   dataset   - the datset to check
%
% Returns:
%   ok        - true if dataset is ok, false otherwise 
%
% Comments:
%
function [ok] = checkDataset(dataset)

  if(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset4D'))
    [ok] = member_checkDataset4D(dataset);
  elseif(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset2D'))
    [ok] = member_checkDataset2D(dataset);
  else
    ok = false;
    disp('CHECK DATASET: Please check the dataset: field "type" is not defined!');
  end
  
end




function [ok] = member_checkDataset4D(dataset)
  ok = true;

  if(~isempty(dataset.data))
    if(size(dataset.data,4) <1)
      disp('CHECK DATASET: Please check the size of field data! The 4th dimension is lower than 1! This dimension should be your number of samples (timesteps, scans)!)');
      ok = false;
    end
  end
 
  %check the size of the mask
  if(~isempty(dataset.data) && ~isempty(dataset.mask))
    sizeData = size(dataset.data);
    sizeMask = size(dataset.mask);
    if(sizeData(1)~=sizeMask(1) || sizeData(2)~=sizeMask(2) || sizeData(3)~=sizeMask(3))
      disp('CHECK DATASET: Please check the size of field data and mask! The dimensions of mask do not fit the space of data!');
      ok = false;
    end
  end
  
  
  %check the size of the actual featureSelectionMap
  if(~isempty(dataset.data) && isfield(dataset, 'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
    sizeData4D = size(dataset.data);
    sizefeatureSelectionMap = size(dataset.featureSelectionMap);
    if(sizeData4D(1)~=sizefeatureSelectionMap(1) || sizeData4D(2)~=sizefeatureSelectionMap(2) || sizeData4D(3)~=sizefeatureSelectionMap(3))
      disp('CHECK DATASET: Please check the size of field data and featureSelectionMap! The dimensions of featureSelectionMap do not fit the space of data!');
      ok = false;
    end
  end 
  
  %check length of chunks vs. nmbSamples
  if(~isempty(dataset.chunks) && ~isempty(dataset.data))
    if(size(dataset.chunks,2) ~= size(dataset.data,4))
      disp('CHECK DATASET: Please check dataset! The number of elements in chunks and the 4th dimension in data must be equal (This is your number of samples (timesteps, scans))!');
      ok = false;
    end

  end
  
  %check length of classIDs vs. nmbSamples
  if(~isempty(dataset.classIDs) && ~isempty(dataset.data))
    if(size(dataset.classIDs,2) ~= size(dataset.data,4))
      disp('CHECK DATASET: Please check dataset! The number of elements in classIDs and the 4th dimension must be equal (This is your number of samples (timesteps, scans))!');
      ok = false;
    end
  end

  if(~isempty(dataset.classIDs) && ~isempty(dataset.chunks))
    if(size(dataset.classIDs,2) ~= size(dataset.data,4))
      disp('CHECK DATASET: Please check dataset! The number of elements in classIDs and the number of elements in chunks must be equal (This is your number of samples (timesteps, scans))!');
      ok = false;
    end
  end
end %end function member_checkDataset4D




function [ok] = member_checkDataset2D(dataset)

  ok = true;
  if(~isempty(dataset.data))
    if(size(dataset.data,2) <1)
      disp('CHECK DATASET: Please check the size of field data! The 2 dimension is lower than 1! This dimension should be your number of samples (timesteps, scans)!)');
      ok = false;
    end
  end
  
  %check length of chunks vs. nmbSamples
  if(~isempty(dataset.chunks) && ~isempty(dataset.data))
    if(size(dataset.chunks,2) ~= size(dataset.data,2))
      disp('CHECK DATASET: Please check dataset! The number of elements in "chunks" and the 2nd dimension in data must be equal (This is your number of samples (i.e. timesteps))!');
      ok = false;
    end

  end
  
  %check length of classIDs vs. nmbSamples
  if(~isempty(dataset.classIDs) && ~isempty(dataset.data))
    if(size(dataset.classIDs,2) ~= size(dataset.data,2))
      disp('CHECK DATASET: Please check dataset! The number of elements in "classIDs" and the 2nd dimension must be equal (This is your number of samples (i.e. timesteps))!');
      ok = false;
    end
  end
  
 
end %end function member_checkDataset2D