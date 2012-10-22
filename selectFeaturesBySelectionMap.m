% Apply a featureSelection map that is stored in a dataset to this dataset. 
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% [dataset, data2D] = selectFeaturesBySelectionMap(dataset)
%
% Description:
% Lowlevel function for extracting training data from 
%  (1) a 4D data matrix and a given 3D selection map (for dataset4D)
% or from 
%  (2)a 2D matrix and a featureSelectionMap1D (for dataset2D).
%
% This Function uses the featureSelection info of the actual dataset! 
% If dataset.featureSelectionMap3D or dataset.featureSelectionMap1D is not set it returns an error!
%
% The returned data is 2D and already excluded features which are zero in selectionMap. That means
% the returned matrix has not the same number of elements as input data. For projecting back into
% original feature space one needs to set the elements according to the featureMap which is still 
% in original feature space:  arrayInOrigSpace(featureMap>0)= data2D.
%
% Parameters:
%   dataset     - the dataset from that features should be selected
%
% Returns:
%   dataset     - the datset with included 'featureSelectionMap3D' or 'featureSelectionMap1D'
%   data2D      - the selected features as 2D matrix, the backward mapping can be done using the selectionMap
%
% Comments:
%
function [dataset, data2D] = selectFeaturesBySelectionMap(dataset)
  
  if(dataset.is4D)
    
     if(~isfield(dataset,'featureSelectionMap') || isempty(dataset.featureSelectionMap))
      error('The field dataset.featureSelectionMap is not set. It is not possible to select features without this field.');
     end
     
     numelFselMap = numel(dataset.featureSelectionMap(dataset.featureSelectionMap > 0));
     
     %more memory less computing time...
     selMap   = repmat(dataset.featureSelectionMap, [1, 1, 1, size(dataset.data,4)]);
     data2D   = double(reshape(dataset.data(selMap>0), numelFselMap, size(dataset.data,4))');
     
     %less memory more computing time
%      data2D   = zeros(size(dataset.data,4), numelFselMap); 
%      for i=1:size(dataset.data,4)
%        tmp            = dataset.data(:,:,:,i);
%        data2D(i,:)    = tmp(dataset.featureSelectionMap>0);
%      end
     
  elseif(dataset.is2D)
    
     if(~isfield(dataset,'featureSelectionMap') || isempty(dataset.featureSelectionMap))
      error('The field dataset.featureSelectionMap is not set. It is not possible to select features without this field.');
     end
     
     numelFselMap = numel(dataset.featureSelectionMap(dataset.featureSelectionMap > 0));

     %more memory less computing time...
     selMap   = repmat(dataset.featureSelectionMap, [1, 1, size(dataset.data,2)]);
     data2D   = double(reshape(dataset.data(selMap>0), numelFselMap, size(dataset.data,2))');
    
     %less memory more computing time
%      data2D   = zeros(size(dataset.data, 2), numel(find(dataset.featureSelectionMap > 0)));
%      for i=1:size(dataset.data,2)
%        tmp            = dataset.data(:,i);
%        data2D(i,:)    = tmp(dataset.featureSelectionMap>0);
%      end
     
  else
    error('selectFeaturesBySelectionMap: Please check the dataset: "type" is not corretly defined!');
  end
  
end