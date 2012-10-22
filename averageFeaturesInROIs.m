% Set the mask (1D or 3D) field of a dataset by a given matrix.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% [dataset] = averageFeaturesInROIs(dataset, roiMap)
%
% Description:
%   Averages the features according to given roiMap. The number of samples is not
%   influenced. The number of feature space dimensions is reduced to the number of 
%   unique indices in the given roiMap.
%
%
% Parameters:
%   dataset   - the original datset
%   roiMap    - 1D or 3D matrix with integers (0=Background, 1=roi1, 2= roi2, ...)
%
% Returns:
%   dataset   - the original input dataset
%   dataset2D - a 2D dataset with feature values representing the average signal of rois
%
% Comments:
%
function [dataset, dataset2D] = averageFeaturesInROIs(dataset, roiMap)

  sizeRoiMap = size(roiMap);

  %4D dataset
  if(dataset.is4D)
    %check dimensions
    sizeInputData = size(dataset.data);
    if(sizeInputData(1)~=sizeRoiMap(1) || sizeInputData(2)~=sizeRoiMap(2) || sizeInputData(3)~=sizeRoiMap(3))
      error('The dimensions of given roiMap (4D) do not fit the space of dataset.data!');
    end
    
    %variable featureData is 2D: [nmbFeatures, nmbSamples]
    featureData = reshape(dataset.data, sizeInputData(1)*sizeInputData(2)*sizeInputData(3),sizeInputData(4));
    
  elseif(dataset.is2D)
    %check dimensions
    sizeInputData = size(dataset.data);
    if(sizeInputData(1)~= length(sizeRoiMap))
      error('The given roiMap (1D) does not fit the nmb of features in dataset.data!');
    end
    
    %variable featureData is 2D: [nmbFeatures, nmbSamples]
    featureData = dataset.data;
    
  end
  
  %reshape roiMap to 2D
  roiMap  = reshape(roiMap, sizeRoiMap(1)*sizeRoiMap(2)*sizeRoiMap(3), 1);  
  
  %extract number of ROIs in dataset
  roiLabels = unique(roiMap(:));
  
  nmbROIs = numel(roiLabels)-1;
  disp('Number of ROIs:');
  disp(nmbROIs);
  
  newFeatureData = zeros(nmbROIs, sizeInputData(end));
  roiIndex = 1;
  
  %loop over roi indices
  for i=1:numel(roiLabels)
    
    if(roiLabels(i)~=0)
      %select all elements with given roi index
      tmpArr = featureData(roiMap==roiLabels(i),:);
      newFeatureData(roiIndex, :) = mean(tmpArr);
      roiIndex = roiIndex+1;
    end
    
  end
  
  %this one will be returned
  dataset2D = getEmpty2DDataset();
  dataset2D = setDataset_chunks_ByVector(dataset2D, dataset.chunks);
  dataset2D = setDataset_classIDs_ByVector(dataset2D, dataset.classIDs);
  dataset2D = setDataset_data_ByMatrix(dataset2D, newFeatureData);
  
  printDatasetInfo(dataset2D);
   
end



