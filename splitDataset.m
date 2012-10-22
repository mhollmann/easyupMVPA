% Splits a dataset by selecting the elements defined by given vectors for the two result datasets.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
% Lowlevel function for splitting a dataset by giving two non-overlaping vectors. The vectors contain zeros and ones, where 
% one means the sample is included and zero it is not.
%
%
% Parameters:
%   dataset     - the dataset to split
%   vectorDS1   - a vector of the length of number of samples containing zeros and ones, elements with ones will be included in ds1
%   vectorDS2   - a vector of the length of number of samples containing zeros and ones, elements with ones will be included in ds2
%
% Returns:
%   dataset    - the filtered datset
%
% Comments:
%
function [dataset1, dataset2] = splitDataset(dataset, vectorDS1, vectorDS2)
  
  %check if the vectors are overlapping
  
  sumV = vectorDS1+vectorDS2;
  if(~(length(vectorDS1)==length(vectorDS2)) || sum(sumV>1))
   error('The vectors describing the split must be of equal length and non-overlapping!');
  end  
    
  
  if(dataset.is4D)
    
    dataset1 = getEmpty4DDataset();
    dataset2 = getEmpty4DDataset();

    vectorDS1 = logical(vectorDS1);
    vectorDS2 = logical(vectorDS2);

    dataset1.data                  = dataset.data(:,:,:, vectorDS1);
    dataset1.chunks                = dataset.chunks(vectorDS1);
    dataset1.classIDs              = dataset.classIDs(vectorDS1);
    dataset1.mask                  = dataset.mask;
    dataset1.featureSelectionMap   = dataset.featureSelectionMap;
    dataset1.dataFilelist          = dataset.dataFilelist;

    dataset2.data                  = dataset.data(:,:,:, vectorDS2);
    dataset2.chunks                = dataset.chunks(vectorDS2);
    dataset2.classIDs              = dataset.classIDs(vectorDS2);
    dataset2.mask                  = dataset.mask;
    dataset2.featureSelectionMap   = dataset.featureSelectionMap;
    dataset2.dataFilelist          = dataset.dataFilelist;
    
  elseif(dataset.is2D)
    
    dataset1 = getEmpty2DDataset();
    dataset2 = getEmpty2DDataset();

    vectorDS1 = logical(vectorDS1);
    vectorDS2 = logical(vectorDS2);

    dataset1.data                = dataset.data(:,vectorDS1);
    dataset1.chunks              = dataset.chunks(vectorDS1);
    dataset1.classIDs            = dataset.classIDs(vectorDS1);
    dataset1.mask                = dataset.mask;
    dataset1.featureSelectionMap = dataset.featureSelectionMap;

    dataset2.data                = dataset.data(:, vectorDS2);
    dataset2.chunks              = dataset.chunks(vectorDS2);
    dataset2.classIDs            = dataset.classIDs(vectorDS2);
    dataset2.mask                = dataset.mask;
    dataset2.featureSelectionMap = dataset.featureSelectionMap;
    
  else
    disp('SPLIT DATASET: Please check the dataset: "type" is not correctly defined!');
  end
  
  
end