% Permutes the classIDs in dataset but preserves class sizes.
%
% Author: Maurice Hollmann
% Date  : 11/10
%
% Description:
% This methods permutes the classIDs in the dataset. Thus the function
% destroys the connection between data inherent structure and classes. This
% may be used for permutation tests to determine the guessing level of a classifier.
%
%
% Parameters:
%   dataset   - the datset to permute the classIDs for
% Returns:
%   dataset   - the datset with permuted classIDs
%
% Comments:
%
function [dataset] = permuteClassIDs(dataset)

  if( ~exist('dataset','var'))
     error('Usage of permuteClassIDs: [dataset] = permuteClassIDs(dataset)');
  end
  
  savedVec = dataset.classIDs;
  randIndices = randperm(size(dataset.classIDs,2));
  dataset.classIDs(randIndices) = savedVec;
  
end