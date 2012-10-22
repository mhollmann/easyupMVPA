% Select a sub-dataset out of a dataset according to given criteria (string or binary vector).
%
% Author: Maurice Hollmann
% Date  : 04/11
%
% Description:
% This methods selects samples from the dataset. The selection criterion 
% is given as a string or a binary vector. It refers always to index, chunks 
% or classIDs if given as a string and always to index if given as a vector.
%
%  Selection given as string:
%  ds = selectSamples(ds, 'index <= 500')
%  ds = selectSamples(ds, 'chunks > 0')
%  ds = selectSamples(ds, 'chunks > 0 & classIDs < 2')
%  ds = selectSamples(ds, 'chunks > 0 | classIDs < 2')
%  ds = selectSamples(ds, 'chunks == 1')
%
%  Selection given as vector:
%  ds = selectSamples(ds, [0 1 1 0 0 0 1]) !!Here the vector must have the length of the number of samples in ds!!
%
%
% Parameters:
%   dataset       - the datset to set the classIDs for
%   criterion     - selection criterion (string or vector) that will be evaluated as given 
%                   (i.e. 'index <= 500' or 'dataset.chunks >= 30' or a vector of ones and zeros)
%
% Returns:
%   dataset   - the datset with excluded samples
%
% Comments:
%
function [dataset] = selectSamples(dataset, criterion)
   
   if( ~exist('dataset','var') || ~exist('criterion','var'))
     error('Usage of selectSamples: [dataset] = selectSamples(dataset, criterion [i.e. "chunks > 0"])');
   end

   index = 1:length(dataset.chunks);
   
   
   %criterion is a vector
   if(isnumeric(criterion))
     
     if(size(criterion,1) > size(criterion,2))
       criterion = criterion';
     end
     
     sizeData = size(dataset.data);
     
     if(length(criterion) ~= sizeData(end) || size(criterion,1) ~= 1)
       disp('ssssssssERROR: selectSamples: Selection criterion not valid!')
       disp('Examples for the usage of selectSamples using a string as criterion:')
       disp('ds = selectSamples(ds, `index <= 500`)');
       disp('ds = selectSamples(ds, `chunks > 0`)');
       disp('ds = selectSamples(ds, `chunks == 2`)');
       disp('ds = selectSamples(ds, `chunks > 0 & classIDs >= 1`)');
       disp('ds = selectSamples(ds, `chunks > 0 | classIDs >= 1`)');
       disp(' ');
       disp('Examples for the usage of selectSamples using a vector as criterion:')
       disp('ds = selectSamples(ds, [0 1 1 0 0 0 1]); %length of selection vector must fit numbers of samples in ds !');
       error('Error selecting samples!');
     else
       criterion = ['[',num2str(criterion),']'];
     end
   
   %criterion is a string
   else
     criterion = strrep(criterion, 'classIDs', 'dataset.classIDs');
     criterion = strrep(criterion, 'chunks', 'dataset.chunks');
   end
   
   try
       sel = eval(criterion);
   catch
       disp('ERROR: selectSamples: Selection criterion not valid!')
       disp('Examples for the usage of selectSamples using a string as criterion:')
       disp('ds = selectSamples(ds, `index <= 500`)');
       disp('ds = selectSamples(ds, `chunks > 0`)');
       disp('ds = selectSamples(ds, `chunks == 2`)');
       disp('ds = selectSamples(ds, `chunks > 0 & classIDs >= 1`)');
       disp('ds = selectSamples(ds, `chunks > 0 | classIDs >= 1`)');
       disp(' ');
       disp('Examples for the usage of selectSamples using a vector as criterion:')
       disp('ds = selectSamples(ds, [0 1 1 0 0 0 1]); %length of selection vector must fit numbers of samples in ds !');
       error('Error selecting samples!');
   end
      
   dataset.chunks(~sel)        = [];
   dataset.classIDs(~sel)      = [];
   
   if(dataset.is4D)
     dataset.data(:,:,:, ~sel) = [];
   elseif(dataset.is2D)
     dataset.data(:, ~sel) = [];
   end
   
   if(isempty(dataset.data) || isempty(dataset.classIDs) || isempty(dataset.chunks))
     error('Selection results in an empty dataset, please check selection criteria.');
   end

   
end