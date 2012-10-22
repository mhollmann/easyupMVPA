% Does Z-Scoring in sample dimension (i.e. per voxel over samples in fMRI-time series). 
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%    [dataset] = doZScoring(dataset, [chunkwise])
%
%    This high-level function does the z-scoring of the data. 
%    This means the data in sample domain is set to x = x - mean(x)./ stdev(x).
%    Sample domain means that for each element in feature-space (i.e. voxel in fMRI) 
%    the course of this element over the samples is used for Z-Scoring.
%
%    This is necessary to compare the results of several subjects.
%
% Parameters:
%   dataset    - the dataset to z-score 
%   chunkwise  - [optional] true if z-scoring should be done per chunk, default = false
%
% Returns:
%   dataset    - the filtered datset
%
% Comments:
%
function [dataset] = doZScoring(dataset, chunkwise)
  
  if ~exist('dataset','var')
      error('Usage of doZScoring: [dataset] = doZScoring(dataset, chunkwise [optional, true or false, default = false])');
  end

  if (~exist('chunkwise','var'))
    chunkwise = false;
  end
  
  if(~easyupMVPA_getGlobals('quietMode'))
    disp('Running Z-Scoring. Please wait ...');
    tic;
  end
  
  
  %*********** 4D ***********
  if(dataset.is4D)
  
    if(chunkwise)

      %how many different unique chunks are in the data
      uniqueChunkIDs  = unique(dataset.chunks);
      
      qmTmp = easyupMVPA_getGlobals('quietMode');
      easyupMVPA_setGlobals('quietMode', true);

      %loop over uniqe chunks
      for i=1:length(uniqueChunkIDs)
        tmpDS = selectSamples(dataset, ['chunks==', num2str(uniqueChunkIDs(i))]);
        tmpDS = doZScoring(tmpDS, false);
        if(i == 1)
          dsAll = tmpDS;
        else
          dsAll = concatDatasets(dsAll, tmpDS);
        end
      end
      dataset = dsAll;
      easyupMVPA_setGlobals('quietMode', qmTmp);
      
    else
      sizeData = size(dataset.data);
      data2D     = reshape(dataset.data, sizeData(1)*sizeData(2)*sizeData(3),sizeData(4));
      if(~strcmp(class(data2D), 'single'))
        data2D = single(data2D);
      end
      %Column wise Z-Scoring
      data2D = zscore(data2D')';
      dataset.data = reshape(data2D,sizeData(1),sizeData(2),sizeData(3),sizeData(4));
      %set just data that is included in mask if mask is defined
      if(isfield(dataset, 'mask') && ~isempty(dataset.mask))
        for i=1:sizeData(4)
          ds = dataset.data(:,:,:,i);
          ds(dataset.mask==0) = 0;
          dataset.data(:,:,:,i) = ds;
        end
      end
    end%end if chunkwise
    
  %******** 2D *******
  elseif(dataset.is2D)
    
    if(chunkwise)
      
      %how many different unique chunks are in the data
      uniqueChunkIDs  = unique(dataset.chunks);
      
      qmTmp = easyupMVPA_getGlobals('quietMode');
      easyupMVPA_setGlobals('quietMode', true);
      %loop over uniqe chunks
      for i=1:length(uniqueChunkIDs)
        tmpDS = selectSamples(dataset, ['chunks==', num2str(uniqueChunkIDs(i))]);
        tmpDS = doZScoring(tmpDS, false);
        if(i == 1)
          dsAll = tmpDS;
        else
          dsAll = concatDatasets(dsAll, tmpDS);
        end
      end
      dataset = dsAll;
      easyupMVPA_setGlobals('quietMode', qmTmp);
      
    else
      
      sizeData = size(dataset.data);
      data2D = dataset.data;
      if(~strcmp(class(data2D), 'single'))
        data2D = single(data2D);
      end
      %Column wise Z-Scoring
      data2D = zscore(data2D')';
      dataset.data = data2D;

    end%end if chunkwise
    
  else
    disp('Z-SCORING: Please check the dataset: type is not correctly defined!');
  end

  if(~easyupMVPA_getGlobals('quietMode'))
    t = toc;
    disp(['Done Z-Scoring (Time needed was: ',num2str(t), ' sec).']);
  end
   
end

