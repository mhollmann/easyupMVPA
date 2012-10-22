% Removes linear trends in sample dimension (i.e. per voxel over samples in fMRI-time series). 
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = doLinearDetrending(dataset, [breakpoints])
%
%   This method does a linear Detrending of the data. Every 3D element is processed for itself.
%   Fitting is done by a least squares fit of a linear function. The mean of the timecourse is
%   again added at the end of  the algorithm.
%
%   Usage of breakpoints: 
%   If i.e. your dataset is from 2 different runs (distinct in time) set a breakpoint at the
%   end of session 1 (i.e. sample 100) and at the begining of session 2 (i.e. sample 101):  
%   [ds] = doLinearDetrending(ds, [100 101]);
%                         
%
% Parameters:
%   dataset     - the dataset to detrend
%   breakpoints - [optional] scalar or vector of breakpoints in signal i.e. if two independent datasets are combined
%
% Returns:
%   dataset    - the detrended datset
%
% Comments:
%
function [dataset] = doLinearDetrending(dataset, breakpoints)
  
  if(~exist('dataset','var'))
    error('Usage of doLinearDetrending: [dataset] = doLinearDetrending(dataset, breakpoints [optional - scalar or vector of breaks in signal])');
  end
  
  if(exist('breakpoints','var') && ~isnumeric(breakpoints))
    error('Usage of doLinearDetrending: [dataset] = doLinearDetrending(dataset, breakpoints [optional - scalar or vector of breaks in signal])');
  end
  
  if(dataset.is4D)
    sizeData = size(dataset.data);
    data2D   = reshape(dataset.data, sizeData(1)*sizeData(2)*sizeData(3),sizeData(4));
  elseif(dataset.is2D)
    sizeData = size(dataset.data);
    data2D  = dataset.data;
  end
  
  
  if(~easyupMVPA_getGlobals('quietMode'))
    disp('Running linear detrending. Please wait ...');
  end
  
  tic;
  if(~strcmp(class(data2D), 'single') )
   data2D = single(data2D);
  end
  
  %calculate the mean for adding it afterwards again
  meanData2D = mean(data2D,2);
  
  if(exist('breakpoints','var'))
    data2D = detrend(data2D', 'linear', breakpoints)';
    histString = ['Linear Detrending. Breakpoints: ', num2str(breakpoints)];
  else
    data2D = detrend(data2D', 'linear')';
    histString = 'Linear Detrending. Breakpoints: none';
  end
  
  %add the mean again
  if(dataset.is4D)
    meanData3D = reshape(meanData2D, sizeData(1),sizeData(2),sizeData(3));
    dataset.data = reshape(data2D, sizeData(1), sizeData(2), sizeData(3), sizeData(4));
    for i=1:sizeData(4)
      dataset.data(:,:,:,i) = dataset.data(:,:,:,i)+meanData3D;
    end
  elseif(dataset.is2D)
    dataset.data = data2D;
    for i=1:sizeData(2)
      dataset.data(:,i) = dataset.data(:,i)+meanData2D;
    end
  end
  
  t = toc;
  
  if(isfield(dataset,'processingHistory') && ~isempty(dataset.processingHistory))  
    dataset.processingHistory(size(dataset.processingHistory,2)+1) = {histString};
  else
    dataset.processingHistory = {histString};
  end
  
  if(~easyupMVPA_getGlobals('quietMode'))
    disp(['Done linear detrending (Time needed was: ',num2str(t), ' sec).']);
  end
  
end