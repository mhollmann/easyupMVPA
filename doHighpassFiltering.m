% Does the high-pass filtering in sample dimension (i.e. per voxel in fMRI-time series).
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = doHighpassFiltering(dataset,  samplingFreqency, [cutoffFrequency])
%
%   This method does the highpass-filtering in time domain. Input is a dataset with 4D Data.
%   The 4th dimension is used for filtering of every single element of the first 3 dimensions.
%   If cutoff frequency is not given 1/128  Hz (128 seconds cycle) is used as cutoff. 
%   The mean of the timecourse is again added at the end of  the algorithm.
%
% Parameters:
%   dataset            - the dataset to filter 
%   samplingFrequency  - in Hz (i.e. for fMRI with TR 2 sec, the samplingFreq is 0.5 )
%   cutoffFrequency    - the cuttoff frequency in Hz (all frequencies below will be removed) [optional, default = 1/128]
%
% Returns:
%   dataset    - the filtered datset
%
% Comments:
%
function [dataset] = doHighpassFiltering(dataset, samplingFrequency, cutoffFreq)
  
  %use Nyquist Theorem? 
  %nyq = 1/(2*1);% 1/ TR*2;
  %[b,a] = butter(2,cutoffFreq/nyq,'high');
  
  if ~exist('dataset','var')
      error('Usage of doHighpassFiltering: [dataset] = doHighpassFiltering(dataset,  samplingFreqency, cutoffFreq[optional, default = 1/128])');
  end
  
  if ~exist('samplingFrequency','var')
      error('Usage of doHighpassFiltering: [dataset] = doHighpassFiltering(dataset,  samplingFreqency, cutoffFreq[optional, default = 1/128])');
  end
  
  if ~exist('cutoffFreq','var')
      cutoffFreq = 1/128;
      disp('INFO: doHighpassFiltering: Setting cuttoff frequency to 1/128 Hz.');
  end
  
  
  % +++ HighPass Filter +++
  %Normalize Frequency (freq values are specified normalized between 0 and 1. 1.0 corresponds to half the sampling frequency f)
  %For 300 Hz and sampling freq 1000 the normalized cutoff is 300/(1000/2)
  normFactor = samplingFrequency*0.5;
  
  % Define Filter
  [b,a] = butter(2,cutoffFreq/normFactor,'high');

  if(~easyupMVPA_getGlobals('quietMode'))
    disp(['Running temporal filtering (high-pass) with cuttoff frequency: ',num2str(cutoffFreq),'. Please wait ...']);
  end
  
  if(dataset.is4D)
    sizeData = size(dataset.data);
    data2D   = reshape(dataset.data, sizeData(1)*sizeData(2)*sizeData(3),sizeData(4));
  elseif(dataset.is2D)
    sizeData = size(dataset.data);
    data2D  = dataset.data;
  end
  
  if(~strcmp(class(data2D), 'double'))
    data2D = double(data2D);
  end

  
  %calculate the mean for adding it afterwards again
  meanData2D = mean(data2D,2);
  
  tic;
  data2D = filtfilt(b,a, data2D')';
  
  
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
  
  dataset.data = single(dataset.data);
  
  t = toc;
  
  histString = ['High-Pass Filtering. SamplingFreq: ', num2str(samplingFrequency), ' CutoffFreq: ', num2str(cutoffFreq)];
  
  if(isfield(dataset,'processingHistory') && ~isempty(dataset.processingHistory))  
    dataset.processingHistory(1,size(dataset.processingHistory,2)+1) = {histString};
  else
    dataset.processingHistory = [{histString}];
  end
  
  if(~easyupMVPA_getGlobals('quietMode'))
    disp(['Done temporal filtering (Time needed was: ',num2str(t), ' sec).']);
  end 
  
end

