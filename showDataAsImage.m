% Show the data that is given as image.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%  Shows image data as simple color-coded mosaic-view. 
%  If the data is 2D it is shown as image as given.
%  If the data is 3D it is expected that z are slices and these are shown as mosaic in one figure.
%  If the data is 4D and there is a z-Index given this 3D block will be choosen and shown as mosaic in one figure.
%  If the data is 4D and there is no z-Index given the first image in Z will be choosen for showing.
%
% Parameters:
%   dataIn     - the dataset from that features should be selected
%   infoString - a string that appears in the title of the figure window 
%   tIndex     - (optional) the index in t-dimension (for 4D data)
%
% Returns:
%
% Comments:
%
function showDataAsImage(dataIn, infoString, tIndex)
  
  if( ~exist('dataIn','var') || ~exist('infoString','var') ) 
    error('Usage of showDataAsImage: showDataAsImage(dataIn [2D, 3D, or 4D matrix], infoString, tIndex [optional - if 4D is given this can give the index in the 4th dim])');
  end

  sizeDataIn = size(dataIn);
  if(numel(sizeDataIn) < 2 || numel(sizeDataIn) >4)
     error('Usage of showDataAsImage: showDataAsImage(dataIn [2D, 3D, or 4D matrix], infoString, tIndex [optional - if 4D is given this can give the index in the 4th dim])');
  end
  
  if(isempty(dataIn))
    error('Trying to display an empty matrix!');
  end
  
  %if the user has given an info string but no tIndex
  if(exist('tIndex', 'var') && ~exist('infoString','var') && ischar(tIndex))
    infoString = tIndex;
    tIndex = 1;
  end
  
  if( ~exist('infoString','var') ) 
    infoString = '';
  end
  
  titleAddText = '';
  
  %simply show the data if it is 2D
  if(numel(sizeDataIn)==2)
    if(isempty(findobj('Tag', infoString)))    
      h = figure('Tag' , infoString, 'Name' ,infoString);
    else
      figure(findobj('Tag', infoString));
    end
    imagesc(dataIn(:,:));
    return;
  end
  
  %select the data to view if 4D
  if(numel(sizeDataIn)==4)
    if(exist('tIndex','var') && tIndex <= sizeDataIn(4) && tIndex >= 1)
      dataIn = dataIn(:,:,:,tIndex);
    else
      disp('INFO: showDataAsImage: Choosing first element in 4th dimension for viewing.');
      dataIn = dataIn(:,:,:,1);
    end
  end
  
  %show 3D data
  sizeDataIn     = size(dataIn);
  nmbTilesPerDim = ceil(sqrt(sizeDataIn(3)));
  
  if(isempty(findobj('Tag', infoString)))
    h = figure('Tag', infoString, 'Name' ,infoString);
  else
    figure(findobj('Tag', infoString));
  end
  
  subplot(nmbTilesPerDim,nmbTilesPerDim,1);
  tmpIndex = ceil(sizeDataIn(3)/2);
  
  imagesc(dataIn(:,:,tmpIndex));
  climScaled = get(gca,'CLim');
  
  for i=1:sizeDataIn(3)
    subplot(nmbTilesPerDim,nmbTilesPerDim,i);
    imagesc(dataIn(:,:,i));
    set(gca, 'CLim', climScaled)
    set(gca, 'XTick', [], 'YTick', []);
  end
  
end