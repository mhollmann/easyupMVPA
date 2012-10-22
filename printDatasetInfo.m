% Prints the content of the dataset on the screen.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   printDatasetInfo(dataset)
%
%   This methods prints the whole info about dataset at the screen.
%
% Parameters:
%   dataset  - the dataset to print info 
%
% Returns:
%
%
% Comments:
%
function printDatasetInfo(dataset)

  if(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset4D'))
    member_printDataset4DInfo(dataset);
  elseif(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset2D'))
    member_printDataset2DInfo(dataset);
  else
    disp('PRINT DATASET INFO: Please check the dataset: field "type" is not defined!');
  end
end %end function printDatasetInfo(dataset)


function member_printDataset4DInfo(dataset)
  if( ~isfield(dataset, 'data') )
    dataset.data = [];
  end
  if( ~isfield(dataset, 'dataFilelist') )
    dataset.dataFilelist = [];
  end
   if( ~isfield(dataset, 'data4D_3DNiftiHdr') )
    dataset.data_3DNiftiHdr = [];
  end
  if( ~isfield(dataset, 'mask') )
    dataset.mask = [];
  end
  
  if( ~isfield(dataset, 'featureSelectionMap') )
    dataset.featureSelectionMap = [];
  end
  if( ~isfield(dataset, 'chunks') )
    dataset.chunks = [];
  end
  if( ~isfield(dataset, 'classIDs') )
    dataset.classIDs = [];
  end

  sizeData   = size(dataset.data);
  sizeMask   = size(dataset.mask);
  
  disp('----------------------------------------')
  disp('*** Info Dataset ***')
  disp(' ');
  disp(['Dataset Type: ', dataset.type]);
  disp(' ');
  if(isfield(dataset,'processingHistory') && ~isempty(dataset.processingHistory))
    disp('processing history:');
    disp(' ');
    for i=1:size(dataset.processingHistory,2)
      disp([' * ',char(dataset.processingHistory(1,i))]);
    end
    disp(' ');
  end
  disp('field data (4D):');
  disp(['              size : ', num2str(sizeData)]);
  disp(['              class: ', class(dataset.data)]);
  disp('field dataFilelist:');
  disp(['              size : ', num2str(size(dataset.dataFilelist))]);
  disp(['              class: ', class(dataset.dataFilelist)]);
  disp('field data_3DNiftiHdr:');
  disp(['              size : ', num2str(size(dataset.data_3DNiftiHdr))]);
  disp(['              class: ', class(dataset.data_3DNiftiHdr)]);
  disp('field mask (3D):');
  if(~isempty(dataset.mask))
    disp(['              size                 : ', num2str(sizeMask)]);
    disp(['              nmb non-zero elements: ', num2str(num2str(sum(dataset.mask(:) > 0)))]);
    disp(['              class                : ', class(dataset.mask)]);
  else
    disp('               No mask is defined for dataset.');
  end
  disp('field featureSelectionMap (3D):');
  if(~isempty(dataset.featureSelectionMap))
    disp(['              size                 : ', num2str(size(dataset.featureSelectionMap))]);
    disp(['              nmb non-zero elements: ', num2str(num2str(sum(dataset.featureSelectionMap(:) > 0)))]);
    disp(['              class                : ', class(dataset.featureSelectionMap)]);
  else
    disp('               No featureSelectionMap is defined for dataset.');
  end
  disp('field chunks:');
  disp(['              length: ', num2str(length(dataset.chunks))]);
  disp(['              class : ', class(dataset.chunks)]);
  disp(num2str(dataset.chunks));
  disp('field classIDs:');
  disp(['              length: ', num2str(length(dataset.classIDs))]);
  disp(['              class : ', class(dataset.classIDs)]);
  disp(num2str(dataset.classIDs));
  disp('----------------------------------------')
  disp(' ');

end %end function member_printDatasetInfo(dataset)




function member_printDataset2DInfo(dataset)

  if( ~isfield(dataset, 'data') )
    dataset.data = [];
  end
  if( ~isfield(dataset, 'mask') )
    dataset.mask = [];
  end
  
  if( ~isfield(dataset, 'featureSelectionMap') )
    dataset.featureSelectionMap1D = [];
  end
  if( ~isfield(dataset, 'chunks') )
    dataset.chunks = [];
  end
  if( ~isfield(dataset, 'classIDs') )
    dataset.classIDs = [];
  end

  sizeData   = size(dataset.data);
  sizemask   = size(dataset.mask);
  
  disp('----------------------------------------')
  disp('*** Info Dataset ***')
  disp(' ');
  disp(['Dataset Type: ', dataset.type]);
  disp(' ');
  if(isfield(dataset,'processingHistory') && ~isempty(dataset.processingHistory))
    disp('processing history:');
    disp(' ');
    for i=1:size(dataset.processingHistory,2)
      disp([' * ',char(dataset.processingHistory(1,i))]);
    end
    disp(' ');
  end
  disp('field data (2D):');
  disp(['              size : ', num2str(sizeData)]);
  disp(['              class: ', class(dataset.data)]);
  disp('field mask (1D):');
  if(~isempty(dataset.mask))
    disp(['              size                 : ', num2str(sizemask)]);
    disp(['              nmb non-zero elements: ', num2str(num2str(sum(dataset.mask > 0)))]);
    disp(['              class                : ', class(dataset.mask)]);
  else
    disp('               No mask is defined for dataset.');
  end
  disp('field featureSelectionMap (1D):');
  if(~isempty(dataset.featureSelectionMap))
    disp(['              size                 : ', num2str(size(dataset.featureSelectionMap))]);
    disp(['              nmb non-zero elements: ', num2str(num2str(sum(dataset.featureSelectionMap > 0)))]);
    disp(['              class                : ', class(dataset.featureSelectionMap)]);
  else
    disp('               No featureSelectionMap1D is defined for dataset.');
  end
  disp('field chunks:');
  disp(['              length: ', num2str(length(dataset.chunks))]);
  disp(['              class : ', class(dataset.chunks)]);
  disp(num2str(dataset.chunks));
  disp('field classIDs:');
  disp(['              length: ', num2str(length(dataset.classIDs))]);
  disp(['              class : ', class(dataset.classIDs)]);
  disp(num2str(dataset.classIDs));
  disp('----------------------------------------')
  disp(' ');
end %end function member_printDatasetInfo(dataset)
