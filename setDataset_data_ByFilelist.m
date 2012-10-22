% Set the data field of a dataset by a given list of nifti-files.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = setDataset_data_ByFilelist(dataset, fileList)
%
%   This methods sets the data (timeseries of 3D images) for the dataset. The given fileList
%   can be a character array or a cellstr array. The fileList is expected to contain the .hdr or .nii files
%   in nifti-format.
%
%   The fileList can contain either 3D images or 4D images (but not in a mixed fashion).
%
% Parameters:
%   dataset   - the datset to set the data field and data_3DNiftiHdr for
%   fileList  - the filelist as char array or cellstr array (3D or 4D images as .nii or .hdr)
%
% Returns:
%   dataset   - the datset with included data4D and data_3DNiftiHdr - struct
%
% Comments:
%
function [dataset] = setDataset_data_ByFilelist(dataset, fileList)

  if( ~exist('dataset','var') || ~exist('fileList','var') ) 
    error('Usage of setDataset_data_ByFilelist: [dataset] = setDataset_data_ByFilelist(dataset, fileList)');
  end
  
  if(isempty(fileList))
    error('Given filelist is empty, please check your file-selection!');
  end
  
  %check which type the filelist is
  %2 types are supported: 1. char array and 2. cellstr array
  if(ischar(fileList))
    fileList = cellstr(fileList);
  end
  
  if(~iscellstr(fileList))
    error('Check input type for filelist! Supported types are: character array and cellstr array!');
  end
  
  dataset.data = [];
  
   
  if(~easyupMVPA_getGlobals('quietMode'))
    disp('INFO: setDataset_data_ByFilelist: Loading image files ...');
  end
  
  
  %Load the first file and check for 3D or 4D
  %It is expected that if the first is 3D all the rest is 3D too and vice versa
  fileString = fileList{1};
  %check existance of the file
  if(~exist(fileString, 'file') )
    error(['Could not read image file: ', fileString]);
  end
  if(~easyupMVPA_getGlobals('quietMode'))
        disp(['Loading file: ', fileString]);
  end
  %Read the data (nii or hdr)
  %the nii loading function is slightly changed for better memory performance
  dataNii = load_untouch_nii(fileString);
  sizeIMG = size(dataNii.img);
  
  
  %CASE 3D
  if(length(sizeIMG) == 3)
    
    %use preallocation, because of speed and memory performance
    dataset.data = zeros(sizeIMG(1), sizeIMG(2), sizeIMG(3), size(fileList,1),'int16');
    dataset.data(:,:,:,1)   = dataNii.img;
    dataset.data_3DNiftiHdr = dataNii.hdr;
    
    %load all data and store in data
    for i = 2:size(fileList,1)
      fileString = fileList{i};
      if(~easyupMVPA_getGlobals('quietMode'))
        disp(['Loading file: ', fileString]);
      end

      %check existance of the file
      if(~exist(fileString, 'file') )
        error(['Could not read image file: ', fileString]);
      end

      %Read the data (nii or hdr)
      %the nii loading function is slightly changed for better memory performance
      dataNii = load_untouch_nii(fileString);
      dataset.data(:,:,:,i)   = dataNii.img;
      
    end%endfor
    
  %CASE 4D  
  elseif(length(sizeIMG) == 4)
    
    %set the data of the image loaded already
    dataset.data            = dataNii.img;
    dataset.data_3DNiftiHdr = dataNii.hdr;
       
    %Load all the data in the filelist and concatenate in dataset.data
    for i=2:size(fileList,1)
      fileString = fileList{i};

      if(~easyupMVPA_getGlobals('quietMode'))
        disp(['Loading file: ', fileString]);
      end

      %check existance of the file
      if(~exist(fileString, 'file') )
        error(['Could not read image file: ', fileString]);
      end

      %Read the data (nii or hdr)
      %the nii loading function is slightly changed for better memory performance
      dataNii = load_untouch_nii(fileString);
      dataset.data = cat(4, dataset.data, dataNii.img);
    
    end %endfor
    
  else
      error('Unrecognized dimensions for image data (3D or 4D expected)!');  
  end
    
    
  %check 4th dimension of given data
  if(size(dataset.data, 4) < 1)
    disp('ERROR: setDataset_data_ByFilelist: Please check the dimensions! The , 4th dimension is bel 1 ! This dimension should be your number of samples (timesteps, scans)!)');
    dataset.data = [];
    return;
  end
  
  dataset.dataFilelist = fileList;
  
end