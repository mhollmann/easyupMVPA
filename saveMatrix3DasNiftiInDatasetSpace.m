% Save a 3D matrix as nifti-image, according to nifti-hdr info in dataset.
% 
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
% saveMatrix3DasNiftiInDatasetSpace(dataset, matrix3D, hdrFilename, [scaleFactor])
%
% This methods saves a 3D matrix as Nifti-hdr/img in the 3D space of a single volume of the dataset with type dataset4D. 
% The written file will be in Analyze-Nifti HDR/IMG Format. If file already exists it will be overwritten.
%
% This function can just be used with datasets of type: dataset4D.
%
% Parameters:
%   dataset     - the datset that defines the space (by field data_3DNiftiHdr)
%   matrix3D    - the 3D matrix to save
%   hdrFilename - filename of the header
%   scaleFactor - [optional] a scale factor multiplied with dataset (this scaling is directly applied and NOT written in header)
%   dataType    - [optional] one of the char 'int16'(default) or 'float32'
%
% Returns:
%
% Comments:
%
%
function saveMatrix3DasNiftiInDatasetSpace(dataset, matrix3D, hdrFilename, scaleFactor, dataType)

  if( ~exist('dataset','var') || ~exist('matrix3D','var') || ~exist('hdrFilename','var')) 
    error('Usage of saveMatrix3DasNiftiInDatasetSpace: saveMatrix3DasNiftiInDatasetSpace(dataset, matrix3D (3D matrix), hdrFilename (filename of the header), scaleFactor (optional)  )');
  end

  if( ~dataset.is4D)
     error('Function saveMatrix3DasNiftiInDatasetSpace can just be used for 4D datasets! ');
  end

  
  if(exist('scaleFactor','var'))
    %this can mean the dataType is given
    if(ischar(scaleFactor))
      dataType = scaleFactor;
    end
  end
  
  if(exist('dataType', 'var'))
      if(strcmp('scaleFactor', 'int16') || strcmp('scaleFactor', 'float32')) 
        
      end
    
  end
  
  
  %scale if necessary
  if(exist('scaleFactor','var'))
    matrix3D = matrix3D*scaleFactor;
  end
    
  niiHdr = dataset.data_3DNiftiHdr;

  %set the 4th dimension to 1 because it is just a 3D array
  niiHdr.dime.dim(5) = 1;
  
  [pathstr, name, ext, versn] = fileparts(hdrFilename);
  
  if( ~strcmp(ext,'.hdr') )
    disp('INFO: saveMatrix3DasNiftiInDatasetSpace: Extension of header filename was set to: .hdr !' );
    ext = '.hdr';
  end
  fileNameHDR = [pathstr, filesep(),name, ext];
  
  
  if(~exist(pathstr,'dir'))
    error('The given PATH in HDR-filename is not existant. Hdr-filename is not valid.');
  end
  
  %save the hdr data 
  fid_writefile_4d_hdr = fopen(fileNameHDR, 'w');
  
  %if field magic does not exist the data is old analyze-style
  if(isfield(niiHdr.hist, 'magic'))
    niiHdr.hist.magic = 'ni1';
    save_untouch_nii_hdr(niiHdr, fid_writefile_4d_hdr);
  else
    save_untouch0_nii_hdr(niiHdr, fid_writefile_4d_hdr);
  end
    
  fclose(fid_writefile_4d_hdr);

  %save the img data
  fid_writefile_4d_img = fopen([pathstr, filesep(), name, '.img'], 'w');
  fwrite(fid_writefile_4d_img, matrix3D, 'int16');
  fclose(fid_writefile_4d_img);
  

end