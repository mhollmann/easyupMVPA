% Set the chunks of a dataset by a given attribute-file.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = setDataset_chunks_ByAttribFile(dataset, attribFile)
%
%   This methods sets the chunks of the given dataset. Chunks can be used to
%   subdivide a dataset. One may for example set a separate chunk for all transition 
%   scans for easy removing after preprocessing. 
%   The given attributes file must be in ASCII-format.
%   An attribute file has the content (classID  chunk):
%
%     0       0
%     0       1
%     0       1
%     1       0
%     1       2
%     1       2
%     .       .
%     .       .
%
%
% Parameters:
%   dataset     - the datset to set the data4D for
%   attribFile  - ASCII-formatted file holding chunk and classID information
%
% Returns:
%   dataset   - the datset with included chunks
%
% Comments:
%
function [dataset] = setDataset_chunks_ByAttribFile(dataset, attribFile)

  if( ~exist('dataset','var') || ~exist('attribFile','var') ) 
    error('Usage of setDataset_chunks_ByAttribFile: [dataset] = setDataset_chunks_ByAttribFile(dataset, attribFile [ASCII-formatted file])');
  end

  if(~exist(attribFile, 'file'))
    error(['setDataset_chunks_ByAttribFile: Attributes file: ', attribFile, ' does not exist!']);
  end
  
  fid = fopen(attribFile);
  attribs = textscan(fid, '%s', 'delimiter', '\n');
  fclose(fid); 
  
  dataset.chunks = zeros(1, size(attribs{1},1), 'uint16');
  
  %loop over lines
  for i=1:size(attribs{1},1)
    attribLine  = str2num(char(attribs{1}(i)));
    
    if(numel(attribLine)<2)
      error('ERROR: It seems as there are empty or incomplete lines in your attribute file, please fix that first!');
    end

    dataset.chunks(i) = attribLine(2);
  end  
  
  
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_chunks_ByAttribFile: In the current state the dataset is not suitable for further processing, please see messages before!');
  end
  
end