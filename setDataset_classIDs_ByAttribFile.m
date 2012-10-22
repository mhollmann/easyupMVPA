% Set the class-IDs of a dataset by a given attribute-file.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset] = setDataset_classIDs_ByAttribFile(dataset, attribFile)
%
%   This methods sets the classIDs of the given dataset. 
%   The given attributes file must be in ASCII-format.
%   An attributse file has the content (classID  chunk):
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
%   dataset     - the datset to set the classIDs for
%   attribFile  - ASCII-formatted file holding chunk and classID information
%
% Returns:
%   dataset   - the datset with included classIDs
%
% Comments:
%
function [dataset] = setDataset_classIDs_ByAttribFile(dataset, attribFile)

  if(~exist('dataset', 'var') || ~exist('attribFile', 'var'))
    disp('Usage of setDataset_classIDs_ByAttribFile : [dataset] = setDataset_classIDs_ByAttribFile(dataset, attribFile [ASCII-formatted file])');
    return;
  end


  if(~exist(attribFile, 'file'))
    error(['setDataset_classIDs_ByAttribFile: Attributes file: ', attribFile, ' does not exist!']);
  end
  
  fid = fopen(attribFile);
  attribs = textscan(fid, '%s', 'delimiter', '\n');
  fclose(fid); 
  
  dataset.classIDs = zeros(1, size(attribs{1},1), 'uint8');
  
  %loop over lines
  for i=1:size(attribs{1},1)
    attribLine  = str2num(char(attribs{1}(i)));
    dataset.classIDs(i) = attribLine(1);
  end  
  
  if(~checkDataset(dataset))
    disp('WARNING: setDataset_classIDs_ByAttribFile: In the current state the dataset is not suitable for further processing, please see messages before!');
  end
    
end