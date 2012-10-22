% Set the data field of a dataset by a given myDataset.configParameters variable.
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:
%
%   [dataset] = setDataset_data_ByconfigParameters(myDataset,configParameters)

%   This methods sets the data (timeseries of 3D images) for the dataset.
%   The input variable myDataset.configParameters is determined by the
%   easyupMVPA_configurator.

%
% Returns:
%   dataset   - the datset with included data4D and data_3DNiftiHdr - struct
%
% Comments:
%
function [myDataset] = setDataset_data_byconfigParameters(myDataset)


for r=1:numel(myDataset.configParameters.runs)
    fileListEntry = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.runs(r).directory,myDataset.configParameters.dataFileFormat);
    fileList(r,1:numel(fileListEntry)) = fileListEntry;
end

[myDataset] = setDataset_data_ByFilelist(myDataset, fileList);