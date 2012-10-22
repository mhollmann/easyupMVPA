% Returns the standardized filename for a dataset, given the preprocessing
% parameters specified in the configurator file
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [DataSetName] = getDataSetName(configParameters)



function [DataSetName] = getDataSetname(configParameters)

% strip off the filname extension .nii .hdr .img
DataSetName = configParameters.dataFileFormat;
if strfind(DataSetName,'.nii');
    DataSetName(strfind(DataSetName,'.nii'):end) = [];
end
if strfind(DataSetName,'.img');
    DataSetName(strfind(DataSetName,'.img'):end) = [];
end
if strfind(DataSetName,'.hdr');
    DataSetName(strfind(DataSetName,'.hdr'):end) = [];
end

% put a _ in front
DataSetName = strcat('_',DataSetName);

% check the preprocessing and add letters
if configParameters.detrending
    DataSetName = strcat('d',DataSetName); 
end
if configParameters.highpassFilter
   DataSetName = strcat('h',DataSetName);
end
if configParameters.z_scoring
   DataSetName = strcat('z',DataSetName);
end

DataSetName = strcat(configParameters.dataset_name,DataSetName,'.mat');








