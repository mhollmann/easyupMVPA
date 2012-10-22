% Returns a filled dataset given the input of the configurator
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [dataset] = getFilledDataset()

%   Returned sructure:
%
%   dataset.type                  = 'dataset4D';
%   dataset.is4D = true;
%   dataset.is2D = false;
%   dataset.dataFilelist          = files specified;
%   dataset.data                  = ...;
%   dataset.chunks                = ...;
%   dataset.classIDs              = ...;
%   dataset.mask                  = ...;
%   dataset.featureSelectionMap   = ...;
%   dataset.data_3DNiftiHdr       = ...;
%   dataset.processingHistory     = ...;
%
% Parameters:
%
% Returns:
%
%   dataset          -> filled dataset
%
% Comments:
%



function [myDataset] = getFilledDataset(configParameters)
%get empty dataset first
[myDataset] = getEmpty4DDataset();

%put the configParameters inside
myDataset.configParameters = configParameters;

%either load in the .nii files OR load in spm beta files
switch myDataset.configParameters.useSpmBetaFiles
    case 0
        myDataset = setDataset_data_byconfigParameters(myDataset);
    case 1
        myDataset = setDataset_data_bySPMbetavalues(myDataset);
end











