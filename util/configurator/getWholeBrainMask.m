% Automatically creates a whole brain roi for voxels that a) exceed a mean
% certain signal intensity and b) a specific std SIMULTANEOUSLY
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%


function [configParameters] = getWholeBrainMask(configParameters);

thresh_mean = 500;
thresh_std = 1;


dir_save = fullfile(configParameters.baseDirectory,configParameters.subjectname,configParameters.dir_analysis_mvpa);
if exist(dir_save)~=7 mkdir(dir_save); end
functionalfile = fullfile(configParameters.baseDirectory,configParameters.subjectname,configParameters.runs(1).directory,configParameters.dataFileFormat);
niifile = load_nii(functionalfile);
data = niifile.img;
data = single(data);
datadim = size(data);
exc_mean = find(mean(data,4) > thresh_mean);
exc_std = find(std(data,1,4) > thresh_std);
exc_both = exc_mean(ismember(exc_mean,exc_std));
mask = zeros(datadim(1),datadim(2),datadim(3));
mask(exc_both) = 1;
niifile.fileprefix = fullfile(dir_save,'wholebrainmask.nii');
niifile.hdr.dime.dim = [3,datadim(1),datadim(2),datadim(3),1,1,1,1];
niifile.img = mask;
save_nii(niifile,fullfile(dir_save,'wholebrainmask.nii'));
%view_nii(niifile)











