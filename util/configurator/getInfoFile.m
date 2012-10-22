% Returns info txt file to be saved later
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [myDataset] = getInfoFile(myDataset)
%


function [myDataset] = getInfoFile(myDataset)

infotxt = char();

tmp = ['Comments: ',myDataset.configParameters.Comments];
infotxt=char(infotxt,tmp);

tmp = ['Basemethod: ',myDataset.configParameters.BaseMethod];
infotxt=char(infotxt,tmp);

tmp = ['averaging: ',num2str(myDataset.configParameters.averaging)];
infotxt=char(infotxt,tmp);

tmp = ['loocvSplitMethod: ',myDataset.configParameters.loocvSplitMethod];
infotxt=char(infotxt,tmp);

tmp = ['removeTailTransitions: ',num2str(myDataset.configParameters.removeTailTransitions)];
infotxt=char(infotxt,tmp);

tmp = ['removeTrailingTransitions: ',num2str(myDataset.configParameters.removeTrailingTransitions)];
infotxt=char(infotxt,tmp);

tmp = ['baseDirectory: ',myDataset.configParameters.baseDirectory];
infotxt=char(infotxt,tmp);

tmp = ['comparisons: ',num2str(myDataset.configParameters.comparisons)];
infotxt=char(infotxt,tmp);

tmp = ['conditions: ',cell2mat(myDataset.configParameters.conditions)];
infotxt=char(infotxt,tmp);

tmp = ['secondLevelSmoothingFWHM: ',num2str(myDataset.configParameters.secondLevelSmoothingFWHM)];
infotxt=char(infotxt,tmp);

if isfield(myDataset.configParameters,'obtainmask')
    tmp = ['mask was obtained automatically'];
else
    tmp = ['mask: ',num2str(myDataset.configParameters.mask)];
end

infotxt=char(infotxt,tmp);

tmp = ['dataFileFormat: ',num2str(myDataset.configParameters.dataFileFormat)];
infotxt=char(infotxt,tmp);



myDataset.configParameters.InfoTxt = infotxt;







