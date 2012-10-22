% dddddddddddddddddddd
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [myDataset] = loadDataFromMatOrNii(configParameters)




function [myDataset] = loadDataFromMatOrNii(configParameters)


loadnii = 1;

currentAnalysisDirectory = fullfile(configParameters.baseDirectory,configParameters.subjectname,configParameters.dir_analysis_mvpa);


%check whether .mat file exists at correct place
fileString = fullfile(currentAnalysisDirectory,configParameters.DataSetName);
if exist(fileString)
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Loading dataset MAT file: ', fileString]);
    end
    tic
    tmp=load(fileString);
    time_load = toc;
    %if file existed - is the number of runs the same?
    if numel(tmp.myDataset.dataFilelist) == numel(configParameters.runs)
        loadnii = 0;
        %if number of runs the same, are the nii filenames the same?
        for k=1:numel(tmp.myDataset.dataFilelist)
            if strcmp(tmp.myDataset.dataFilelist{k},fullfile(configParameters.baseDirectory,configParameters.subjectname,configParameters.runs(k).directory,configParameters.dataFileFormat))
            else loadnii = 1;
                if(~easyupMVPA_getGlobals('quietMode'))
                    disp('INFO: Loaded dataset MAT does NOT correspond to your settings. Loading .nii files');
                end
            end
        end
    end 
end


%load in files from nii (load the image files) or mat (load an existing
%easyupMVPA dataset)

if loadnii
    myDataset = getFilledDataset(configParameters);
    
else
    myDataset = tmp.myDataset;
    clear tmp
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Loaded dataset MAT file corresponds to your settings. Dataset loaded in ',num2str(round(time_load)),'seconds.']);
    end
end


%do preprocessing if data was loaded
if loadnii
    if configParameters.detrending; myDataset = doLinearDetrending(myDataset); end
    if configParameters.highpassFilter; myDataset = doHighpassFiltering(myDataset, 1/configParameters.TR, 1/configParameters.highpassFilter); end
    if configParameters.z_scoring; myDataset = doZScoring(myDataset); end
end
   