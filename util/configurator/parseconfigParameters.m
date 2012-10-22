

function myDataset = parseconfigParameters(configParameters);

%-- start matlabpool
% if exist('cpu_cores')
%     if matlabpool('size') == 0
%         matlabpool('open', cpu_cores)
%     end
% end


for subj = 1:numel(configParameters.subjectlist);
    clearvars -except configParameters subj
    configParameters.subjectname = configParameters.subjectlist{subj};
    if(~easyupMVPA_getGlobals('quietMode'))
        disp('')
        disp('')
        disp('===============================================');
        disp(['STARTING MVPA PROCEDURE FOR SUBJECT ',configParameters.subjectname]);
        disp('===============================================');
    end

    %--load default parameters into configParameters STRUCT
    configParameters = getDefaultConfigParameters(configParameters);

    
    %--determine the current analysis directory
    configParameters.currentAnalysisDirectory = fullfile(configParameters.baseDirectory,configParameters.subjectname,configParameters.dir_analysis_mvpa);
    
    
    %--determine the filename given the preprocessing options
    configParameters.DataSetName = getDataSetname(configParameters);
    
    
    %--load the .nii data OR .mat data (depending whether the correct ones
    %do exist. preprocessing takes place here!
    myDataset = loadDataFromMatOrNii(configParameters);
    
    
    %-- override eventually OLD configParameters, which were loaded, with new configParameters
    myDataset.configParameters = configParameters;
    
    
    %--get the info .txt file
    myDataset = getInfoFile(myDataset);
    
    
    %-- determine the analysis SUBFOLDER 
    
    myDataset.configParameters.currentAnalysisSubfolder=fullfile(myDataset.configParameters.currentAnalysisDirectory,myDataset.configParameters.analysisSubfolder);
   

    %--set a brain mask. check whether one should be created and if it is
    %already in place. if so, use it, otherwise create it.
    myDataset = setBrainMaskbyconfigParameters(myDataset);
    
    
    %--get the chunks and classIDs, dont do that for SPM beta files
    if myDataset.configParameters.useSpmBetaFiles %if SPM beta files are selected, a splitscheme is generated automatically
    else
        myDataset = setDataset_chunks_classID_byconfigParameters(myDataset);
    end
    
    
    
    %--save preprocessing to .mat file (create folder if it doesn't exist
    if configParameters.savePreProcessing;
        if exist(configParameters.currentAnalysisDirectory)~=7 mkdir(configParameters.currentAnalysisDirectory); end
        fileString =fullfile(configParameters.currentAnalysisDirectory,configParameters.DataSetName);
        if(~easyupMVPA_getGlobals('quietMode'))
            disp(['INFO: Saving file: ', fileString]);
        end
        tic
        save(fileString,'myDataset','-v7.3');
        time_save = toc;
        if(~easyupMVPA_getGlobals('quietMode'))
            disp(['INFO: File saved in ', num2str(round(time_save)),' ','seconds.']);
        end
    end
    
    %-- AVERAGING (using what is specified in averagingSplitMethod
    if configParameters.averaging
        myDataset = getSplitScheme(myDataset,'AVERAGING');
        myDataset = averageOverSplitScheme(myDataset);
    end
    
    
    %--GET AN LOOCV SCHEME? ONLY FOR ROI AND 
    if strcmp(configParameters.BaseMethod,'sl') || strcmp(configParameters.BaseMethod,'SL') ||  strcmp(configParameters.BaseMethod,'roi')||  strcmp(configParameters.BaseMethod,'ROI')||  strcmp(configParameters.BaseMethod,'cella')
        if myDataset.configParameters.useSpmBetaFiles %if SPM beta files are selected, a splitscheme is generated automatically
        else
            myDataset = getSplitScheme(myDataset,'LOOCV');
        end
    end
   
    
    
    
    %-- CHECK ALL POSSIBLE BASEMETHODS AND COMPUTE THEM
    %-- OPTION BASIC SEARCHLIGHT
        if strcmp(configParameters.BaseMethod,'sl') || strcmp(configParameters.BaseMethod,'SL')
            if(~easyupMVPA_getGlobals('quietMode'))
                disp(['INFO: Starting standard searchlight procedure with a diameter of ',num2str(myDataset.configParameters.SearchLightDiameter),'voxels']);
            end
            doSearchlightBatch(myDataset);
        end

        %-- OPTION ROI DECODING
        if strcmp(configParameters.BaseMethod,'roi') || strcmp(configParameters.BaseMethod,'ROI')
            disp(['INFO: Startin ROI based decoding']);
            doROIdecodingBatch(myDataset);
        end
    
    
    
end % of subject loop

%-- SECOND LEVEL STUFF
%-- BASIC SEARCHLIGHT
if strcmp(configParameters.BaseMethod,'sl') || strcmp(configParameters.BaseMethod,'SL')
    doSearchlightGroupMap(myDataset)
end

%-- ROI DECODING
if strcmp(configParameters.BaseMethod,'roi') || strcmp(configParameters.BaseMethod,'ROI')
    
    doROIdecodingGroupStats(myDataset)
end



%--close matlabpool
if exist('configParameters.cpu_cores')
    if matlabpool('size') ~= 0
        matlabpool close
    end
end