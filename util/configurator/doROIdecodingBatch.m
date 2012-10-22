function [myDataset] = doROIdecodingBatch(myDataset)

mkdir(myDataset.configParameters.currentAnalysisSubfolder);

for c = 1:numel(myDataset.configParameters.comparisons)
    tmp = num2str(myDataset.configParameters.comparisons(c));
    class1 = str2double(tmp(1));
    class2 = str2double(tmp(2));

    class1_name = myDataset.configParameters.conditions{class1};
    class2_name = myDataset.configParameters.conditions{class2};
    
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: ROI decoding conditions ',class1_name,' vs. ',class2_name]);
    end
    
    tic
    accuracy = doROIdecoding(myDataset,class1,class2);
    
    time_roidecoding = toc;
    savename = ['ROIaccuracy','_',class1_name,'-',class2_name,'_',myDataset.configParameters.subjectname,'.txt'];

    %save comments as filename of empty txt file
    dlmwrite(fullfile(myDataset.configParameters.currentAnalysisSubfolder,savename),accuracy,'')
    
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: ROI decoding procedure completed in ',num2str(time_roidecoding),' seconds.']);
        disp(['INFO: Accuracy for this ROI (and participant) is ',num2str(accuracy),'% over chance. File saved:']);
        disp(fullfile(myDataset.configParameters.currentAnalysisSubfolder,savename));
    end

end



%get configParameters
configParameters=myDataset.configParameters;

%save as MAT
save(fullfile(myDataset.configParameters.currentAnalysisSubfolder,'configParameters.mat'),'configParameters');

%save as TXT
dlmwrite(fullfile(myDataset.configParameters.currentAnalysisSubfolder,'InfoText.txt'),myDataset.configParameters.InfoTxt,'')

%save comments as filename of empty txt file
dlmwrite(fullfile(myDataset.configParameters.currentAnalysisSubfolder,['_',configParameters.Comments,'.txt']),'','')

