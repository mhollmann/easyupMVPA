% Takes all accuracymaps which were created and averages them...
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [myDataset] = doSearchlightGroupMap(myDataset)

 
function [] = doROIdecodingGroupStats(myDataset)

if nargin == 0
    tmpdir = uigetdir('Select subfolder in _MVPA2nd');
    load(fullfile(tmpdir,'configParameters.mat'));
    myDataset.configParameters = configParameters;
    folderextension = input('Folder-Extension? ');
    subjectsel = input('Subject-Selection? ');
    
    if numel(folderextension) > 0
        mvpa_2ndleveldir = [mvpa_2ndleveldir,'-',folderextension];
        mkdir(mvpa_2ndleveldir);
    end
    
    if numel(subjectsel) > 0
        myDataset.configParameters.subjectlist = subjectsel;
    end
    
    mvpa_2ndleveldir = fullfile(myDataset.configParameters.baseDirectory,'_MVPA2nd',myDataset.configParameters.timeString);
    mkdir(mvpa_2ndleveldir);

else
    mvpa_2ndleveldir = fullfile(myDataset.configParameters.baseDirectory,'_MVPA2nd',myDataset.configParameters.timeString);
    mkdir(mvpa_2ndleveldir);
 
end

firstlevelcollectiondir = fullfile(mvpa_2ndleveldir,'firstlevel');
mkdir(firstlevelcollectiondir);



for c = 1:numel(myDataset.configParameters.comparisons)
    tmp = num2str(myDataset.configParameters.comparisons(c));
    condition1 = str2double(tmp(1));
    condition2 = str2double(tmp(2));
    condition1_name = myDataset.configParameters.conditions{condition1};
    condition2_name = myDataset.configParameters.conditions{condition2};

    for subj=1:numel(myDataset.configParameters.subjectlist)
        
        tmpdir = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectlist{subj},myDataset.configParameters.dir_analysis_mvpa,myDataset.configParameters.timeString);
        
        savename = ['ROIaccuracy','_',condition1_name ,'-',condition2_name ,'_',myDataset.configParameters.subjectlist{subj},'.txt'];
        groupacc(subj) = load(fullfile(tmpdir,savename));
    end
    
    groupacc_mean = mean(groupacc)
    groupacc_std = std(groupacc)

    groupsavename = ['ROIaccuracy','_',condition1_name,'-',condition2_name,'_','groupaverage','.txt'];
    %dlmwrite(fullfile(mvpa_2ndleveldir,groupsavename),groupacc,',')
    csvwrite(fullfile(mvpa_2ndleveldir,groupsavename),groupacc)
    %tmp.hdr.dime.dim(5) = 
    


    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Group Accuracy Values saved']);
        disp(fullfile(mvpa_2ndleveldir,groupsavename));
    end
    
    

        

end

%save configParameters
configParameters=myDataset.configParameters;
save(fullfile(mvpa_2ndleveldir,'configParameters.mat'),'configParameters');
dlmwrite(fullfile(mvpa_2ndleveldir,'InfoText.txt'),myDataset.configParameters.InfoTxt,'')
