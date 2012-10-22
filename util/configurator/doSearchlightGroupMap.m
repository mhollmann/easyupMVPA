% A very simple second level (group level) function, 
% which averages over the individual accuracy maps to create a group
% accuracy map. 
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
% SPM is needed for this function! The individual accuracy maps are
% smoothed (inter-subject variance...) using SPM and a FWHM specified in
% secondLevelSmoothingFWHM (by default 4mm). The resulting group accuracy
% maps is then saved in baseDirectory / _MVPA2nd.
% All specifications are taken from the dataset. Alternatively this
% function may be run without input, in that case any of the configParameters.mat 
% files will work (which are found in each subjects directory / analysis_mvpa / TIME )

%   [myDataset] = doSearchlightGroupMap(myDataset)

 
function [] = doSearchlightGroupMap(myDataset)

if nargin == 0
    tmpdir = uigetdir('Select subfolder in _MVPA2nd');
    load(fullfile(tmpdir,'configParameters.mat'));
    myDataset.configParameters = configParameters;
    folderextension = input('Folder-Extension? ');
    subjectsel = input('Subject-Selection? ');
    
    if numel(folderextension) > 0
        mvpa_2ndleveldir = [tmpdir,'-',folderextension];
        mkdir(mvpa_2ndleveldir);
    end
    
    if numel(subjectsel) > 0
        myDataset.configParameters.subjectlist = subjectsel;
    end
    
     mvpa_2ndleveldir = fullfile(myDataset.configParameters.baseDirectory,'_MVPA2nd',myDataset.configParameters.analysisSubfolder);
     mkdir(mvpa_2ndleveldir);

else
    mvpa_2ndleveldir = fullfile(myDataset.configParameters.baseDirectory,'_MVPA2nd',myDataset.configParameters.analysisSubfolder);
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
        
        tmpdir = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectlist{subj},myDataset.configParameters.dir_analysis_mvpa,myDataset.configParameters.analysisSubfolder);
        
        if myDataset.configParameters.secondLevelSmoothingFWHM > 0
            tmpname_unsmoothed = [myDataset.configParameters.accuracymapsavename,'_',condition1_name,'-',condition2_name,'_',myDataset.configParameters.subjectlist{subj},'.nii'];
            tmpfilename = ['s', tmpname_unsmoothed];
            tmpname = fullfile(tmpdir,['s', tmpname_unsmoothed]);
            tmpname_unsmoothed = fullfile(tmpdir,tmpname_unsmoothed);
            
            spm_smooth(tmpname_unsmoothed,tmpname,myDataset.configParameters.secondLevelSmoothingFWHM);
        else
            tmpname = [myDataset.configParameters.accuracymapsavename,'_',condition1_name,'-',condition2_name,'_',myDataset.configParameters.subjectlist{subj},'.nii'];
            tmpfilename = tmpname;
            tmpname = fullfile(tmpdir,tmpname);
        end
        
        copyfile(fullfile(tmpdir,tmpfilename),fullfile(firstlevelcollectiondir,tmpfilename));
        
        
        tmp = load_nii(tmpname);
        clear configParameters
        load(fullfile(tmpdir,'configParameters.mat'));
        
        
        switch subj
            case 1
                groupmap = zeros(size(tmp.img,1),size(tmp.img,2),size(tmp.img,3),numel(myDataset.configParameters.subjectlist));
        end
        groupmap(:,:,:,subj) = tmp.img;
    end
    
    groupmean = mean(groupmap,4);
    tmp.img = groupmean;
    tmp.hdr.dime.datatype = 16;
    tmp.hdr.dime.bitpix = 32;
    groupsavename = [myDataset.configParameters.accuracymapsavename,'_',condition1_name,'-',condition2_name,'_','groupaverage','.nii'];
    save_nii(tmp,fullfile(mvpa_2ndleveldir,groupsavename));

    %tmp.hdr.dime.dim(5) = 
    


    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Group Accuracy Map saved']);
        disp(fullfile(mvpa_2ndleveldir,groupsavename));
    end
    
    

        

end

%save configParameters
configParameters=myDataset.configParameters;
save(fullfile(mvpa_2ndleveldir,'configParameters.mat'),'configParameters');
dlmwrite(fullfile(mvpa_2ndleveldir,'InfoText.txt'),myDataset.configParameters.InfoTxt,'')
