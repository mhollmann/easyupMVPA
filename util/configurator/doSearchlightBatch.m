% Runs a searchlight procedure over all comparisons specified in the variable comparisons. 
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
% The resulting accuracy maps are saved in the mvpa analysis folder of each subject with the naming
% convention: accmap_comparison_subject (if you compare class1 vs
% class2 on subject ABCD the accuracy map would be named
% accmap_class1-class2_ABCD.nii. The values of this accuracymap are
% the percent of correct decoding against chance, i.e. all possible values
% lie between -50 and +50. This means a value of some location of '10'
% corresponds to 10% over chance classification. The .nii files are saved
% in a subfolder of the mvpa_analysis directory which is named by the DATE
% of the script execution.
%
%   [myDataset] = doSearchlightBatch(myDataset,subjectname)

 
function [] = doSearchlightBatch(myDataset)


mkdir(myDataset.configParameters.currentAnalysisSubfolder);

for c = 1:numel(myDataset.configParameters.comparisons)
    tmp = num2str(myDataset.configParameters.comparisons(c));
    class1 = str2double(tmp(1));
    class2 = str2double(tmp(2));
    class1_name = myDataset.configParameters.conditions{class1};
    class2_name = myDataset.configParameters.conditions{class2};
    
    tic
    
    clear accuracymap niifile
    accuracymap = doSearchlight_p(myDataset,class1,class2);

    time_searchlight = toc;
    savename = [myDataset.configParameters.accuracymapsavename,'_',class1_name,'-',class2_name,'_',myDataset.configParameters.subjectname,'.nii'];
    
    niifile = make_nii(accuracymap);
    niifile.hdr = myDataset.data_3DNiftiHdr;
    niifile.hdr.dime.dim(5) = 1;
    niifile.hdr.dime.datatype = 16;
    niifile.hdr.dime.bitpix = 32;
    
    %save nii
    save_nii(niifile,fullfile(myDataset.configParameters.currentAnalysisSubfolder,savename));
    
    
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Seachlight procedure completed in ',num2str(time_searchlight),' seconds. File saved:']);
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

clearvars -except myDataset

