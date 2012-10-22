% This script explicitly writes down .nii files for the preprocessing. All
% settings are analogue to the easyupMVPA configurator and have to be
% specified.
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:

% Comments:
clear all

% -- BASICS -----------

baseDirectory = '/scr/magnesium1/7T_PPI';
subjectlist = mpc_getsubjects(baseDirectory,4);%{'LH4T','MT3T','NM6T'};%
conditions = {'p1','p2','i1','i2'};

cpu_cores = 1;
global propertyStruct

propertyStruct.quietMode = 0;

% --  SCANS -----------
dataFileFormat = 'rdata.nii';
spm_betaformat = 0;
fruns(1).directory = 'functional_onenii';
% fruns(1).directory = 'functional_c/s1';
% fruns(2).directory = 'functional_c/s2';
% fruns(3).directory = 'functional_c/s3';
% fruns(4).directory = 'functional_c/s4';
% fruns(5).directory = 'functional_c/s5';
% fruns(6).directory = 'functional_c/s6';
TR = 3;
highpassFilter = 90;
detrending = 0;








timeVector = clock;

% --  DON'T EDIT BELOW -----------

%-- loop over all subjects
for subj = 1:numel(subjectlist);
    clear myDataset
    subjectname = subjectlist{subj}
    clear myDataset niifile
    prefix = '';
    for j = 1:numel(fruns)
        runs(1).directory = fruns(j).directory;
        %-- collect all specified variables and put them in a structure
        configSpecifiedVariables=who;
        for i=1:numel(configSpecifiedVariables)
            if strcmp(configSpecifiedVariables{i},'configParameters') && strcmp(configSpecifiedVariables{i},'configSpecifiedVariables')
            else
                eval(['configParameters.',configSpecifiedVariables{i},'=',configSpecifiedVariables{i},';']);
            end
        end
        
        
       
        %--load default parameters
        configParameters = getDefaultConfigParameters(configParameters);
        
        
        
        
        %--load the raw data CHECK FOR EXISTING FILES!! AND SAVE STUFF IF NOT
        myDataset = getFilledDataset(configParameters);
        if detrending
            [myDataset] = doLinearDetrending(myDataset)
            prefix = ['d',prefix];
        end
        
        
        myDataset = doHighpassFiltering(myDataset, 1/TR, 1/highpassFilter);
        myDataset = doZScoring(myDataset);
        prefix = ['zf',num2str(highpassFilter),prefix];
        niifile = load_untouch_nii(fullfile(baseDirectory,subjectname,fruns(j).directory,dataFileFormat));
        niifile.img = myDataset.data;
        niifile.hdr.dime.datatype = 16;
        niifile.hdr.dime.bitpix = 32;
        save_untouch_nii(niifile,fullfile(baseDirectory,subjectname,fruns(j).directory,[prefix,dataFileFormat]));
        
    end
    
end % of subject loop



