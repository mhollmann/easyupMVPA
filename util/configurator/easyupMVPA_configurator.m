% ----------- CONFIGURATOR -----------
clear all



% -- BASICS -----------
Comments = 'FULL BRAIN, LAST SHOT, 7VX';
baseDirectory = '/scr/magnesium3/7T_PPI';
subjectlist = mpc_getsubjects(baseDirectory,4); 
%{'BAIT','EM4T','LH4T','NM6T','NS6T','OD2T','SAST'}
%{'MT3T','PF5T','QC1T','RSDT','SL7T','SMWT','WC8T'};
conditions = {'p1','p2','i1','i2'}; 
comparisons = [12];


cpu_cores = 8;



% --  SCANS -----------

dataFileFormat = 'zfwrdata.nii'; 
spm_betaformat = 0; 
runs(1).directory = 'functional_c/s1'; 
runs(2).directory = 'functional_c/s2';
runs(3).directory = 'functional_c/s3';
runs(4).directory = 'functional_c/s4';
runs(5).directory = 'functional_c/s5';
runs(6).directory = 'functional_c/s6';
TR = 3; %...
mask = 'functional_c/rgreymatterroi.nii'; 

% --  PREPROCESSING -----------
detrending = 0; 
highpassFilter = 0; 
z_scoring = 0; 
savePreProcessing = 0;

% --  ONSETS, DURATIONS, DATA GROUPING -----------
useSpmBetaFiles = 0; 
sameTrialDurations = 4;
removeTrailingTransitions =3; 
removeTailTransitions = 0;
unitofdesign = 'scans'; 
spm_onsets_specification = 1;
averaging = 0; 
averagingSplitMethod = 'chunks'; 



% -- CLASSIFICATION -----------
%BaseMethod = 'SL';
%BaseMethod = 'SL-swap';

BaseMethod = 'SL';
%BaseMethod = 'parcellation';
%BaseMethod = 'RFE';
SearchLightDiameter = 7;
loocvSplitMethod = 'runs';


% -- GROUP STATS
secondLevelSmoothingFWHM = 6; %nach der classification wird fuer die gruppenstat gesmootht.
accuracymapsavename = 'accmap';




% --  ---------------- -----------
% --  DON'T EDIT BELOW -----------
% --  ---------------- -----------

%-- init

%easyupMVPA_init('nmbCores', cpu_cores, 'quietMode', false);

%-- get the time (needed later for the analysis subfolder)
timeVector = clock;

%-- create configParameters STRUCT (collect all specified variables and
%put them in a structure)
configSpecifiedVariables=who;
for i=1:numel(configSpecifiedVariables)
    if strcmp(configSpecifiedVariables{i},'configParameters') && strcmp(configSpecifiedVariables{i},'configSpecifiedVariables')
    else
        if numel(eval(configSpecifiedVariables{i})) > 1000
        else
            eval(['configParameters.',configSpecifiedVariables{i},'=',configSpecifiedVariables{i},';']);
        end
    end
end

%-- put the configParameters structure into the parser
myDataset = parseconfigParameters(configParameters);


%PPPPPPPAAAAAAAAARRRRRRRRRSSSSSSSSSSSSSEEEEEEEEEERRRRRRRRRRRRRRRRRRRr
