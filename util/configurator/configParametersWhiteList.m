function [in_whitelist] = configParametersWhiteList(variablename);

allowed_vars = {'Comments'
    'baseDirectory'
    'subjectlist'
    'conditions'
    'comparisons'
    'cpu_cores'
    'dataFileFormat'
    'TR'
    'mask'
    'detrending'
    'highpassFilter'
    'z_scoring'
    'savePreProcessing'
    'useSpmBetaFiles'
    'sameTrialDurations'
    'removeTrailingTransitions'
    'removeTailTransitions'
    'unitofdesign'
    'spm_onsets_specification'
    'averaging'
    'averagingSplitMethod'
    'timeVector'
    'runs'
    
    'SearchLightDiameter'
    'BaseMethod'
    'loocvSplitMethod'
    'secondLevelSmoothingFWHM'
    'accuracymapsavename'
    'customAnalysisSubfolder'
    
    'Searchlightcachesize'
    
    
    
    
    
    };
% 
% for i=1:100
%     allowed_vars{numel(allowed_vars)+i} = ['runs(',num2str(i),').directory'];
% end

in_whitelist= sum(strcmp(allowed_vars,variablename));




