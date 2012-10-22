% dddddddddddddddddddd
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [myDataset] = setBrainMaskbyconfigParameters(myDataset)




function [myDataset] = setBrainMaskbyconfigParameters(myDataset)



currentAnalysisDirectory = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.dir_analysis_mvpa);





if myDataset.configParameters.obtainmask %automated mask generation? (mask will be loadad if generated earlier)
    if(~easyupMVPA_getGlobals('quietMode'))
        disp('INFO: Brain mask: Using an automatically generated mask...');
    end
    if exist(fullfile(currentAnalysisDirectory,'wholebrainmask.nii')) == 0
        if(~easyupMVPA_getGlobals('quietMode'))
            disp('INFO: Brain mask: Obtaining the mask ...');
        end
        getWholeBrainMask(myDataset.configParameters);
        myDataset = setDataset_mask_ByImageFile(myDataset, fullfile(currentAnalysisDirectory,'wholebrainmask.nii'));
        
    else
        myDataset = setDataset_mask_ByImageFile(myDataset, fullfile(currentAnalysisDirectory,'wholebrainmask.nii'));

    end
else
    
    
    

    %use a global mask for ALL participants or a personal one for each
    %subject?
    if numel(strfind(myDataset.configParameters.mask,'scr/')) || numel(strfind(myDataset.configParameters.mask,'SCR/')) || numel(strfind(myDataset.configParameters.mask,':\'))
        maskfilename = myDataset.configParameters.mask;
        if(~easyupMVPA_getGlobals('quietMode'))
            disp(['INFO: Brain mask: Using *the same* mask for *all* subjects: '])
        end
        
    else
        maskfilename = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.mask);
        if(~easyupMVPA_getGlobals('quietMode'))
            disp(['INFO: Brain mask: Using mask: ']);
        end
    end
    myDataset = setDataset_mask_ByImageFile(myDataset,maskfilename);
end
    
