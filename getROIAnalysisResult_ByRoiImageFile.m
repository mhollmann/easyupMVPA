% DEPRECATED!
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%
% Parameters:
%
% Returns:
%
% Comments:
%
%
function [roiAnalysisResult] = getROIAnalysisResult_ByRoiImageFile(dataset, hdrFilename)

 roiAnalysisResult = struct('meanValuesClasses', []);
 
 nmbSamples = length(dataset.classIDs);
 
 %Read the data (nii or hdr)
 %the nii loading function is slightly changed for better memory performance
 dataNii = load_untouch_nii(hdrFilename);
 roiMask = dataNii.img;
 
 showDataAsImage(roiMask, 'Mask');
 
 nmbVoxelsInROI = numel(find(roiMask>0));
 
 %ROI infos
 disp(['Roi mask with ',num2str(nmbVoxelsInROI),' elements...']);
 
 uniqueClasses        = unique(dataset.classIDs);
 nmbOfClasses         = length(uniqueClasses);
 meanValuesClassesTmp = zeros(nmbOfClasses, nmbSamples);
 
 meanValuesClasses = {};
 
 for j=1:nmbOfClasses
   
   ds = selectSamples(dataset, ['classIDs==', num2str(uniqueClasses(j))]);
   nmbSamples = length(ds.classIDs);
      
   for i=1:nmbSamples
     masked_funcData3D_ROI = ds.data4D(:,:,:,i);
     masked_funcData3D_ROI(roiMask == 0) = 0;
     meanROI     = sum(masked_funcData3D_ROI(:)) / nmbVoxelsInROI;
     meanValuesClassesTmp(j,i) = meanROI;
     disp(meanROI);
   end
   
   meanValuesClasses{j} = meanValuesClassesTmp(:);
 end
 
 roiAnalysisResult.meanValuesClasses = meanValuesClasses;
 
 keyboard;
end