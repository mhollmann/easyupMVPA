% Prints the content of the result struct (result of prediction, LOOCV, RFE) on the screen.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%
%   printResultStruct(resultStruct)
%
%   This methods prints the whole info about dataset at the screen.
%
% Parameters:
%   dataset  - the dataset to print info 
%
% Returns:
%
%
% Comments:
%
function printResultStruct(resultStruct)

     disp('*** RESULT STRUCT:  ***');
     disp(['Number of Tests: ', num2str(resultStruct.nmbTests)])
     disp(['Accuracy:        ', num2str(resultStruct.accuracy), ' %'])
     disp(['True Positives:  ', num2str(resultStruct.TP)])
     disp(['True Negatives:  ', num2str(resultStruct.TN)])
     disp(['False Positives: ', num2str(resultStruct.FP)])
     disp(['False Negatives: ', num2str(resultStruct.FN)])
     disp(['Sensitivity:     ', num2str(resultStruct.sensitivity),'  (TP/TP+FN = Proportion of true positives to all positives).'])
     disp(['Specificity:     ', num2str(resultStruct.specificity),'  (TN/TN+FP = Proportion of true negatives to all negatives).'])
     if(isfield(resultStruct, 'infoString'))
       disp(['Additional Info: ', resultStruct.infoString])
     end
     disp('************************');
     
     if(isfield(resultStruct, 'innerResultStruct') && ~isempty(resultStruct.innerResultStruct))
       disp(' ');
       disp('*** INNER RESULT STRUCT HOLDING ADDITIONAL INFORMATION ***');
       printResultStruct(resultStruct.innerResultStruct);
     end
     
end