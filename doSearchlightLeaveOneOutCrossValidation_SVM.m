% Does a searchlight classification using a SVM and Leave One Out Cross Validation.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%  This high-level function implements a simple prediction using a Support Vector Machine (SVM).
%
%
% Parameters:
%   dataset     - the datset holding the test set to do prediction for
%   svmModel    - the model learned via train_SVM
%
% Returns:
%   resultStruct   - The struct holding the classification results:
%                    resultStruct.nmbTests     (the number of samples tested for this result)
%                    resultStruct.accuracy     (percentual value of correct predictions (correct * 100 / nmbSamples))
%                    resultStruct.sensitivity  (TP/TP+FN = Proportion of true positives to all positives)
%                    resultStruct.specificity  (TN/TN+FP = Proportion of true negatives to all negatives)
%                    resultStruct.TP           (True positives = all correct predicted in class 1)
%                    resultStruct.TN           (True negatives = all correct predicted in class 2)
%                    resultStruct.FP           (False positives = all incorrect predicted in class 1)
%                    resultStruct.FN           (False negatives = all incorrect predicted in class 2)
%
%   probEstimates  - Probability estimates if model was learned with Option 'probEstimates' = 1, otherwise empty matrix
%
% Comments:
%
function [dataset, resultStruct, avgWeights] = doSearchlightLeaveOneOutCrossValidation_SVM(dataset, dataSplitter, kernelMode, costParam, paramStruct)

end