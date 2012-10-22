% Predicts classes of the samples in dataset using the given SVM-model.
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
function [resultStruct, probEstimates] = predict_SVM(dataset, svmModel)

  nmbSamples    = length(dataset.chunks);
  probEstimates = []; 
  
  
  %get the data as 2D array by using the selection-map or the mask 
  sizeData   = size(dataset.data);
  if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
    [dataset, testData2D] = selectFeaturesBySelectionMap(dataset);
  elseif(isfield(dataset,'mask3D') && ~isempty(dataset.mask))
    dataset = setDataset_featureSelectionMap_ByMatrix(dataset, dataset.mask);
    [dataset, testData2D] = selectFeaturesBySelectionMap(dataset);
  else
    %all features are used
    if(dataset.is2D)
     dataset = setDataset_featureSelectionMap_ByMatrix(dataset, ones(sizeData(1:end-1),1));
    else
     dataset = setDataset_featureSelectionMap_ByMatrix(dataset, ones(sizeData(1), sizeData(2), sizeData(3))); 
    end
    [dataset, testData2D] = selectFeaturesBySelectionMap(dataset);
  end

   %predict the class ID of the test data
   %check if the model was trained with probability estimates
   if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
     [predicted_labels, accuracy, probEstimates] = svmpredict(double(dataset.classIDs)', testData2D, svmModel, '-b 1');
   else
     [predicted_labels, accuracy, decisionValues] = svmpredict(double(dataset.classIDs)', testData2D, svmModel); 
   end

   cVec = predicted_labels' == dataset.classIDs;
   
 
   nmbTruePos  = sum(cVec(dataset.classIDs==1));
   nmbTrueNeg  = sum(cVec(dataset.classIDs==0));
   nmbFalsePos = sum(dataset.classIDs) - sum(predicted_labels);
   if(nmbFalsePos < 0)
     nmbFalsePos = abs(nmbFalsePos);
   else
     nmbFalsePos = 0;
   end
   nmbFalseNeg = nmbSamples - nmbTrueNeg - nmbFalsePos - nmbTruePos;
     
   accuracy    = sum(cVec)/nmbSamples*100;
   sensitivity = nmbTruePos/(nmbTruePos+nmbFalseNeg); 
   specificity = nmbTrueNeg/(nmbTrueNeg+nmbFalsePos);

   resultStruct             = {};
   resultStruct.nmbTests    = nmbSamples;
   resultStruct.accuracy    = accuracy;
   resultStruct.sensitivity = sensitivity;
   resultStruct.specificity = specificity;
   resultStruct.TP          = nmbTruePos;
   resultStruct.TN          = nmbTrueNeg;
   resultStruct.FP          = nmbFalsePos;
   resultStruct.FN          = nmbFalseNeg;
   
end