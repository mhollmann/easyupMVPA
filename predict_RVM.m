% (UNDER CONSTRUCTION) Predicts classes of the samples in dataset using the given SVM-model.
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
% Comments:
%
function [resultStruct, probEstimates] = predict_RVM(dataset, rvmModel)


  resultStruct = {};  

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
     dataset = setDataset_featureSelectionMap_ByMatrix(dataset, ones(sizeData(1:end-1))); 
    end
    [dataset, testData2D] = selectFeaturesBySelectionMap(dataset);
  end

  
  %size testData :  nmbSamples, nmbFeatures
  %rvmModel.X(rvmModel.relevantIndices,:) results in [nmbUsed, nmbFeatures]

  
  nmbSamples = length(dataset.classIDs);
    
  %compute the basis set
  basisSet	= sbl_kernelFunction(testData2D,rvmModel.sparseTrainData,rvmModel.kernelType,rvmModel.kernelWidth);

  %Compute the inferred prediction function
  %y		= trainData2D*weightVals;
  %y_l	= double(SB2_Sigmoid(y)>0.5);

   %y_rvm	  = basisSet*rvmModel.weights;
   %y_class	= double(SB2_Sigmoid(y_rvm)>0.5);
   

   
   probEstimates = 1./(1+exp(-basisSet*rvmModel.weights));
   y_class	     = double(probEstimates>0.5);

   %disp(y_class');
   %disp(dataset.classIDs);
   
   cVec = y_class' == dataset.classIDs;
%    disp(['correct: ', num2str(sum(cVec)), ' of ', num2str()]);
   
 
   nmbTruePos  = sum(cVec(dataset.classIDs==1));
   nmbTrueNeg  = sum(cVec(dataset.classIDs==0));
   nmbFalsePos = sum(dataset.classIDs) - sum(y_class);
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
   
   
  % keyboard;
   
   %y_prob		= 1./(1+exp(-basisSet*rvmModel.weights));
   
   % apply sigmoid for probabilities
   %y_prob		= 1./(1+exp(-PHI*RVM.weights)); 
  
  


end