% !JUST FOR INTERNAL USE! Implements a very fast Leave One Out Cross Validation (LOOCV) just for 2D case using a Support Vector Machine Classifier.
%
% Author: Maurice Hollmann
% Date  : 05/11
%
% Description:
%
%   [dataset, resultStruct] = doLeaveOneOutCrossValidation_SVM_2DforceQuiet(dataset, dataSplitter, svmCommandString)
%
%   Low level - just for internal use!  
%   Just 2D input!
%   No weights are computed!
%   Forced to be quiet!
%   No parallel execution, because it is expected that this one is called from inside 
%   a parfor loop and a nested parfor would slow down the execution.
%
% Parameters:
%   dataset          - The dataset to work on  (all samples are included in LOOCV)
%   dataSplitter     - describes the splitting of the data in LOOCV
%   svmCommandString - the svm command string for exmaple '-t 0 -c 0.5' 
%
% Returns:
%   dataset        - the datset that has been the input 
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
function [dataset, resultStruct] = doLeaveOneOutCrossValidation_SVM_2DforceQuiet(dataset, dataSplitter, svmCommandString)
   
   nmbCorrect       = 0;
   nmbTruePosAll    = 0;
   nmbTrueNegAll    = 0;
   nmbFalsePosAll   = 0;
   nmbFalseNegAll   = 0; 
   
   %extract the number of splits that are used
   nmbSplits = size(dataSplitter.splitMatrix,1);
   nmbTests = 0;
   splitMatrix = dataSplitter.splitMatrix;
   
   %parfor i=1:nmbSplits
   for i=1:nmbSplits
          
     %at first split the dataset according to given splitting
     ds1Indices = splitMatrix(i,:) == 1; %test data
     ds2Indices = splitMatrix(i,:) == 2; %train data

     [ds1, ds2] = splitDataset(dataset, ds1Indices, ds2Indices );
           
     %train on all samples with one left out
     svmModel  = svmtrain(double(ds2.classIDs)', double(ds2.data)', svmCommandString);
     
      %predict the class ID of the test data
     %check if the model was trained with probability estimates
     if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
       [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', double(ds1.data)', svmModel, '-b 1');
       %disp(probEstimatesVector);
     else
       [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', double(ds1.data)', svmModel); 
     end
     
     cVec = predicted_labels' == ds1.classIDs;
     
     nmbTruePos  = sum(cVec(ds1.classIDs==1));
     nmbTrueNeg  = sum(cVec(ds1.classIDs==0));
     nmbFalsePos = sum(ds1.classIDs) - sum(predicted_labels);
     if(nmbFalsePos < 0)
       nmbFalsePos = abs(nmbFalsePos);
     else
       nmbFalsePos = 0;
     end
     nmbFalseNeg = length(ds1.classIDs) - nmbTrueNeg - nmbFalsePos - nmbTruePos;

     nmbCorrect = nmbCorrect+nmbTruePos+nmbTrueNeg;

     nmbTruePosAll    = nmbTruePosAll+nmbTruePos;
     nmbTrueNegAll    = nmbTrueNegAll+nmbTrueNeg;
     nmbFalsePosAll   = nmbFalsePosAll+nmbFalsePos;
     nmbFalseNegAll   = nmbFalseNegAll+nmbFalseNeg; 
     nmbTests         = nmbTests+length(ds1.classIDs);

   end % end parfor
   
   accuracy    = nmbCorrect/nmbTests*100;
   sensitivity = nmbTruePosAll/(nmbTruePosAll+nmbFalseNegAll); 
   specificity = nmbTrueNegAll/(nmbTrueNegAll+nmbFalsePosAll);
   
   resultStruct             = {};
   resultStruct.nmbTests    = nmbTests;
   resultStruct.accuracy    = accuracy;
   resultStruct.sensitivity = sensitivity;
   resultStruct.specificity = specificity;
   resultStruct.TP          = nmbTruePosAll;
   resultStruct.TN          = nmbTrueNegAll;
   resultStruct.FP          = nmbFalsePosAll;
   resultStruct.FN          = nmbFalseNegAll;
   
end
