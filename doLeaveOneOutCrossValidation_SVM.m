% Implements a Leave One Out Cross Validation (LOOCV) using a Support Vector Machine Classifier (parallel execution).
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   [dataset, resultStruct, avgWeights] = doLeaveOneOutCrossValidation_SVM(dataset, dataSplitter, svmType, kernelMode, costParam, [paramStruct])
%
%   This highlevel function does the Leave One Out Cross Validation (LOOCV) of a complete dataset(DS) using a Support Vector Machine.
%   LOOCV means a repeated classification in the dataset to get a general performance measure.
%   Which samples are used in every repetition as training and testing set is described in the given dataSplitter-struct.
%
%   For a "oneSampleOut"-splitter LOOCV means that for every single sample in DS the DS without this sample is trained and the test is 
%   done on the excluded sample. This means the classifier is trained n times for n samples, which may be very time consuming.
%
%   In every new training the weights used in the model are saved and the function returns the average weights of all single models trained.
%
%   Explanation of resultStruct fields:
%   Accuracy:    Percentual value of correct predictions (correct * 100 / nmbSamples)
%   Sensitivity: TP/TP+FN (Proportion of true positives to all positives)
%   Specificity: TN/TN+FP (Proportion of true negatives to all negatives)
%
%   If possible, the function executes parallelized depending on settings defined via easyupMVPA_init().
%
%
% Parameters:
%   dataset        - The dataset to work on  (all samples are included in LOOCV)
%   dataSplitter   - describes the splitting of the data in LOOCV
%   svmType        - Types:
%                     ['classification', 'regression_epsilon', 'regression_nu']
%   kernelMode     - Kernels: ['linear', 'polynomial', 'radial', 'sigmoid']
%   costParam      - The slack variable C in SVM (range 0 to 1  0 = low cost, 1 = highest costs). 
%                    It defines the costs for misclassification (How strongly are outliers punished?).
%   paramStruct    - [optional] - i.e. {"degree", 3, "gamma", 0.5, "coeff0", 0.5, "nu", 0.3, "epsilon", 0.01, "probEstimates", 1}
%                     possible fields:
%                    'degree'        : default=3     Describes the exponent in polynomial kernel function
%                    'gamma'         : default=1/k   The gamma-factor in kernel function
%                    'coef0'         : default=0     The coefficient summand in kernel function
%                    'nu'            : default=0.5   The nu parameter of regression_nu
%                    'epsilon'       : default=0.001 The epsilon parameter of regression_epsilon
%                    'probEstimates' : default=0     1 if probabilistic estimates should be computed, 0 if not

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
%   avgWeights     - a 3D or 1D map containing the average weights used in the LOOCV (See Description of this function)
%
%
% Comments:
%
function [dataset, resultStruct, avgWeights] = doLeaveOneOutCrossValidation_SVM(dataset, dataSplitter, svmType, kernelMode, costParam, paramStruct)
 

  if( ~exist('dataset','var') || ~exist('dataSplitter','var') || ~exist('svmType','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    keyboard;
    error('Usage of doLeaveOneOutCrossValidation_SVM: [dataset, classificationResult, avgWeights] = doLeaveOneOutCrossValidation_SVM(dataset, dataSplitter, svmType - [classification, regression_epsilon, regression_nu], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %extract the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, paramStruct);
  end
  if( ~paramStructIsValid)
    error('Usage of doLeaveOneOutCrossValidation_SVM: [dataset, classificationResult, avgWeights] = doLeaveOneOutCrossValidation_SVM(dataset, dataSplitter, svmType - [classification, regression_epsilon, regression_nu], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
 
  
  %use quiet mode (no outputs)
  cmdString = [cmdString, ' -q '];
  
   if(dataset.is4D)
     sizeData   = size(dataset.data);
     avgWeights = zeros(sizeData(1),sizeData(2),sizeData(3),1);
     tmpWeightsSize = [sizeData(1),sizeData(2),sizeData(3),1];
   elseif(dataset.is2D)
     sizeData   = size(dataset.data);
     avgWeights = zeros(sizeData(1),1);
     tmpWeightsSize = [sizeData(1),1];
   else
    error('LOOCV: Please check the dataset: field "type" is not defined!');
   end
   
   nmbSamples = length(dataset.chunks);
   
   nmbCorrect    = 0;
   
   nmbTruePosAll    = 0;
   nmbTrueNegAll    = 0;
   nmbFalsePosAll   = 0;
   nmbFalseNegAll   = 0; 
   
   %extract the number of splits that are used
   nmbSplits = size(dataSplitter.splitMatrix,1);
   nmbTests  = numel(find(dataSplitter.splitMatrix(:)==1));
   
   %variable to store predicted labels
   predictedLabels = {};
   
   localQuietMode = easyupMVPA_getGlobals('quietMode');
   
   progressIndices = [1 1:nmbSplits];
   
   if(~localQuietMode)
     disp(['Running SVM Leave One Out Cross Validation with command string: ',cmdString,' ...']);
     % create a progress display that works also for parallel loops 
     if(nmbSplits <50)
       disp(['0%', num2str(repmat(' ',1,nmbSplits-1)),'100%']);
       progressIndices = [1 1:nmbSplits];
     else
       disp(['0%', num2str(repmat(' ',1,50)),'100%']);
       %create a vector with floored indicees for repetitions
       progressIndices = [1 1:nmbSplits];
       progressIndices = floor(progressIndices*(50/nmbSplits));
     end
     fprintf('    ');
   end

   nmbTests = 0;
   
   splitMatrix = dataSplitter.splitMatrix;
   
   
   parfor i=1:nmbSplits
     
     trainData2D = [];
     testData2D  = [];  
     
     %at first split the dataset according to given splitting
     ds1Indices = splitMatrix(i,:) == 1; %test data
     ds2Indices = splitMatrix(i,:) == 2; %train data


     [ds1, ds2] = splitDataset(dataset, ds1Indices, ds2Indices );     
     
     if(dataset.is4D)
       [ds2, trainData2D] = member_getMasked2DData_from4D(ds2);
     elseif(dataset.is2D)
       [ds2, trainData2D] = member_getMasked2DData_from2D(ds2);
     end
     
     %train on all samples with one left out
     svmModel  = svmtrain(double(ds2.classIDs)', trainData2D, cmdString);
     
     %extract the weights
     weights = svmModel.SVs' * svmModel.sv_coef;
     
     tmpWeights = zeros(tmpWeightsSize);
     
     if(dataset.is4D)
       
       tmpWeights(ds2.featureSelectionMap >0) = weights;
       avgWeights = avgWeights + tmpWeights*(1/nmbSamples);
       
       %select the data for the testset
       [ds1, testData2D] = member_getMasked2DData_from4D(ds1);
       
     elseif(dataset.is2D)
       
       tmpWeights(ds2.featureSelectionMap >0) = weights;
       avgWeights = avgWeights + tmpWeights*(1/nmbSamples);
       
        %select the data for the testset
       [ds1, testData2D] = member_getMasked2DData_from2D(ds1);
       
     end
     
     %predict the class ID of the test data
     %check if the model was trained with probability estimates
     if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
       [predicted_labels, accuracy, probEst] = svmpredict_modifiedMH(double(ds1.classIDs)', testData2D, svmModel, '-b 1');
       %disp(probEstimatesVector);
     else
       [predicted_labels, accuracy, probEst] = svmpredict_modifiedMH(double(ds1.classIDs)', testData2D, svmModel); 
     end
     
     %remember predicted labels
     predictedLabels(i).predLabels = predicted_labels;
     
     if(strcmp(svmType, 'classification'))

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
       
     end
     
     
     nmbTests         = nmbTests+length(ds1.classIDs);

     if(~localQuietMode)
       if(progressIndices(i)<progressIndices(i+1))
         fprintf('\b\b*');
         disp(['' 0]);
       end
     end

   end % end parfor
   
   if(~localQuietMode)
     fprintf('\n');
   end
   
   predLabelsVec = [];
   for i=1:length(predictedLabels)
     %sorry it is growing but hard to avoid 
     predLabelsVec = horzcat(predLabelsVec, predictedLabels(i).predLabels);
   end
   
   if(strcmp(svmType, 'classification'))

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
     resultStruct.predictedClassIDs = predLabelsVec;
     
   else %must be regression
     
     resultStruct             = {};
     resultStruct.nmbTests    = nmbTests;
     resultStruct.accuracy    = NaN;
     resultStruct.sensitivity = NaN;
     resultStruct.specificity = NaN;
     resultStruct.TP          = NaN;
     resultStruct.TN          = NaN;
     resultStruct.FP          = NaN;
     resultStruct.FN          = NaN;
     resultStruct.predictedClassIDs = predLabelsVec;
     
   end
   
end



function [ds, maskedData2D] = member_getMasked2DData_from4D(ds)
  sizeMaskedData2D = size(ds.data);
     
  %get the data as 2D array by using the 3D mask or selection map
  %the function selectFeaturesBySelectionMap sets the field dataset.featureSelectionMap3D
  if(isfield(ds,'featureSelectionMap') && ~isempty(ds.featureSelectionMap))
    [ds, maskedData2D] = selectFeaturesBySelectionMap(ds);
  elseif(isfield(ds,'mask') && ~isempty(ds.mask))      
    ds = setDataset_featureSelectionMap_ByMatrix(ds, ds.mask);
    [ds, maskedData2D] = selectFeaturesBySelectionMap(ds);
  else
    ds = setDataset_featureSelectionMap_ByMatrix(ds, ones(sizeMaskedData2D(1)*sizeMaskedData2D(2)*sizeMaskedData2D(3), 'uint8'));
    sizeData = size(ds.data);
    maskedData2D = reshape(ds.data, sizeData(1)*sizeData(2)*sizeData(3), sizeData(4))';    
  end
 
end % end function  member_getMasked2DData_from4D


function [ds, maskedData2D] = member_getMasked2DData_from2D(ds)
  sizeMaskedData2D = size(ds.data);
     
  %numel(find(dataset.featureSelectionMap1D > 0)
  %get the data as 2D array by using the 1D mask1D
  %the function selectFeaturesBySelectionMap sets the field dataset.featureSelectionMap3D
  if(isfield(ds,'featureSelectionMap') && ~isempty(ds.featureSelectionMap))
    [ds, maskedData2D] = selectFeaturesBySelectionMap(ds);
  elseif(isfield(ds,'mask') && ~isempty(ds.mask))      
    ds = setDataset_featureSelectionMap_ByMatrix(ds, ds.mask);
    [ds, maskedData2D] = selectFeaturesBySelectionMap(ds);
  else
    ds = setDataset_featureSelectionMap_ByMatrix(ds, ones(sizeMaskedData2D(1),1, 'uint8'));
    maskedData2D = ds.data';
  end

end % end function  member_getMasked2DData_from2D
