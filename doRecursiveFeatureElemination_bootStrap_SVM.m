% [DO NOT USE: UNDER DEVELOPEMENT!] Recursively removes features according to SVM classification weights obtained by bootstrapping.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%   
%   [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut, kernelMode, costParam, [paramStruct])
%
%   This high-level function implements a recursive feature elemination (RFE) using a Support Vector Machine (SVM).
%   This means that voxels that hold no or just little information for classification are removed to recursively improve
%   classification performance. 
%   The algorithm is basically a LeaveOneOutCrossValidation and the given dataSplitter defines the test/training 
%   sets in every single run. Which  means that the  test  results  returned  with "resultStruct" can  be  used to  estimate 
%   the overall classification performance of the given SVM-classifier.
%
%   For every iteration the following is done:
%     For all splits(nmbSplits) in "dataSplitter" the data is splitted into training and test sets. In every training set the 
%     features are selected SEPERATELY according to the input parameter "thresholdPercentOfFeaturesOut". And the test set is 
%     classified with the trained model.
%
%   That means there are nSplits models trained and every model has its own INDEPENDENT feature selection. For all this single 
%   selections a selectionMap is the result and these are returnd by the parameter "avg_rfe_featureSelectionMap". The same holds
%   for the weights trained in every single model. This ensures that no information between test and training data is transferred.
%
%   The average maps for weights and featureSelectionMaps are returned also:
%   BUT BE AWARE: These may not be used as feature selection maps for further classification on the same dataset, because of 
%   information transfer!
%
%   To see the classification-performance during the RFE see the output on command line (If "quietMode" is switched off).
%   
%   The returned classification results may not be the best because always the LAST iteration sets the returned resultStruct. If
%   an intermediate iteration is best run this function again with this iteration-nmb-1 as input and then use the results.
%
%   The best performing iteration having the lowest number of features and all in between will be stored in additional info string
%   in resultStruct (May be displayed using "printResultStruct(...)") and as capsuled resultStruct in field resultStruct.innerResultStruct.
%
% Parameters:
%   dataset                       - the datset to set the classIDs for
%   nmbIterations                 - iterations (how often is basic feature set reduced by "thresholdPercentOfFeaturesOut")
%   thresholdPercentOfFeaturesOut - percentual value of non-zero elements in weight map that should be cut (suggestion btw. 10 and 50)
%   dataSplitter                  - describes the splitting of the data in the background LOOCV
%   kernelMode                    - Kernels: ['linear', 'polynomial', 'radial', 'sigmoid']
%   costParam                     - The slack variable C in SVM (range 0 to 1  0 = low cost, 1 = highest costs). 
%                                   It defines the costs for misclassification (How strongly are outliers punished?).
%   paramStruct                   - [optional] i.e. {"degree", 3}
%
% Returns:
%   dataset                     - In this dataset the field "featureSelectionMap" is updated to the result of this RFE
%   resultStruct                - The struct holding the classification results: 
%                                 resultStruct.nmbTests     (the number of samples tested for this result)
%                                 resultStruct.accuracy     (percentual value of correct predictions (correct * 100 / nmbSamples))
%                                 resultStruct.sensitivity  (TP/TP+FN = Proportion of true positives to all positives)
%                                 resultStruct.specificity  (TN/TN+FP = Proportion of true negatives to all negatives)
%                                 resultStruct.TP           (True positives = all correct predicted in class 1)
%                                 resultStruct.TN           (True negatives = all correct predicted in class 2)
%                                 resultStruct.FP           (False positives = all incorrect predicted in class 1)
%                                 resultStruct.FN           (False negatives = all incorrect predicted in class 2)
%   avg_rfe_weightMap           - Average map of the weights of the models trained (always determined in the last iteration - Dimension: Featurespace)
%   avg_rfe_featureSelectionMap - Average map of the selected features of the models trained (always determined in the last iteration - Dimension: Featurespace)
%   rfe_weightMaps              - All  all weight maps that are used in the last iteration (Dimension: Featurespace x nmbSplits)
%   rfe_featureSelectionMaps    - All feature-SelectionMaps (elements just 0 or 1) that are used in the last iteration (Dimension: Featurespace x nmbSplits)
%
% Comments:
%
function [dataset, resultStruct] = doRecursiveFeatureElemination_bootStrap_SVM(dataset, nmbBootSteps)

  if( ~exist('dataset','var') )%|| ~exist('nmbIterations','var') || ~exist('thresholdPercentOfFeaturesOut','var') || ~exist('dataSplitter','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    error('Usage of doRecursiveFeatureElemination_bootStrap_SVM: [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  sizeData   = size(dataset.data);
  if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
    [dataset, trainData2D] = selectFeaturesBySelectionMap(dataset);
  elseif(isfield(dataset,'mask3D') && ~isempty(dataset.mask))
    dataset = setDataset_featureSelectionMap_ByMatrix(dataset, dataset.mask);
    [dataset, trainData2D] = selectFeaturesBySelectionMap(dataset);
  else
    %all features are used
    if(dataset.is2D)
     dataset = setDataset_featureSelectionMap_ByMatrix(dataset, ones(sizeData(1:end-1),1));
    else
     dataset = setDataset_featureSelectionMap_ByMatrix(dataset, ones(sizeData(1:end-1))); 
    end
    [dataset, trainData2D] = selectFeaturesBySelectionMap(dataset);
  end
  
  resultStruct = {};
  
  %keyboard;
  
 [descrim_vec, descrimsMean, boot_se, boot_CI, boot_ipred, boot_gpred] = ... 
    svm_632boot_unique( trainData2D, dataset.classIDs, nmbBootSteps, 1, true );
  
  featureSet = prod(boot_CI,2)>0;
  disp('Number of selected features boot_CI:');
  disp(sum(featureSet));
  
%   for i=1:numel(descrimsMean)
%       disp(['Mean: ', num2str(descrimsMean(i)), ' SE: ', num2str(boot_se(i)), ' CIL: ',num2str(boot_CI(i,1)),' CIH: ',num2str(boot_CI(i,2))]);
%   end
   

  thresholdPercentOfFeaturesOut = 20;

  %remove the % features with lowest mean discriminability 
  [mapOut, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(descrimsMean, thresholdPercentOfFeaturesOut);
  
  featureSet = mapOut ~= 0;
  disp('Number of selected features:');
  disp(totalNmbOut);
  
  if(sum(featureSet) == 0)
    warning('No significant features selected!');
  end
  
  
  ds2D = getEmpty2DDataset();
  ds2D = setDataset_data_ByMatrix(ds2D, trainData2D');
  ds2D = setDataset_chunks_ByVector(ds2D, ones(size(trainData2D,1),1));
  ds2D = setDataset_classIDs_ByVector(ds2D, dataset.classIDs);
  ds2D = setDataset_featureSelectionMap_ByMatrix(ds2D, featureSet);
  
  printDatasetInfo(ds2D);
  
  splitterOSO = getDataSplitter(ds2D, 'oneSampleOut');
  [datasetTest, resultStruct, avgWeights] = doLeaveOneOutCrossValidation_SVM(ds2D, splitterOSO, 'linear', 0.5);
  
  printResultStruct(resultStruct);
  %disp(feature)
  %keyboard;
  %[descrim, SV, misclass]=svm_kfcv(featurespace(:,featureset==1), dataset.classIDs, K, C, autoscale);
  
 %keyboard;
  
end


%member function to select a certain amount of elements from a 3D map
function [mapOut, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(mapIn, thresholdPercentOfElementsOut)

  mapOut = mapIn;
  
  %overall nmb of non-zero-elements
  totalNmbIn = sum(mapIn(:)~=0);
  
  %avgWeights3D = avgWeights3D + tmpWeights3D*(1/nmbSamples);
  %now set the rfe selection for this split to the according features
  %nmb of voxels to throw out
  nmbOut = floor(totalNmbIn*thresholdPercentOfElementsOut/100);
  
  if(nmbOut == 0)
    if(totalNmbIn<2)
      warning('Could not select features by given threshold. Not enough input features.');
      return;
    else
      nmbOut = 2;
    end  
  end
  
  %a sorted array to get the value for selection
  tmp = sort(abs(mapIn(mapIn(:)~=0)));
    
  selThresh = tmp(nmbOut);
  
  tmpMap = abs(mapIn);
  tmpMap(tmpMap<selThresh) = 0;
  mapOut(tmpMap==0) = 0;
  
  totalNmbOut = sum(mapOut(:)~=0);
  
end %end member_getSelectedMapByThresholdPercentOut