% ! DEPRECATED ! Recursively removes features according to SVM classification weights from LOOCV.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%   
%   [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut, kernelMode, costParam, [paramStruct])
%
%Algorithm:
%
% For splitterSplits 1:N do
%   
%    trainDataset  = select just trainData
%    loocvSplitter = one to loocvSplitNmb
%
%    avgLoocvWeights =  doRecursiveFeatureElemination(trainDataset, loocvSplitter ...)
%    
%    
%
% end
%
%
% Parameters:
%   dataset                       - the datset to set the classIDs for
%   nmbIterations                 - iterations (how often is basic feature set reduced by "thresholdPercentOfFeaturesOut")
%   thresholdPercentOfFeaturesOut - percentual value of non-zero elements in weight map that should be cut (suggestion btw. 10 and 50)
%   dataSplitter                  - describes the splitting of the data in the background LOOCV
%   kernelMode                    - Kernels: ['linear', 'polynomial', 'radial', 'sigmoid']
%   costParam                     - The slack variable C in SVM (range 0 to 1  0 = low cost, 1 = highest costs). 
%                                   It defines the costs for misclassification (How strongly are outliers punished?).
%   paramStruct                   - (optional)s i.e. {"degree", 3}
%
% Returns:
%   dataset                     - In this dataset the field "featureSelectionMap3D" is updated to the result of this RFE
%   resultStruct                - The struct holding the classification results: 
%                                 resultStruct.nmbTests     (the number of samples tested for this result)
%                                 resultStruct.accuracy     (percentual value of correct predictions (correct * 100 / nmbSamples))
%                                 resultStruct.sensitivity  (TP/TP+FN = Proportion of true positives to all positives)
%                                 resultStruct.specificity  (TN/TN+FP = Proportion of true negatives to all negatives)
%                                 resultStruct.TP           (True positives = all correct predicted in class 1)
%                                 resultStruct.TN           (True negatives = all correct predicted in class 2)
%                                 resultStruct.FP           (False positives = all incorrect predicted in class 1)
%                                 resultStruct.FN           (False negatives = all incorrect predicted in class 2)
%   avg_rfe_weightMap           - Average map of the weights of the models trained (always the last selection)
%   avg_rfe_featureSelectionMap - Average map of the selected features of the models trained (always the last selection)
%
%
%
% Comments:
%
function [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM_nestedLOOCV(dataset, nmbIterations, thresholdPercentOfFeaturesOut, nmbSplitsLOOCV, dataSplitter, kernelMode, costParam, paramStruct)
  
  if( ~exist('dataset','var') || ~exist('nmbIterations','var') || ~exist('thresholdPercentOfFeaturesOut','var') || ~exist('dataSplitter','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  
  loocvNmbSplits = nmbSplitsLOOCV;
  
  
  
  %extractt the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, paramStruct);
  end
  if( ~paramStructIsValid)
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %use quiet mode (no outputs)
  cmdString = [cmdString, ' -q '];
  
  if(thresholdPercentOfFeaturesOut < 0 || thresholdPercentOfFeaturesOut > 100)
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
      
  if(~easyupMVPA_getGlobals('quietMode'))
    disp('Running Recursive Feature Elemination (This may take a while!) ...');
  end
  
  %extract the number of splits that are used
  nmbSplits = size(dataSplitter.splitMatrix,1);
  
  %check if there is a global feature selection map is given, if yes all 
  %the following steps are computed solely on these selected features
  if(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset4D'))
    
    %4D case
    if(isfield(dataset,'featureSelectionMap3D') && ~isempty(dataset.featureSelectionMap3D))
      globalSelectionMap = dataset.featureSelectionMap3D;
    elseif(isfield(dataset,'mask3D') && ~isempty(dataset.mask3D))
      globalSelectionMap = dataset.mask3D;
    else
      globalSelectionMap = ones(sizeTrainData4D(1)*sizeTrainData4D(2)*sizeTrainData4D(3));
    end
    
    %this array will hold all selection maps used for every single model 
    %in the different splits
    sizeData   = size(dataset.data4D);
    rfe_featureSelectionMaps     = zeros(sizeData(1),sizeData(2),sizeData(3),nmbSplits);
    rfe_weightMaps               = zeros(sizeData(1),sizeData(2),sizeData(3),nmbSplits);
    
  elseif(isfield(dataset, 'type') &&  strcmp(dataset.type,'dataset2D'))
    
    %2D case
    if(isfield(dataset,'featureSelectionMap1D') && ~isempty(dataset.featureSelectionMap1D))
      globalSelectionMap = dataset.featureSelectionMap1D;
    elseif(isfield(dataset,'mask1D') && ~isempty(dataset.mask1D))
      globalSelectionMap = dataset.mask1D;
    else
      globalSelectionMap = ones(sizeTrainData4D(1)*sizeTrainData4D(2)*sizeTrainData4D(3));
    end
    
    %this array will hold all selection maps used for every single model 
    %in the different splits
    sizeData   = size(dataset.data2D);
    rfe_featureSelectionMaps     = zeros(sizeData(1),nmbSplits);
    rfe_weightMaps               = zeros(sizeData(1),nmbSplits);
   
  else
    error('RFE: Please check the dataset: field "type" is not defined!');
  end
  
    
  totalNmbIn  = 0;
  totalNmbOut = 0;
   
  if(~easyupMVPA_getGlobals('quietMode'))
    disp(['Running Recursive Feature Elemination with command string: ',cmdString,' ...']);
  end
    
  
  %run over iterations (+1 is for the prediction of the last choice of features)
  for i = 1:nmbIterations+1
    if(~easyupMVPA_getGlobals('quietMode'))
      if(i<=nmbIterations)
        disp(['Started iteration ', num2str(i),' ...']);
        h = waitbar(0, ['Running Recursive Feature Elemination Iteration: ', num2str(i)]);
      else
        disp('Started final Prediction  ...');
        h = waitbar(0, 'Running Recursive Feature Elemination Final Prediction...');
      end 
    end
    
    %every iteration has its own tests over
    %nmbOfSplits examples
    nmbCorrect    = 0;
    nmbTruePos    = 0;
    nmbTrueNeg    = 0;
    nmbFalsePos   = 0;
    nmbFalseNeg   = 0;
    nmbTests      = 0; 

    %get the selected voxels of this iteration to show them
    actTotalNmbIn  = totalNmbIn;
    actTotalNmbOut = totalNmbOut;
    
    
    for j=1:nmbSplits
      
     if(~easyupMVPA_getGlobals('quietMode'))
       waitbar(j/nmbSplits,h);
     end
      
     %at first split the dataset according to given splitting
     ds1Indices = dataSplitter.splitMatrix(j,:) == 1; %test data
     ds2Indices = dataSplitter.splitMatrix(j,:) == 2; %train data
     
     [ds1, ds2] = splitDataset(dataset, ds1Indices, ds2Indices );     

     %the first iteration, that means there is no 
     %rfe selection map available and data can be restricted to
     %the globalSelectionMap
     if(i == 1)
       ds1 = setDataset_featureSelectionMap_ByMatrix(ds1, globalSelectionMap);
       ds2 = setDataset_featureSelectionMap_ByMatrix(ds2, globalSelectionMap);
       [ds1, testData2D]  = selectFeaturesBySelectionMap(ds1);
       [ds2, trainData2D] = selectFeaturesBySelectionMap(ds2);
     else
       %not the first iteration, we expect a rfe selection to be set for this split
       %so restrict data to this features
       if(strcmp(dataset.type,'dataset4D'))
         ds1 = setDataset_featureSelectionMap_ByMatrix(ds1, rfe_featureSelectionMaps(:,:,:,j));
         ds2 = setDataset_featureSelectionMap_ByMatrix(ds2, rfe_featureSelectionMaps(:,:,:,j));
       elseif(strcmp(dataset.type,'dataset2D'))
         ds1 = setDataset_featureSelectionMap_ByMatrix(ds1, rfe_featureSelectionMaps(:,j));
         ds2 = setDataset_featureSelectionMap_ByMatrix(ds2, rfe_featureSelectionMaps(:,j));
       end
       [ds1, testData2D]  = selectFeaturesBySelectionMap(ds1);
       [ds2, trainData2D] = selectFeaturesBySelectionMap(ds2);
     end
     
     %*** TRAINING ***
     svmModel  = svmtrain(double(ds2.classIDs)', trainData2D, cmdString);
     
     %*** TESTING ***
     %predict the class ID of the test data
     %check if the model was trained with probability estimates
     if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
       [predicted_label, accuracy, probEstimatesVector] = svmpredict(double(ds1.classIDs)', testData2D, svmModel, '-b 1');
       %disp(probEstimatesVector);
     else
       [predicted_label, accuracy, decisionValues] = svmpredict(double(ds1.classIDs)', testData2D, svmModel); 
       %disp(decisionValues);
     end
     
     %count correct predictions for statistics
     for k=1:size(predicted_label,1)
       if(predicted_label(k)==1) 
         if(ds1.classIDs(k) == 1)
           nmbCorrect = nmbCorrect+1;
           nmbTruePos = nmbTruePos+1;
         else
           nmbFalsePos = nmbFalsePos+1;
         end
       elseif(predicted_label(k)==0)
         if(ds1.classIDs(k) == 0)
           nmbCorrect = nmbCorrect+1;
           nmbTrueNeg = nmbTrueNeg+1;
         else
           nmbFalseNeg = nmbFalseNeg+1;
         end
       end
       nmbTests = nmbTests+1;
     end%endfor size predicted label
    
     
     %Do a new selection based on the weights 
     %but just for the nmb of iterations
     if(i <= nmbIterations)
       %*** SET THE SELECTION MAP ***
       if(strcmp(ds2.type,'dataset4D'))
%          %extract the weights
%          tmpWeights = zeros(sizeData(1),sizeData(2),sizeData(3),1);
%          weights      = svmModel.SVs' * svmModel.sv_coef;
%          tmpWeights(ds2.featureSelectionMap3D >0) = weights;
%          
%          %remove the features according to percentual selection threshold
%          [weightsRes, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(tmpWeights, thresholdPercentOfFeaturesOut);
% 
%          featSelMap = weightsRes;
%          featSelMap(weightsRes>0)=1;
%          featSelMap(weightsRes<0)=1;
%         
%          rfe_featureSelectionMaps(:,:,:,j) = featSelMap;
%          
%          %save the weight maps
%          rfe_weightMaps(:,:,:,j) = weightsRes;        


           %to extract the weights now a RFE is done
           %first select the actual training data

           %create a splitter that does an loocvNmbSplits x nmbSamples split
           nmbDS2Samples = length(ds2.classIDs);
           splitMatrix = ones(loocvNmbSplits,nmbDS2Samples);
           
           stepLength = floor(nmbDS2Samples/loocvNmbSplits);
           
           for k = 1:loocvNmbSplits
             startIndex = (k-1)*stepLength+1;
             splitMatrix(k,startIndex:startIndex+stepLength-1) = 2;
           end
           
           %a problem may arise if in one split not all class-types are in train data
           loocvCustomSplitter = getDataSplitter(ds2, 'custom', splitMatrix);
           
           [dstmp, resultStruct, loocvAvgWeights] = doLeaveOneOutCrossValidation_SVM(ds2, loocvCustomSplitter, kernelMode, costParam);
           
           
           %remove the features according to percentual selection threshold
           [weightsRes, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(loocvAvgWeights, thresholdPercentOfFeaturesOut);

           featSelMap = weightsRes;
           featSelMap(weightsRes>0)=1;
           featSelMap(weightsRes<0)=1;
           rfe_featureSelectionMaps(:,:,:,j) = featSelMap;
           
           
           showDataAsImage(featSelMap,'RFE - Nested LOOCV: Feature Selection Map');
           
           %save the weight maps
           rfe_weightMaps(:,:,:,j) = weightsRes;        
           
           
           %keyboard;
           
           
         
       elseif(strcmp(ds2.type,'dataset2D'))
         
%          %extract the weights
%          tmpWeights = zeros(sizeData(1),1);
%          weights      = svmModel.SVs' * svmModel.sv_coef;
%          tmpWeights(ds2.featureSelectionMap1D >0) = weights;
% 
%          %remove the features according to percentual selection threshold
%          [weightsRes, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(tmpWeights, thresholdPercentOfFeaturesOut);
% 
%          featSelMap = weightsRes;
%          featSelMap(weightsRes>0)=1;
%          featSelMap(weightsRes<0)=1;
%         
%          rfe_featureSelectionMaps(:,j) = featSelMap;
%          %save the weight maps
%          rfe_weightMaps(:,j) = weightsRes;        
         
       end
       
      end % endif nmb iteration
    end % endfor nmbSplits
    
    %For every iteration we get stat results
    accuracy    = nmbCorrect/nmbTests*100;
    sensitivity = nmbTruePos/(nmbTruePos+nmbFalseNeg); 
    specificity = nmbTrueNeg/(nmbTrueNeg+nmbFalsePos);
   
    resultStruct             = {};
    resultStruct.nmbTests    = nmbTests;
    resultStruct.accuracy    = accuracy;
    resultStruct.sensitivity = sensitivity;
    resultStruct.specificity = specificity;
    resultStruct.TP          = nmbTruePos;
    resultStruct.TN          = nmbTrueNeg;
    resultStruct.FP          = nmbFalsePos;
    resultStruct.FN          = nmbFalseNeg; 
    
    if(~easyupMVPA_getGlobals('quietMode'))
      close(h);
      if(i==1)
        disp(['RFE with ', num2str(nmbSplits),' Split(s). Selected all features for classification:']);
      else
        disp(['RFE with ', num2str(nmbSplits),' Split(s). Selected ' ,num2str(actTotalNmbOut),' out of ',num2str(actTotalNmbIn),' features for classification:']);
      end
      printResultStruct(resultStruct);
    end
    
  end %enfor iterations
    
  
  %THE LAST SUPPER
  %create the mean images and estimate deviation
  if(strcmp(dataset.type,'dataset4D'))
    avg_rfe_featureSelectionMap = zeros(sizeData(1),sizeData(2),sizeData(3),1);
    avg_rfe_weightMap           = zeros(sizeData(1),sizeData(2),sizeData(3),1);

    for i=1:nmbSplits
      avg_rfe_featureSelectionMap = avg_rfe_featureSelectionMap + rfe_featureSelectionMaps(:,:,:,i);
      avg_rfe_weightMap           = avg_rfe_weightMap + rfe_weightMaps(:,:,:,i);
    end
  elseif(strcmp(dataset.type,'dataset2D'))
    avg_rfe_featureSelectionMap = zeros(sizeData(1),1);
    avg_rfe_weightMap           = zeros(sizeData(1),1);

    for i=1:nmbSplits
      avg_rfe_featureSelectionMap = avg_rfe_featureSelectionMap + rfe_featureSelectionMaps(:,i);
      avg_rfe_weightMap           = avg_rfe_weightMap + rfe_weightMaps(:,i);
    end
  end
  
  avg_rfe_featureSelectionMap = avg_rfe_featureSelectionMap/nmbSplits;
  avg_rfe_weightMap           = avg_rfe_weightMap/nmbSplits;
  
end %end function
 




%member function to select a certain amount of elements from a 3D map
function [mapOut, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(mapIn, thresholdPercentOfElementsOut)

  mapOut = mapIn;
  
  %overall nmb of non-zero-elements
  totalNmbIn = sum(mapIn(:)~=0);
  
  %avgWeights3D = avgWeights3D + tmpWeights3D*(1/nmbSamples);
  %now set the rfe selection for this split to the according features
  %nmb of voxels to throw out
  nmbOut = floor(totalNmbIn*thresholdPercentOfElementsOut/100);
  
  %a sorted array to get the value for selection
  tmp = sort(abs(mapIn(mapIn(:)~=0)));
  selThresh = tmp(nmbOut);
  
  tmpMap = abs(mapIn);
  tmpMap(tmpMap<selThresh) = 0;
  mapOut(tmpMap==0) = 0;
  
  
  totalNmbOut = sum(mapOut(:)~=0);
  
end %end member_getSelectedMapByThresholdPercentOut

