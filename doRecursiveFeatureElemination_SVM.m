% Recursively removes features according to SVM classification weights (parallel execution).
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%   
%   [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut, dataSplitter, kernelMode, costParam, [paramStruct])
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
%   If possible, the function executes parallelized depending on settings defined via easyupMVPA_init().
%
% Parameters:
%   dataset                       - the datset to set the classIDs for
%   nmbIterations                 - iterations (how often is basic feature set reduced by "thresholdPercentOfFeaturesOut")
%   thresholdPercentOfFeaturesOut - percentual value of non-zero elements in weight map that should be cut (suggestion btw. 10 and 50)
%   dataSplitter                  - describes the splitting of the data in the background LOOCV
%   svmType                       - Types: Just 'classification' supported!
%   kernelMode                    - Kernels: ['linear', 'polynomial', 'radial', 'sigmoid']
%   costParam                     - The slack variable C in SVM (range 0 to 1  0 = low cost, 1 = highest costs). 
%                                   It defines the costs for misclassification (How strongly are outliers punished?).
%   paramStruct                   - example: {"degree", 3, "gamma", 0.5, "coeff0", 0.5, "nu", 0.3, "epsilon", 0.01, "probEstimates", 1}
%                                   possible fields:
%                                   'degree'        : default=3     Describes the exponent in polynomial kernel function
%                                   'gamma'         : default=1/k   The gamma-factor in kernel function
%                                   'coef0'         : default=0     The coefficient summand in kernel function
%                                   'nu'            : default=0.5   The nu parameter of regression_nu
%                                   'epsilon'       : default=0.001 The epsilon parameter of regression_epsilon
%                                   'probEstimates' : default=0     1 if probabilistic estimates should be computed, 0 if not%
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
function [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut, dataSplitter, svmType, kernelMode, costParam, paramStruct)
  
  if( ~exist('dataset','var') || ~exist('nmbIterations','var') || ~exist('thresholdPercentOfFeaturesOut','var') || ~exist('dataSplitter','var') || ~exist('svmType','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, svmType - [classification], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %extractt the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, paramStruct);
  end
  if( ~paramStructIsValid)
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, svmType - [classification], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  if(~strcmp(svmType, 'classification'))
    error('Sorry, the RFE function just supports C-classification in this version.');
  end
  
  %use quiet mode (no outputs)
  cmdString = [cmdString, ' -q '];
  
  if(thresholdPercentOfFeaturesOut < 0 || thresholdPercentOfFeaturesOut > 100)
    error('Usage of doRecursiveFeatureElemination_SVM: [dataset] = doRecursiveFeatureElemination_SVM(dataset, nmbIterations, thresholdPercentOfFeaturesOut - [0-100%], dataSplitter, svmType - [classification], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
      
  localQuietMode = easyupMVPA_getGlobals('quietMode');
  
  if(~localQuietMode)
    disp('Running Recursive Feature Elemination (This may take a while!) ...');
  end
  
  %extract the number of splits that are used
  nmbSplits = size(dataSplitter.splitMatrix,1);
  
  resultStruct   = {};
  
  %check if there is a global feature selection map is given, if yes all 
  %the following steps are computed solely on these selected features
  if(dataset.is4D)
    
    sizeData   = size(dataset.data);
    
    %4D case
    if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
      globalSelectionMap = dataset.featureSelectionMap;
    elseif(isfield(dataset,'mask') && ~isempty(dataset.mask))
      globalSelectionMap = dataset.mask;
    else
      globalSelectionMap = ones(sizeData(1)*sizeData(2)*sizeData(3));
    end
    
    %this array will hold all selection maps used for every single model 
    %in the different splits
    rfe_featureSelectionMaps     = zeros(sizeData(1),sizeData(2),sizeData(3),nmbSplits);
    rfe_weightMaps               = zeros(sizeData(1),sizeData(2),sizeData(3),nmbSplits);
    
    
    %if nmb iterations is higher than number of features set it to nmbFeatures-1
    if(nmbIterations >= sizeData(1)*sizeData(2)*sizeData(3))
      warning('Number of iterations is higher than number of features. Setting to nmbFeatures - 1 !');
      nmbIterations = sizeData(1)*sizeData(2)*sizeData(3) - 1;
    end
    
  elseif(dataset.is2D)
    
    %2D case
    if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
      globalSelectionMap = dataset.featureSelectionMap;
    elseif(isfield(dataset,'mask') && ~isempty(dataset.mask))
      globalSelectionMap = dataset.mask;
    else
      globalSelectionMap = ones(size(dataset.data,1),1);
    end
    
    %this array will hold all selection maps used for every single model 
    %in the different splits
    sizeData   = size(dataset.data);
    rfe_featureSelectionMaps     = zeros(sizeData(1),nmbSplits);
    rfe_weightMaps               = zeros(sizeData(1),nmbSplits);
  
    %if nmb iterations is higher than number of features set it to nmbFeatures-1
    if(nmbIterations >= sizeData(1))
      warning('Number of iterations is higher than number of features. Setting to nmbFeatures - 1 !');
      nmbIterations = sizeData(1) - 1;
    end
  else
    error('RFE: Please check the dataset: field "type" is not defined!');
  end
        
  if(~localQuietMode)
    disp(['Running Recursive Feature Elemination with command string: ',cmdString,' ...']);
  end
    
  
  %get the selected voxels of this iteration to show them
  actTotalNmbIn  = zeros(1, nmbSplits, 'int16');
  actTotalNmbOut = zeros(1, nmbSplits, 'int16');
  
  
  %run over iterations (+1 is for the prediction of the last choice of features)
  for i = 1:nmbIterations+1
    
   if(~localQuietMode)
      if(i<=nmbIterations)
        disp(['Started iteration ', num2str(i),' ...']);
        h = waitbar((i-1)/(nmbIterations+1), ['Running Recursive Feature Elemination Iteration: ', num2str(i)]);
      else
        disp('Started final Prediction  ...');
        h = waitbar((i-1)/(nmbIterations+1), 'Running Recursive Feature Elemination Final Prediction...');
      end 
    end     
    
    %every iteration has its own tests over
    %nmbOfSplits examples
    nmbCorrect     = 0;
    nmbTests       = 0; 
    nmbTruePosAll  = 0;
    nmbTrueNegAll  = 0;
    nmbFalsePosAll = 0;
    nmbFalseNegAll = 0; 
   
    
    
    %distinct parfors for 4D and 2D because of parallelization
    if(dataset.is4D)
      
      %******** 4D Case ************
      if(~localQuietMode)
        disp(['Running Leave One Out Cross Validation with command string: ',cmdString,' ...']);
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
      
       parfor j=1:nmbSplits
               
         if(~localQuietMode)
           if(progressIndices(j)<progressIndices(j+1))
             fprintf('\b\b*');
             disp(['' 0]);
           end
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
             ds1 = setDataset_featureSelectionMap_ByMatrix(ds1, rfe_featureSelectionMaps(:,:,:,j));
             ds2 = setDataset_featureSelectionMap_ByMatrix(ds2, rfe_featureSelectionMaps(:,:,:,j));
             [ds1, testData2D]  = selectFeaturesBySelectionMap(ds1);
             [ds2, trainData2D] = selectFeaturesBySelectionMap(ds2);
         end

         %*** TRAINING ***
         svmModel  = svmtrain(double(ds2.classIDs)', trainData2D, cmdString);

         %*** TESTING ***
         %predict the class ID of the test data
         %check if the model was trained with probability estimates
         if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
           [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', testData2D, svmModel, '-b 1');
           %disp(probEstimatesVector);
         else
           [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', testData2D, svmModel); 
           %disp(decisionValues);
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
         nmbCorrect  = nmbCorrect+nmbTruePos+nmbTrueNeg;

         nmbTruePosAll    = nmbTruePosAll+nmbTruePos;
         nmbTrueNegAll    = nmbTrueNegAll+nmbTrueNeg;
         nmbFalsePosAll   = nmbFalsePosAll+nmbFalsePos;
         nmbFalseNegAll   = nmbFalseNegAll+nmbFalseNeg; 
         nmbTests         = nmbTests+length(ds1.classIDs);
 
         %Do a new selection based on the weights 
         %but just for the nmb of iterations
         if(i <= nmbIterations)
           %*** SET THE SELECTION MAP ***
             %extract the weights
             tmpWeights = zeros(sizeData(1),sizeData(2),sizeData(3),1);
             weights      = svmModel.SVs' * svmModel.sv_coef;
             tmpWeights(ds2.featureSelectionMap >0) = weights;

             %remove the features according to percentual selection threshold
             [weightsRes, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(tmpWeights, thresholdPercentOfFeaturesOut);

             actTotalNmbIn(j)  = totalNmbIn;
             actTotalNmbOut(j) = totalNmbOut;
             
             featSelMap = weightsRes;
             featSelMap(weightsRes>0)=1;
             featSelMap(weightsRes<0)=1;

             rfe_featureSelectionMaps(:,:,:,j) = featSelMap;

             %save the weight maps
             rfe_weightMaps(:,:,:,j) = weightsRes;        

          end % endif nmb iteration
          
       end % endfor nmbSplits (parfor)
    
       if(~localQuietMode)
         fprintf('\n');
       end
       
    
    else  %******** 2D Case ************
        
      if(~localQuietMode)
        disp(['Running Leave One Out Cross Validation with command string: ',cmdString,' ...']);
        disp(['0%', num2str(repmat(' ',1,nmbSplits)),'100%']);
        disp('   ');
      end
      
       parfor j=1:nmbSplits
               
         if(~localQuietMode)           
           fprintf('\b\b*');
           disp(['' 0]);
         end
         
         %at first split the dataset according to given 
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
           ds1 = setDataset_featureSelectionMap_ByMatrix(ds1, rfe_featureSelectionMaps(:,j));
           ds2 = setDataset_featureSelectionMap_ByMatrix(ds2, rfe_featureSelectionMaps(:,j));
           [ds1, testData2D]  = selectFeaturesBySelectionMap(ds1);
           [ds2, trainData2D] = selectFeaturesBySelectionMap(ds2);
         end

         %*** TRAINING ***
         svmModel  = svmtrain(double(ds2.classIDs)', trainData2D, cmdString);

         %*** TESTING ***
         %predict the class ID of the test data
         %check if the model was trained with probability estimates
         if(~isempty(svmModel.ProbA) && ~isempty(svmModel.ProbB))
           [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', testData2D, svmModel, '-b 1');
           %disp(probEstimatesVector);
         else
           [predicted_labels, accuracy, probEst] = svmpredict(double(ds1.classIDs)', testData2D, svmModel); 
           %disp(decisionValues);
         end

         %count correct predictions for statistics
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
         nmbCorrect  = nmbCorrect+nmbTruePos+nmbTrueNeg;

         nmbTruePosAll    = nmbTruePosAll+nmbTruePos;
         nmbTrueNegAll    = nmbTrueNegAll+nmbTrueNeg;
         nmbFalsePosAll   = nmbFalsePosAll+nmbFalsePos;
         nmbFalseNegAll   = nmbFalseNegAll+nmbFalseNeg; 
         nmbTests         = nmbTests+length(ds1.classIDs);


         %Do a new selection based on the weights 
         %but just for the nmb of iterations
         if(i <= nmbIterations)
           %*** SET THE SELECTION MAP ***
           
             %extract the weights
             tmpWeights = zeros(sizeData(1),1);
             weights      = svmModel.SVs' * svmModel.sv_coef;
             tmpWeights(ds2.featureSelectionMap >0) = weights;

             %remove the features according to percentual selection threshold
             [weightsRes, totalNmbIn, totalNmbOut] = member_getSelectedMapByThresholdPercentOut(tmpWeights, thresholdPercentOfFeaturesOut);

             actTotalNmbIn(j)  = totalNmbIn;
             actTotalNmbOut(j) = totalNmbOut;

             featSelMap = weightsRes;
             featSelMap(weightsRes>0)=1;
             featSelMap(weightsRes<0)=1;

             rfe_featureSelectionMaps(:,j) = featSelMap;
             %save the weight maps
             rfe_weightMaps(:,j) = weightsRes;        

          end % endif nmb iteration

       end % endfor nmbSplits
       
       if(~localQuietMode)
           fprintf('\n');
       end

       
    end %endif 4D/2D case
        
    
    %For every iteration we get stat results
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
    
    %if(isfield(resultStruct, 'innerResultStruct') && ~isempty(resultStruct.innerResultStruct))
    if(~isfield(resultStruct, 'innerResultStruct') || isempty(resultStruct.innerResultStruct) || accuracy >= resultStruct.innerResultStruct.accuracy)
      resultStruct.innerResultStruct.nmbTests    = nmbTests;
      resultStruct.innerResultStruct.accuracy    = accuracy;
      resultStruct.innerResultStruct.sensitivity = sensitivity;
      resultStruct.innerResultStruct.specificity = specificity;
      resultStruct.innerResultStruct.TP          = nmbTruePosAll;
      resultStruct.innerResultStruct.TN          = nmbTrueNegAll;
      resultStruct.innerResultStruct.FP          = nmbFalsePosAll;
      resultStruct.innerResultStruct.FN          = nmbFalseNegAll;
      resultStruct.innerResultStruct.infoString  = ['This struct holds RFE maximum performance results (iteration ',num2str(i),' used nmb of features: ',num2str(actTotalNmbOut(end)),').'];
    end
    
    if(isfield(resultStruct, 'infoString'))
      resultStruct.infoString  = [resultStruct.infoString, ' * ', 'RFE:Iteration ', num2str(i), ' Accuracy: ', num2str(accuracy), '% NmbFeatures: ', num2str(actTotalNmbOut(end))];
    else
      resultStruct.infoString  = ['RFE:Iteration ', num2str(i), 'Accuracy: ', num2str(accuracy), '% NmbFeatures: ', num2str(actTotalNmbOut(end))];
    end
    
    
    if(~easyupMVPA_getGlobals('quietMode'))
      close(h);
      if(i==1)
        disp(['RFE with ', num2str(nmbSplits),' Split(s). Selected all features for classification:']);
      else
        disp(['RFE with ', num2str(nmbSplits),' Split(s). Selected ' ,num2str(actTotalNmbOut(end)),' out of ',num2str(actTotalNmbIn(end)),' features for classification:']);
      end
      printResultStruct(resultStruct);
    end
    
  end %enfor iterations
    
  
  %THE LAST SUPPER
  %create the mean images and estimate deviation
  if(dataset.is4D)
    avg_rfe_featureSelectionMap = zeros(sizeData(1),sizeData(2),sizeData(3),1);
    avg_rfe_weightMap           = zeros(sizeData(1),sizeData(2),sizeData(3),1);

    parfor i=1:nmbSplits
      avg_rfe_featureSelectionMap = avg_rfe_featureSelectionMap + rfe_featureSelectionMaps(:,:,:,i);
      avg_rfe_weightMap           = avg_rfe_weightMap + rfe_weightMaps(:,:,:,i);
    end
    
    %set the last rfe selection map for returned dataset
    dataset = setDataset_featureSelectionMap_ByMatrix(dataset, rfe_featureSelectionMaps(:,:,:,end));
         
  elseif(dataset.is2D)
    avg_rfe_featureSelectionMap = zeros(sizeData(1),1);
    avg_rfe_weightMap           = zeros(sizeData(1),1);

    parfor i=1:nmbSplits
      avg_rfe_featureSelectionMap = avg_rfe_featureSelectionMap + rfe_featureSelectionMaps(:,i);
      avg_rfe_weightMap           = avg_rfe_weightMap + rfe_weightMaps(:,i);
    end
    
    %set the last rfe selection map for returned dataset
    dataset = setDataset_featureSelectionMap_ByMatrix(dataset, rfe_featureSelectionMaps(:,end));  
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

