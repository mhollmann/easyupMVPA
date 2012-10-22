% !UNDER CONSTRUCTION! Implements a Searchlight Classification using One Out Cross Validation (LOOCV) with a SVM Classifier (parallel execution).
%
% Author: Maurice Hollmann
% Date  : 05/11
%
% Description:
%
%   [dataset, resultStruct] = doLeaveOneOutCrossValidation_SVM_2DforceQuiet(dataset, dataSplitter, svmCommandString)
%
%   Just 4D dataset input!
%   It is recommended to use a mask including just grey matter to reduce computation time!
%
%
% Parameters:
%   dataset          - The dataset to work on  (all samples are included in LOOCV)
%   dataSplitter     - describes the splitting of the data in LOOCV
%   svmCommandString - the svm command string for exmaple '-t 0 -c 0.5' 
%
% Returns:
%   dataset            - the datset that has been the input
%   resultAccuracyMap  - a 3D map with the results of searchlight classification for each voxel
%   resultStruct       - The struct holding the classification results: 
%                        resultStruct.nmbTests     (the number of samples tested for this result)
%                        resultStruct.accuracy     (percentual value of correct predictions (correct * 100 / nmbSamples))
%                        resultStruct.sensitivity  (TP/TP+FN = Proportion of true positives to all positives)
%                        resultStruct.specificity  (TN/TN+FP = Proportion of true negatives to all negatives)
%                        resultStruct.TP           (True positives = all correct predicted in class 1)
%                        resultStruct.TN           (True negatives = all correct predicted in class 2)
%                        resultStruct.FP           (False positives = all incorrect predicted in class 1)
%                        resultStruct.FN           (False negatives = all incorrect predicted in class 2)
%
% Comments:
%
function [datasetOut, resultAccuracyMap, resultStruct] = doSearchlightLOOCV_SVM(dataset, searchlightDiameter, dataSplitter, kernelMode, costParam, paramStruct)

  
  %extractt the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, paramStruct);
  end
  if( ~paramStructIsValid)
    error(['Usage of doLeaveOneOutCrossValidation_SVM: [datasetOut, resultAccuracyMap3D, resultStruct] = doSearchlightLOOCV_SVM(dataset, searchlightDiameter, dataSplitter,', ... 
          'kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])']);
  end
  
  tic
  
  %for each voxel extract all voxels in the searchlight volume
  %and do a LOOCV on these...
  datasetOut = [];
  resultStruct = [];
  
  if(dataset.is4D)
    sizeData   = size(dataset.data);
    %4D case
    if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
      globalSelectionMap = dataset.featureSelectionMap;
    elseif(isfield(dataset,'mask') && ~isempty(dataset.mask))
      globalSelectionMap = dataset.mask;
    else
      globalSelectionMap = ones(sizeData(1),sizeData(2),sizeData(3));
    end

      %get the adjacency matrix
      adjMatrix = getSearchlightAdjacencyMatrix(dataset, searchlightDiameter);
      sl_numEl  = length(adjMatrix(adjMatrix>0));
      adjMatrix = repmat(adjMatrix, [1, 1, 1, sizeData(4)]);
  elseif(dataset.is2D)
  end
  
  %now move through all voxels in feature selection mask
  indVec = find(globalSelectionMap>0);
  
  %performanceMatrices = zeros(sizeData(1),sizeData(2),sizeData(3));
  searchlightRadius = floor(searchlightDiameter/2);
  dataDim4 = sizeData(4);
  
  %use as 1D because of parfor
  resultAccuracyVec = zeros(1,length(indVec));
   
  
   localQuietMode = easyupMVPA_getGlobals('quietMode');
   
   if(~localQuietMode)
     disp(['Running Searchlight Leave One Out Cross Validation with command string: ',cmdString,' ...']);
     
     % create a progress display that works also for parallel loops 
     if(length(indVec) <50)
       disp(['0%', num2str(repmat(' ',1,length(indVec)-1)),'100%']);
       progressIndices = [1 1:length(indVec)];
     else
       disp(['0%', num2str(repmat(' ',1,50)),'100%']);
       %create a vector with floored indicees for repetitions
       progressIndices = [1 1:length(indVec)];
       progressIndices = floor(progressIndices*(50/length(indVec)));
     end
     fprintf('    ');   
   end
  
  
  
  
  
  
  %disp(length(indVec));
  %for i=1:length(indVec)
  parfor i=1:length(indVec)
    
    %disp(i)
    
    
    data3D = [];

    if(~localQuietMode)
      if(progressIndices(i)<progressIndices(i+1))
        fprintf('\b\b*');
        disp(['' 0]);
      end
    end
    
    %place the searchlight
    [indX,indY,indZ] = ind2sub(size(globalSelectionMap), indVec(i));
    
    %try because of possible invalid indicees
    try
      %extract the data (at first the whole cubic volume)
      data3D = dataset.data(indX-searchlightRadius:indX+searchlightRadius,...
                            indY-searchlightRadius:indY+searchlightRadius,...
                            indZ-searchlightRadius:indZ+searchlightRadius,:);
    catch
      continue;
    end
    
    %data2D = reshape(data3D(adjMatrix>0), sl_numEl, dataDim4);
    
    ds2D = getEmpty2DDataset();
    ds2D.data       = reshape(data3D(adjMatrix>0), sl_numEl, dataDim4);
    ds2D.chunks     = dataset.chunks;
    ds2D.classIDs   = dataset.classIDs;
    
    %printDatasetInfo(ds2D);
    %tic
    %disp(size(data2D));
    [ds, sl_resStruct] = doLeaveOneOutCrossValidation_SVM_2DforceQuiet(ds2D, dataSplitter, cmdString);
    %toc
    %disp(sl_resStruct.accuracy);
    
    resultAccuracyVec(i) = sl_resStruct.accuracy;
    
  end
  
  toc
  
  resultAccuracyMap = zeros(sizeData(1), sizeData(2), sizeData(3));
  resultAccuracyMap(find(globalSelectionMap>0)) = resultAccuracyVec;  
  
  
end