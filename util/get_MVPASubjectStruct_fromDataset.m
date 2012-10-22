% Returns the Princeton MVPA subject struct from datasetInput.
% IMPORTANT: The path must contain the princeton mvpa toolbox.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%
%   This method assumes that samples coded with 0 in chunks are NOT part of any class (will be set to 0 in mvpas-regs struct).
%   The chunk numbers will be translated into run-numbers!
%
% Parameters:
%
% Returns:
%
% Comments:
%
function [mvpaSubjectStruct] = get_MVPASubjectStruct_fromDataset(datasetIn, experimentID, datasetID)
  
  mvpaSubjectStruct = struct();

  if(~datasetIn.is4D)
    error('Sorry conversion just usable for 4D datasets...');
    return
  end

  %subject is the basic structure
  mvpaSubjectStruct = init_subj(experimentID,datasetID);
  
  
  %*** setting mask *** 
  if(~isempty(datasetIn.mask))
    mvpaSubjectStruct = init_object(mvpaSubjectStruct,'mask','dataMask'); 
    mvpaSubjectStruct = set_mat(mvpaSubjectStruct,'mask', 'dataMask', datasetIn.mask);
  end 
  
  
  %*** set the image data as pattern ***
  %all features are used to create a 2D pattern
  sizeData = size(datasetIn.data);
  if(datasetIn.is2D)
   datasetIn = setDataset_featureSelectionMap_ByMatrix(datasetIn, ones(sizeData(1:end-1),1));
  else
   datasetIn = setDataset_featureSelectionMap_ByMatrix(datasetIn, ones(sizeData(1:end-1))); 
  end
  [datasetIn, data2D] = selectFeaturesBySelectionMap(datasetIn);

  %we have to init the object by ourselves if we want to use nifti
  mvpaSubjectStruct = init_object(mvpaSubjectStruct,'pattern','data2D');
  mvpaSubjectStruct = set_mat(mvpaSubjectStruct,'pattern', 'data2D',  data2D');
  mvpaSubjectStruct = set_objfield(mvpaSubjectStruct, 'pattern', 'data2D',  'masked_by', 'dataMask');
  
  
  
  %*** set regressors and conditions ***
  %regs : unique classNumbers refer to own vector in x dimension (1==belongs to regressor, 0 otherwise)
  uniqueClasses = unique(datasetIn.classIDs);
  
  regs = zeros(length(uniqueClasses), length(datasetIn.classIDs));
  runs = datasetIn.chunks;
  cndNamesString = '{';
  for i=1:length(uniqueClasses)
    
    tmpReg = zeros(1, length(datasetIn.classIDs));
    
    tmpReg(datasetIn.classIDs==uniqueClasses(i) & datasetIn.chunks > 0) = 1;
    regs(i,:) = tmpReg;
    cndNamesString = [cndNamesString, '''class',num2str(i), ''','];   
  end
  
  cndNamesString = [cndNamesString(1:end-1),'};'];
  disp(cndNamesString);
  
  %regressors and conditions 
  mvpaSubjectStruct      = init_object(mvpaSubjectStruct,'regressors','conditions');
  mvpaSubjectStruct      = set_mat(mvpaSubjectStruct,'regressors','conditions',regs);
  
  mvpaSubjectStruct      = set_objfield(mvpaSubjectStruct,'regressors','conditions','condnames',eval(cndNamesString));
  
  %a selector over the runs (runs is a vector that holds an integer for every timepoint which is the run nmb)
  mvpaSubjectStruct = init_object(mvpaSubjectStruct,'selector','runs'); 
  mvpaSubjectStruct = set_mat(mvpaSubjectStruct,'selector','runs', runs);
  
  
  disp('Info Input easyupMVPA-dataset:');
  printDatasetInfo(datasetIn);
  disp(' ');
  disp('Info Out matMVPA-dataset:');
  summarize(mvpaSubjectStruct);
  

end