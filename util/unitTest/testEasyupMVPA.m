% Unit Test for basic functions of the toolbox for developement.
%
% Author: Maurice Hollmann
% Date  : 02/11
%
% Description:
function [messageStack] = testEasyupMVPA()

  tic;
  
  %supress graphics output?
  graphOut = false;

  %set scopes to check to true:
  scopesToCheck = struct(...
    'check_dataCreation',            1,...
    'check_preprocessing',           0,...
    'check_sampleSelection',         0,...
    'check_simpleClassificationSVM', 1,...
    'check_simpleClassificationRVM', 0,...
    'check_LOOCV_SVM',               1,...
    'check_RFE_SVM',                 0,...
    'check_RFE_BOOTSTRAP_SVM',       0,...
    'check_Searchlight_SVM',         0,...
    'check_Configurator',            0,...
    'check_regression',              1);
  
  
  messageStack = struct('testScope', '', 'testMessage', '');
  msgIndex     = 1;
  errorCount = 0;
  fs = filesep();
  testFuncPath = fileparts(which('testEasyupMVPA'));
  dir_data = [testFuncPath, fs, 'testData', fs];
  easyupMVPA_init('nmbCores', 2);
    
  
  wPath = which('easyupMVPA_init');
  [pathstr] = fileparts(wPath);
  addpath(genpath(pathstr));
  
  
  
  %% Scope preparation - always excecuted
  
  %creates variable datasetTest
  load([dir_data, 'datasetTest.mat']);
  sizeDataset = size(datasetTest.data);
  
  grayMatterMaskFile = [dir_data, 'grayMatter_mask.nii'];
  grayMatterMaskNii  = load_untouch_nii(grayMatterMaskFile);
  grayMatterMask     = grayMatterMaskNii.img;
  
  if(graphOut)
    showDataAsImage(grayMatterMask, 'gm_mask');
  end
  
  %mask files to use
  rois = {'lDLPFC_post_glm.nii', 'rDLPFC_post_glm.nii', 'lCaudateHead_glm.nii', 'rCaudateHead_glm.nii', ...
           'lIFG_AI_glm.nii', 'rIFG_AI_glm.nii', 'lOFC_glm.nii', 'rOFC_glm.nii',  ...
           'lTPJ_glm.nii', 'rTPJ_glm.nii', 'mPFC_preSMA_glm.nii'};
  
  %concat all rois chosen
  for i=1:size(rois,2)
    
    imgFile = [dir_data, char(rois{i})];
    disp(imgFile);
    dataNii = load_untouch_nii(imgFile);
    sizeIMG = size(dataNii.img);
    if(~exist('roi3D','var'))
      roi3D = zeros(sizeIMG(1), sizeIMG(2), sizeIMG(3), 'int16');
    end
    dataNii.img = int16(dataNii.img);
    roi3D(dataNii.img~=0) = i;
  end
  
  disp('Done loading rois ...');
  if(graphOut)
    showDataAsImage(roi3D, 'roi3D');
  end
  
  [datasetTest, datasetTest2D] = averageFeaturesInROIs(datasetTest, roi3D);
  
  messageStack(msgIndex).testScope   = 'preparation';
  messageStack(msgIndex).testMessage = 'feature averaging features in dataset 4D ... PASSED';
  msgIndex =  msgIndex + 1;
  
  
   
  
  
  
%% Scope dataCreation
  
if(scopesToCheck.check_dataCreation || scopesToCheck.check_preprocessing)  
  %**** CASE nifti 4D ****
  %get the attributes from a txt file
  attribFile = [dir_data,'attributes_testData.txt'];
  %fileList will be a character array in this case
  fileList = spm_get('Files', dir_data, 'nifti_*.hdr');
  %create an empty dataset
  dataset = getEmpty4DDataset();
  %set the field data4D of the dataset
  dataset = setDataset_data_ByFilelist(dataset, fileList);
  %set the chunks 
  dataset = setDataset_chunks_ByAttribFile(dataset, attribFile);
  %set the classIDs
  dataset = setDataset_classIDs_ByAttribFile(dataset, attribFile);
  %set a brain mask
  dataset = setDataset_mask_ByImageFile(dataset, [dir_data,'brainMask_53_63_46.hdr']);
  printDatasetInfo(dataset);
  
  %show the image and the mask
  if(graphOut)
    showDataAsImage(dataset.data, 'First timepoint in data', 1);
    showDataAsImage(dataset.mask, 'Dataset mask');
  end
  
  %check this scope for 4D
  sizeData = size(dataset.data);
  ok =true;
  if(sizeData ~= [53 63 46 507])
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'ERROR: creation 4D by nifti images ... FAILED! Dimension mismatch!';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
    ok = false;
  end
  
  if(~checkDataset(dataset))
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'ERROR: creation 4D by nifti images ... FAILED! Unspecified error while checking dataset!';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
    ok = false;
  end
  
  if(ok)
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'creation 4D by nifti images ... PASSED';
    msgIndex =  msgIndex + 1;
  end

  %****END CASE nifti 4D ****
  
  
  %**** CASE nifti 2D ****
  
  
  %**** END CASE nifti 2D ****
  
  
  
    
  %*** Sample selection in 4D datasets *** 
  ds1 = selectSamples(datasetTest, 'index==1');
  ds2 = selectSamples(datasetTest, 'index>19');
  selVec = zeros(1,40);
  selVec(20:40) = 1;
  ds3 = selectSamples(datasetTest, selVec);
  ds4 = selectSamples(datasetTest, 'classIDs<1');
  scopeSelOK = true;
  if( ~(size(ds1.data,4) == 1) || ~(length(ds1.classIDs)==1) || ~(length(ds1.chunks)==1))
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'ERROR: error selecting one sample in data 4D !';
    msgIndex =  msgIndex + 1;
    scopeSelOK = false;
  end
  
  if( ~sum((ds2.chunks == ds3.chunks))==20 || ~sum((ds2.classIDs == ds3.classIDs))==20 || ~sum(size(ds2.data,4)==size(ds3.data,4))==20  )
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'ERROR: error selecting samples by indices or vectors in data 4D !';
    msgIndex =  msgIndex + 1;
    scopeSelOK = false;
  end
  
  if( ~(size(ds4.data,4) == 18) || ~(length(ds4.classIDs)==18) || ~(length(ds4.chunks)==18))
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'ERROR: error selecting sample by classIDs in data 4D !';
    msgIndex =  msgIndex + 1;
    scopeSelOK = false;
  end
  
  if(scopeSelOK)
    messageStack(msgIndex).testScope   = 'dataCreation';
    messageStack(msgIndex).testMessage = 'sample selection in data 4D ... PASSED';
    msgIndex =  msgIndex + 1;
  end

  clear('dataset','ds1', 'ds2', 'ds3', 'ds4');
  
end %endif check scope
  
   
%% Scope preprocessing
  
if(scopesToCheck.check_preprocessing)
  

  dataset = selectSamples(datasetTest, 'index<=30');  
  a = dataset.data(20,20,20,:);
  
  %linear detrending
  dataset = doLinearDetrending(dataset, [3 4]);
  printDatasetInfo(dataset);
  b = dataset.data(20,20,20,:);
  

  %highpass-filtering (Take a careful look at your design - it is possible that useful frequencies are removed!)
  %Usage: doHighpassFiltering(dataset, cutoffFreq, TR)
  %In this case just a very low freq of more than 300 seconds is the cutoff, because of trial design(look at classIDs)
  dataset = doHighpassFiltering(dataset, 0.5, 1/300);
 
  c = dataset.data(20,20,20,:);
  
  %Plot timecourse of a single voxel
  if(graphOut)
    figureID = 'easyup_unitTest_preprocFig';
    if(isempty(findobj('Tag', figureID)))    
        resFig = figure('Tag', figureID, 'Name', 'Results for preprocessing');
    else
        figure(findobj('Tag', figureID));
    end

    subplot(3,1,1); plot(a(:));
    subplot(3,1,2); plot(b(:));
    subplot(3,1,3); plot(c(:));
  end
  
  printDatasetInfo(dataset);
   
end

  
%% Scope: sampleSelection
  
if(scopesToCheck.check_sampleSelection)
  ok = true;
  ds = selectSamples(datasetTest, 'index==1');
  ds = selectSamples(datasetTest, ['index==', num2str(sizeDataset(4))]);
  
  if(ok)
    messageStack(msgIndex).testScope   = 'sampleSelection';
    messageStack(msgIndex).testMessage = 'sample selection dataset 4D ... PASSED';
    msgIndex =  msgIndex + 1;
  end
end %endif scope   

  
  
%% Scope: simpleClassificationSVM
if(scopesToCheck.check_simpleClassificationSVM)
  disp('********** SCOPE: simpleClassificationSVM ************');
  testScope = 'simpleClassificationSVM';
  
  %** classification in 2D **
  %single samples
  dsTrain = selectSamples(datasetTest2D, 'index~=1');
  dsTest  = selectSamples(datasetTest2D, 'index==1');
  [svmModel, weights] = train_SVM(dsTrain, 'classification', 'linear', 0.5);
  [resultStruct] = predict_SVM(dsTest, svmModel);
  printResultStruct(resultStruct);
  
  %more samples
  dsTrain = selectSamples(datasetTest2D, 'index<=20');
  dsTest  = selectSamples(datasetTest2D, 'index>20');
  [svmModel, weights] = train_SVM(dsTrain, 'classification', 'linear', 0.5);
  [resultStruct] = predict_SVM(dsTest, svmModel);
  printResultStruct(resultStruct);
  
  messageStack(msgIndex).testScope   = testScope;
  messageStack(msgIndex).testMessage = 'feature averaging features in dataset 4D ... PASSED';
  msgIndex =  msgIndex + 1;

  %keyboard;
  
  if(resultStruct.accuracy == 25)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'classification result dataset 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check classification result dataset 2D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end
  
  %** classification in 4D **
  %single samples
  dsTrain = selectSamples(datasetTest, 'index~=1');
  dsTest  = selectSamples(datasetTest, 'index==1');
  [svmModel, weights] = train_SVM(dsTrain, 'classification', 'linear', 0.5);
  [resultStruct] = predict_SVM(dsTest, svmModel);
  printResultStruct(resultStruct);
  
  dsTrain = selectSamples(datasetTest, 'index<=30');
  dsTest  = selectSamples(datasetTest, 'index>30');
  [svmModel, weights] = train_SVM(dsTrain, 'classification', 'linear', 0.5);
  [resultStruct] = predict_SVM(dsTest, svmModel);
  printResultStruct(resultStruct);
  
  if(resultStruct.accuracy == 80)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'classification result dataset 4D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check classification result dataset 4D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end
  
  
  %noMask
  dsTrain = selectSamples(datasetTest, 'index<=30');
  dsTest  = selectSamples(datasetTest, 'index>30');
  
  dsTrain.mask = [];
  dsTest.mask = [];
  [svmModel, weights] = train_SVM(dsTrain, 'classification', 'linear', 0.5);
  [resultStruct] = predict_SVM(dsTest, svmModel);
  printResultStruct(resultStruct);
  
  
end %endif scope
  
 



  %% Scope: simpleClassificationRVM
if(scopesToCheck.check_simpleClassificationRVM)  
  
  disp('********** SCOPE: simpleClassificationRVM ************');
  testScope = 'simpleClassificationRVM';
  
  
  %** classification in 4D **
  datasetTestRVM = setDataset_mask_ByMatrix(datasetTest, roi3D);
  dsTrain = selectSamples(datasetTest, 'index<=25');
  dsTest  = selectSamples(datasetTest, 'index>25');
  printDatasetInfo(datasetTestRVM);
  [rvmModel, weights] = train_RVM(dsTrain, '+poly1', 0.05);
  %[rvmModel, weights] = train_RVM(dsTrain, '+poly1', 0.000005);
  disp(rvmModel);
  [resultStruct] = predict_RVM(dsTest, rvmModel);
  printResultStruct(resultStruct);
  
  %return
  %keyboard
  
  %** classification in 2D **
  dsTrain = selectSamples(datasetTest2D, 'index<=33');
  dsTest  = selectSamples(datasetTest2D, 'index>33');
  [rvmModel, weights] = train_RVM(dsTrain, '+gauss', 0.00001);
  disp(rvmModel);
  [resultStruct] = predict_RVM(dsTest, rvmModel);
  printResultStruct(resultStruct);
  
  
%   [resultStruct] = predict_SVM(dsTest, svmModel);
%   printResultStruct(resultStruct);
  
%   messageStack(msgIndex).testScope   = testScope;
%   messageStack(msgIndex).testMessage = 'feature averaging features in dataset 4D ... PASSED';
%   msgIndex =  msgIndex + 1;
  
  
end %endif scope






%% Scope: LOOCV_SVM
if(scopesToCheck.check_LOOCV_SVM)    
  
  disp('********** SCOPE: LOOCV_SVM ************');
  testScope = 'LOOCV_SVM';
  
  splitterOE  = getDataSplitter(datasetTest, 'oddEvenSamples');
  splitterOSO = getDataSplitter(datasetTest, 'oneSampleOut');
  
  if(isequal(splitterOE.splitMatrix,repmat([1 2], 1, 20)) && isequal(size(splitterOSO.splitMatrix), [40 40]))
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'check splitter generation dataset 4D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check splitter generation dataset 4D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end
  
  
  
%   %LOOCV 4D no mask
%   datasetTestLOOCV = datasetTest;
%   datasetTestLOOCV.mask = [];
%   datasetTestLOOCV.featureSelectionMap = [];
%   
%   tic
%   [datasetX, resultStruct, avgWeights3D] = doLeaveOneOutCrossValidation_SVM(datasetTestLOOCV, splitterOSO, 'linear', 0.5);
%   tN = toc;
%   printResultStruct(resultStruct);
  
  
  
  %LOOCV 4D including mask
  datasetTestLOOCV = setDataset_mask_ByMatrix(datasetTest, roi3D);
  tic
  [datasetX, resultStruct, avgWeights3D] = doLeaveOneOutCrossValidation_SVM(datasetTestLOOCV, splitterOSO, 'classification', 'linear', 0.5);
  tN = toc;
  printResultStruct(resultStruct);
    
  if(resultStruct.accuracy == 80)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = ['LOOCV classification result dataset 4D (Time needed:',num2str(tN),')... PASSED'];
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check LOOCV classification result dataset 4D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end

  
  %splitter 2D
  splitterOE  = getDataSplitter(datasetTest2D, 'oddEvenSamples');
  splitterOSO = getDataSplitter(datasetTest2D, 'oneSampleOut');
  
  if(isequal(splitterOE.splitMatrix,repmat([1 2], 1, 20)) && isequal(size(splitterOSO.splitMatrix), [40 40]))
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'check splitter generation dataset 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check splitter generation dataset 2D ... FAILED';
    msgIndex =  msgIndex + 1; 
    errorCount = errorCount +1;
  end
  
  
  %LOOCV 2D
  [datasetTest2DX, resultStruct, avgWeights1D] = doLeaveOneOutCrossValidation_SVM(datasetTest2D, splitterOSO, 'classification', 'linear', 0.5);
  printResultStruct(resultStruct);
  
  if(single(resultStruct.accuracy) == 55)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'LOOCV classification result dataset 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check LOOCV classification result dataset 2D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end
  
  
  
  %LOOCV 2D with masked 2D
  tDS = datasetTest2D;
  tDS  = setDataset_mask_ByMatrix(tDS, [0 1 1 1 0 0 0 0 1 1 1 ]);
  printDatasetInfo(tDS);
   
  [datasetTest2DX, resultStruct, avgWeights1D] = doLeaveOneOutCrossValidation_SVM(tDS, splitterOSO, 'classification', 'linear', 0.5);
  printResultStruct(resultStruct);
 
  if(single(resultStruct.accuracy) == 42.5)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'LOOCV classification result dataset 2D  inc. mask... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check LOOCV classification result dataset 2D inc. mask ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
  end

  
  
  
  
end %endif scope
  
  



  
  
%% Scope: RFE_SVM

if(scopesToCheck.check_RFE_SVM)
  
  disp('********** SCOPE: RFE_SVM ************');
  testScope = 'RFE_SVM';
  
  splitterOSO = getDataSplitter(datasetTest, 'oneSampleOut');

  tic
  [ds, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = ...
      doRecursiveFeatureElemination_SVM(datasetTest, 3, 45, splitterOSO, 'classification', 'linear', 0.5);
  toc
  printResultStruct(resultStruct);
  
   if(resultStruct.accuracy == 82.5)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'RFE classification result dataset 4D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check RFE classification result dataset 4D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
   end
  
   tic
   [ds, resultStruct, avg_rfe_weightMap, avg_rfe_featureSelectionMap, rfe_weightMaps, rfe_featureSelectionMaps] = ...
   doRecursiveFeatureElemination_SVM(datasetTest2D, 2, 25, splitterOSO, 'classification', 'linear', 0.5);
   toc  
   printResultStruct(resultStruct);
  
   if(single(resultStruct.accuracy) == 55)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'RFE classification result dataset 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: check RFE classification result dataset 2D ... FAILED';
    msgIndex =  msgIndex + 1;
    errorCount = errorCount +1;
   end
end %endif scope  
   
   


%% Scope: RFE_BOOTSTRAP_SVM
if(scopesToCheck.check_RFE_BOOTSTRAP_SVM)
  
  disp('********** SCOPE: RFE_BOOTSTRAP_SVM ************');
  testScope = 'RFE_BOOTSTRAP_SVM';
  
  %[dataset, resultStruct] = doRecursiveFeatureElemination_bootStrap_SVM(datasetTest2D);
  
  datasetTest = setDataset_mask_ByMatrix(datasetTest, roi3D);
  printDatasetInfo(datasetTest);
  
  %[dataset, resultStruct] = doRecursiveFeatureElemination_bootStrap_SVM(datasetTest, 20);
  
  splitterOSO = getDataSplitter(datasetTest, 'oneSampleOut');
  [dataset, resultStruct] = doRecursiveFeatureElemination_bootstrapping_SVM(datasetTest, 100, 2, 30, splitterOSO, 'linear', 0.5);
  
  printResultStruct(resultStruct);
  
end





%% Scope:check_Searchlight_SVM
if(scopesToCheck.check_Searchlight_SVM)
  
  disp('********** SCOPE: check_Searchlight_SVM ************');
  testScope = 'check_Searchlight_SVM';
  
  %[dataset, resultStruct] = doRecursiveFeatureElemination_bootStrap_SVM(datasetTest2D);
  
  datasetTest = setDataset_mask_ByMatrix(datasetTest, grayMatterMask);
  printDatasetInfo(datasetTest);
  
  
  splitterOSO = getDataSplitter(datasetTest, 'oneSampleOut');
  splitterOCO = getDataSplitter(datasetTest, 'oneChunkOut');
  [dsX, resAccuracyMap, resultStruct] = doSearchlightLOOCV_SVM(datasetTest, 3, splitterOCO, 'classification', 'linear', 0.5);
  
  
  showDataAsImage(resAccuracyMap, 'searchlightAccuracy');
  
end





%% Scope: Configurator
  
if(scopesToCheck.check_Configurator)
  
  % ----------- CONFIGURATOR -----------

  % -- BASICS -----------
  Comments = 'ROI TEST';
  baseDirectory = 'C:\maurice\Programming\Matlab\easyupMVPA\util\unitTest\testData\configuratorTest\functional\';
  %subjectlist = mpc_getsubjects(baseDirectory,4); 
  subjectlist ={'BB4T'}
  
  %{'MT3T','PF5T','QC1T','RSDT','SL7T','SMWT','WC8T'};
  
  conditions = {'crave_tasty','ncrave_tasty'};%,'crave_ntasty','ncrave_ntasty' 
  comparisons = [12]; %12 = ct vs nct ; 12 34 = ct vs nct und cnt vs ncnt
  cpu_cores = 2;

  % --  SCANS -----------
  dataFileFormat = 'data.nii'; 
  spm_betaformat = 0; 
  runs(1).directory = 's1'; 
%   runs(2).directory = 'functional_c/s2';
%   runs(3).directory = 'functional_c/s3';
%   runs(4).directory = 'functional_c/s4';
%   runs(5).directory = 'functional_c/s5';
%   runs(6).directory = 'functional_c/s6';
  TR = 2; %...
  
  %relative path means subject specific masks
  mask = 'brainMask_53_63_46.nii'; 

  
  
  % --  PREPROCESSING -----------
  detrending = 0; 
  highpassFilter = 0; 
  z_scoring = 0; 
  savePreProcessing = 0;

  % --  ONSETS, DURATIONS, DATA GROUPING -----------
  useSpmBetaFiles = 0; 
  
  %given in seconds
  sameTrialDurations = 3;
  
  %given in scans
  removeTrailingTransitions = 0; 
  removeTailTransitions = 0;
  
  %at the moment used for trial duration and onsets
  unitofdesign = 'scans'; 
  
  %0== first scan is 1, 1== first scan is 0
  spm_onsets_specification = 1;
  
  
  averaging = 0; 
  averagingSplitMethod = 'chunks'; 



  % -- CLASSIFICATION -----------
  BaseMethod = 'SL';
  %BaseMethod = 'SL-swap';

  %BaseMethod = 'ROI';
  %BaseMethod = 'parcellation';
  %BaseMethod = 'RFE';
  SearchLightDiameter = 5; %voxel space
  
  %chunks, runs, samples, or e.g. split3 
  loocvSplitMethod = 'chunks';

  % -- GROUP STATS
  secondLevelSmoothingFWHM = 4; %nach der classification wird fuer die gruppenstat gesmootht.
  accuracymapsavename = 'accmap';

  %-- init
  easyupMVPA_init('nmbCores', cpu_cores, 'quietMode', false);

  %-- get the time (needed later for the analysis subfolder)
  timeVector = clock;

  %-- create configParameters STRUCT (collect all specified variables and
  %put them in a structure)
  configSpecifiedVariables=who;
  for i=1:numel(configSpecifiedVariables)
      if strcmp(configSpecifiedVariables{i},'configParameters') && strcmp(configSpecifiedVariables{i},'configSpecifiedVariables')
      else
          if numel(eval(configSpecifiedVariables{i})) > 1000
          else
              eval(['configParameters.',configSpecifiedVariables{i},'=',configSpecifiedVariables{i},';']);
          end
      end
  end

  %-- put the configParameters structure into the parser
  myDataset = parseconfigParameters(configParameters);

end %end check configurator




%% Scope: Regression
if(scopesToCheck.check_regression)
  
  testScope = 'REGRESSION';
  
  %load the bodyfat example 
  dataFile_SVMStyle = fullfile(dir_data, 'bodyfat_scaled.txt'); 
  [bfValues,bfInData] = libsvmread(dataFile_SVMStyle);
  
  
  %create 2D dataset
  regDS = getEmpty2DDataset();
  regDS = setDataset_data_ByMatrix(regDS, bfInData');
  
  
  regDS = setDataset_classIDs_ByVector(regDS, bfValues);
  regDS = setDataset_chunks_ByVector(regDS, [1:length(bfValues)]);
  
  printDatasetInfo(regDS);
  
  dsTrain = selectSamples(regDS, 'index <= 200');
  [svmModel, weights] = train_SVM(dsTrain, 'regression_epsilon', 'linear', 1, {'epsilon', 0.01});

  dsTest = selectSamples(regDS, 'index > 200');
  [resultStruct, probEstimates] = predict_SVM(dsTest, svmModel);  
  
  printResultStruct(resultStruct);
  
  %compute correlation of predicted and real results
  cMatrix = corrcoef(resultStruct.predictedClassIDs, dsTest.classIDs);
  if(cMatrix(1,2) > 0.9775)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'Simple Regression 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: Simple Regression 2D ... FAILED';
    msgIndex =  msgIndex + 1;    
    errorCount = errorCount +1;
  end
  
  if(graphOut)
    
    figureID = 'easyup_unitTest_regressionFig';
    if(isempty(findobj('Tag', figureID)))    
        resFig = figure('Tag', figureID, 'Name', 'Results for easyup unit test regression');
    else
        figure(findobj('Tag', figureID));
    end
  
    hold on;
    plot(resultStruct.predictedClassIDs, '-r');
    plot(bfValues(201:end), '-b');
    legend('predicted', 'data', 2);
    hold off;
  end
  
  
  %LOOCV regression
  splitterOSO = getDataSplitter(regDS, 'oneSampleOut');
  [regDSX, resultStruct, avgWeights] = doLeaveOneOutCrossValidation_SVM(regDS, splitterOSO, 'regression_epsilon', 'linear', 1, {'epsilon', 0.01});
  
  printResultStruct(resultStruct);
  
  %compute correlation of predicted and real results
  cMatrix = corrcoef(resultStruct.predictedClassIDs, regDS.classIDs);
  
  disp('LOOCV correlation: ');
  disp(cMatrix(1,2));
  
  if(cMatrix(1,2) > 0.947)
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'LOOCV Regression 2D ... PASSED';
    msgIndex =  msgIndex + 1;
  else
    messageStack(msgIndex).testScope   = testScope;
    messageStack(msgIndex).testMessage = 'ERROR: LOOCV Regression 2D ... FAILED';
    msgIndex =  msgIndex + 1;    
    errorCount = errorCount +1;
  end
  
  if(graphOut)
    
    figureID = 'easyup_unitTest_regressionLOOCVFig';
    if(isempty(findobj('Tag', figureID)))    
        resFig = figure('Tag', figureID, 'Name', 'Results for easyup unit test regression');
    else
        figure(findobj('Tag', figureID));
    end
  
    hold on;
    plot(resultStruct.predictedClassIDs, '-r');
    plot(bfValues, '-b');
    legend('predicted', 'data', 2);
    hold off;
  end


end % end scope regression








  disp(' ');
  disp('************ TEST MESSAGES ***************')
  disp(' ');
  disp('Checked Scopes:');
  fn = fieldnames(scopesToCheck);
  for i=1:size(fn,1)
    fnString = char(fn(i));
    if(eval(['scopesToCheck.',fnString]))
      disp([fnString,'...YES']);
    else
      disp([fnString,'...NO']);
    end
  end
  
  disp(' ');
  for i=1:size(messageStack,2)
    disp(['Scope: ',messageStack(i).testScope,': ',messageStack(i).testMessage]);
  end  
  disp(' ');
  
  if(errorCount==0)
    disp('Finished test with no errors!');
  else
    disp(['Finished test with ',num2str(errorCount), ' error(s)!']);
  end
  disp(' ');
  toc
  
  
end

function data = member_normalize(d)
% scale before svm
% the data is normalized so that max is 1, and min is 0
data = (d -repmat(min(d,[],1),size(d,1),1))*spdiags(1./(max(d,[],1)-min(d,[],1))',0,size(d,2),size(d,2));
end

