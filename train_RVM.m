% (UNDER CONSTRUCTION) Simple training of a dataset based on the given dataset using a Relevance Vector Machine (RVM).
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%  This high-level function trains a Rupport Vector Machine (RVM). The return value rvmModel may be used
%  to predict data based on the model learned here.
%
%
%
function [rvmModel, weights] = train_RVM(dataset, kernelType, kernelWidth, paramStruct)

  if(strcmp(kernelType, 'gauss'))
    kernelType = '+gauss';
  end

  if( ~exist('dataset','var'))% || ~exist('kernelMode','var') ) 
    error('Usage of train_RVM: UPDATE!!!!!! [svmModel, weights3D] = train_SVM(dataset, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end

  
  %get the data as 2D array by using the selection-map or the mask 
  sizeData   = size(dataset.data);
  if(isfield(dataset,'featureSelectionMap') && ~isempty(dataset.featureSelectionMap))
    [dataset, trainData2D] = selectFeaturesBySelectionMap(dataset);
  elseif(isfield(dataset,'mask') && ~isempty(dataset.mask))
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

  
  %keyboard;
  
  
  % *************Modelling using SPARSE BAYES 1.00**********************************
%   maxIts	= 500;
%   monIts  = round(maxIts/10);
%   initAlpha	= (1/N)^2;
%   initBeta	= 0;	% setting to zero for classification
%   %train the RVM 
%   [weights, used, marginal, alpha, beta, gamma] = sbl_rvm(X,t,initAlpha,initBeta,kType,kWidth, maxIts, monIts);
% 
% 
%    %at first apply kernel to the training set
%    %nmbSamples x fspacedim
%    
%    %extract basis set for given kernel
%    basisSet	= sbl_kernelFunction(trainData2D,trainData2D,kType,kWidth); 
%    
%    [weights, used, marginal, alpha, beta, gamma]	= sbl_estimate(basisSet,double(dataset.classIDs'),alpha,beta,maxIts,monIts);
%    keyboard; 
%    
%    rvmModel = {};
%    rvmModel.kernelType      = kType;
%    rvmModel.kernelWidth     = kWidth;
%    rvmModel.relevantIndices = PARAMETER.Relevant;
%    %rvmModel.weights         = weights;
%    rvmModel.weights         = PARAMETER.Value;
%    rvmModel.X               = trainData2D;
%    rvmModel.t               = double(dataset.classIDs');
   
   
  
  
  


  % *************Modelling using SPARSE BAYES 2.00**********************************
   
  % - we set the diagnostics level to 2 (reasonable)
  % - we will monitor the progress every 10 iterations
  %   pSparse		= 0.90;
  % 
   rvmOptions = SB2_UserOptions('iterations',800,...
							                  'diagnosticLevel', 4,...
 							                  'monitor', 10);
  %extract basis set for given kernel
  basisSet	= sbl_kernelFunction(trainData2D,trainData2D,kernelType,kernelWidth); 

  
  %Bernoulli fo 2 class classification
  [PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Bernoulli', basisSet, double(dataset.classIDs'), rvmOptions);
  
  
  %[PARAMETER, HYPERPARAMETER, DIAGNOSTIC] = SparseBayes('Gaussian', basisSet, double(dataset.classIDs'), rvmOptions);
  
%   %
%   % Manipulate the returned weights for convenience later
%   %
   %weights = zeros(length(dataset.classIDs),1);
   %weights(PARAMETER.Relevant)	= PARAMETER.Value;
   weights = [];
   
   
   relevantIndices = PARAMETER.Relevant;
   
   if kernelType(1)=='+'	
    % Take account of bias if originally used ...
    relevantIndices	= relevantIndices - 1;
    if relevantIndices(1)~=0
      % ... and if pruned ...
      kernelType(1)	= [];
    else
      relevantIndices(1)	= [];
    end
  end
   
   
   rvmModel = {};
   rvmModel.kernelType      = kernelType;
   rvmModel.kernelWidth     = kernelWidth;
   rvmModel.relevantIndices = relevantIndices;
   rvmModel.weights         = PARAMETER.Value;
   %rvmModel.X               = trainData2D;
   rvmModel.sparseTrainData = trainData2D(rvmModel.relevantIndices,:);
   rvmModel.trainClassIDs   = double(dataset.classIDs');
   
   
 
   
  % keyboard;
   
  
  % ***********************************************

end