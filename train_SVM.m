% Simple training of a dataset based on the given dataset using a Support Vector Machine (SVM).
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%  This high-level function trains a Support Vector Machine (SVM). The return value svmModel may be used
%  to predict data based on the model learned here.
%
%
% Parameters:
%   dataset        - The dataset to work on 2D or 4D (all samples are included in LOOCV)
%   svmType        - Types:
%                     ['classification', 'regression_epsilon', 'regression_nu']
%
%   kernelMode     - Kernels: 
%                     ['linear', 'polynomial', 'radial', 'sigmoid']
% 	                    linear               : u'*v
% 	                    polynomial           : (gamma*u'*v + coef0)^degree
% 	                    radial basis function: exp(-gamma*|u-v|^2)
% 	                    sigmoid              : tanh(gamma*u'*v + coef0)
%
%   costParam      - The slack variable C in SVM (range 0 to 1  0 = low cost, 1 = highest costs). 
%                    It defines the costs for misclassification (How strongly are outliers punished?).
%   paramStruct    - example: {'degree', 3, 'gamma', 0.02, 'probEstimates', 1}
%                    possible fields:
%                    'degree'        : default=3    Describes the exponent in polynomial kernel function
%                    'gamma'         : default=1/k  The gamma-factor in kernel function
%                    'coef0'         : default=0    The coefficient summand in kernel function
%                    'nu'            : default
%                    'probEstimates' : default=0    1 if probabilistic estimates should be computed, 0 if not
% 
%
% Returns:
%   svmModel       - a struct containing all model data describing the libsvm model
%   weights        - 1D or 3D map of the weights extracted for the trained model
%
% Comments:
%
function [svmModel, weights] = train_SVM(dataset, svmType, kernelMode, costParam, paramStruct)

  if( ~exist('dataset','var') || ~exist('svmType','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    error('Usage of train_SVM: [svmModel, weights] = train_SVM(dataset, svmType - [classification, regression_epsilon, regression_nu], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %extractt the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(svmType, kernelMode, costParam, paramStruct);
  end
  if(~paramStructIsValid)
    error('Usage of train_SVM: [svmModel, weights3D] = train_SVM(dataset, svmType - [classification, regression_epsilon, regression_nu], kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %use quiet mode (no outputs)
  cmdString = [cmdString, ' -q '];
  
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

  weights = zeros(sizeData(1:end-1));

  %train on all samples inside the set
  svmModel  = svmtrain(double(dataset.classIDs)', trainData2D, cmdString);
  %svmModel.svmType = svmType;
  
  %extract the weights
  weights_ = svmModel.SVs' * svmModel.sv_coef;
  weights(find(dataset.featureSelectionMap >0 )) = weights_;  

end



