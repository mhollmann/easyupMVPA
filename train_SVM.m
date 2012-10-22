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
%   dataset        - the datset with included featureSelectionMap3D (this mask just contains zeros and ones)
%   svmModel       - a struct containing all model data describing the libsvm model
%   weights        - 1D or 3D map of the weights extracted for the trained model
%
% Comments:
%
function [svmModel, weights] = train_SVM(dataset, kernelMode, costParam, paramStruct)

  if( ~exist('dataset','var') || ~exist('kernelMode','var') || ~exist('costParam','var')) 
    error('Usage of train_SVM: [svmModel, weights] = train_SVM(dataset, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
  end
  
  %extractt the SVM parameter values from paramStruct
  if( ~exist('paramStruct','var'))
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, {});
  else
    [paramStructIsValid, svmParamInfoStruct, cmdString] = getSVMParamInfo(kernelMode, costParam, paramStruct);
  end
  if(~paramStructIsValid)
    error('Usage of train_SVM: [svmModel, weights3D] = train_SVM(dataset, kernelMode - [linear, polynomial, radial, sigmoid] , costParam [0-1], paramStruct [optional - i.e. {"degree", 3}])');
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
  
  %extract the weights
  weights_ = svmModel.SVs' * svmModel.sv_coef;
  weights(find(dataset.featureSelectionMap >0 )) = weights_;  

end

%   "Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
% 	"libsvm_options:\n"
% 	"-s svm_type : set type of SVM (default 0)\n"
% 	"	0 -- C-SVC\n"
% 	"	1 -- nu-SVC\n"
% 	"	2 -- one-class SVM\n"
% 	"	3 -- epsilon-SVR\n"
% 	"	4 -- nu-SVR\n"
% 	"-t kernel_type : set type of kernel function (default 2)\n"
% 	"	0 -- linear: u'*v\n"
% 	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
% 	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
% 	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
% 	"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"
% 	"-d degree : set degree in kernel function (default 3)\n"
% 	"-g gamma : set gamma in kernel function (default 1/k)\n"
% 	"-r coef0 : set coef0 in kernel function (default 0)\n"
% 	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
% 	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
% 	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
% 	"-m cachesize : set cache memory size in MB (default 100)\n"
% 	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
% 	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
% 	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
% 	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
% 	"-v n : n-fold cross validation mode\n"
% 	"-q : quiet mode (no outputs)\n"


% 	"Usage: [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model, 'libsvm_options')\n"
% 		"Parameters:\n"
% 		"  model: SVM model structure from svmtrain.\n"
% 		"  libsvm_options:\n"
% 		"    -b probability_estimates: whether to predict probability estimates, 0 or 1 (default 0); one-class SVM not supported yet\n"
% 		"Returns:\n"
% 		"  predicted_label: SVM prediction output vector.\n"
% 		"  accuracy: a vector with accuracy, mean squared error, squared correlation coefficient.\n"
% 		"  prob_estimates: If selected, probability estimate vector.\n"

