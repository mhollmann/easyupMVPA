% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%  Low-level function for getting the SVM-command string and infos about the paramStruct.
%
%
% Parameters:
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
%                    'probEstimates' : default=0    1 if probabilistic estimates should be computed, 0 if not%
% Returns:
%   inputIsValid       - true if input was a valid set of SVM parameters, false otherwise
%   svmParamInfoStruct - a struct holding all infos about the parameters given
%   commandString      - the command string describing the input struct, this string can be used as input to SVM interface
%
%
% Comments: 
%   SVM info:
%
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
function [inputIsValid, svmParamInfoStruct, commandString] = getSVMParamInfo(svmType, kernelMode, costParam, paramStruct)
   
  commandString = '';
  svmParamInfoStruct = struct('svmType', '','kernelMode', '', 'degree', [], 'gamma', [], 'coef0', [], 'nu', [], 'probEstimates', []);
  inputIsValid = true;
  
  
  switch svmType
    case 'classification'
      commandString = [commandString, ' -s 0'];
      svmParamInfoStruct.svmType = 'classification';      
    case 'regression_epsilon'
      commandString = [commandString, ' -s 3'];
      svmParamInfoStruct.svmType = 'regression_epsilon';      
    case 'regression_nu'
      commandString = [commandString, ' -s 4'];
      svmParamInfoStruct.svmType = 'regression_nu';          
    otherwise
      inputIsValid = false;
      return;
  end
  
  
  switch kernelMode
    case 'linear'
      commandString = [commandString, ' -t 0'];
      svmParamInfoStruct.kernelMode = 'linear';
    case 'polynomial'
      commandString = [commandString, ' -t 1'];
      svmParamInfoStruct.kernelMode = 'polynomial';
    case 'radial'
      commandString = [commandString, ' -t 2'];
      svmParamInfoStruct.kernelMode = 'radial';
    case 'sigmoid'
      commandString = [commandString, ' -t 3'];
      svmParamInfoStruct.kernelMode = 'sigmoid';
    otherwise
      inputIsValid = false;
      return;
  end
  
  if(exist('costParam', 'var'))
    commandString = [commandString, ' -c ', num2str(costParam)];
  else
    inputIsValid = false;
    return;
  end
  
  %take a look at the paramStruct
  if(~isempty(paramStruct))

    for i=1:size(paramStruct,2)
      
      if(i==size(paramStruct,2)) 
        break; 
      end
      
      %if input is a char array the following value is used
      if(ischar(paramStruct{i}) && isnumeric(paramStruct{i+1}))
        
        optString = paramStruct{i};
        optValue  = paramStruct{i+1};
        
        if(strcmp(optString, 'degree'))
          if(optValue < 1)
            disp('INFO: Ignoring invalid SVM option value for degree!');
          else
            commandString = [commandString, ' -d ', num2str(optValue)];
            svmParamInfoStruct.degree = optValue;
          end
          
        elseif(strcmp(optString, 'gamma'))
          if(optValue > 1)
            disp('INFO: Ignoring invalid SVM option value for gamma!');
          else
            commandString = [commandString, ' -g ', num2str(optValue)];
            svmParamInfoStruct.gamma = optValue;
          end
          
        elseif(strcmp(optString, 'coef0'))
          if(optValue > 1)
            disp('INFO: Ignoring invalid SVM option value for coef0!');
          else
            commandString = [commandString, ' -r ', num2str(optValue)];
            svmParamInfoStruct.coef0 = optValue;
          end

        elseif(strcmp(optString, 'nu'))
          if(optValue > 1)
            disp('INFO: Ignoring invalid SVM option value for nu !');
          else
            commandString = [commandString, ' -n ', num2str(optValue)];
            svmParamInfoStruct.nu = optValue;
          end

        elseif(strcmp(optString, 'epsilon'))
          if(optValue < 0)
            disp('INFO: Ignoring invalid SVM option value for epsilon !');
          else
            commandString = [commandString, ' -p ', num2str(optValue)];
            svmParamInfoStruct.nu = optValue;
          end

        elseif(strcmp(optString, 'probEstimates'))
          if(optValue == 0 || optValue == 1)
            commandString = [commandString, ' -b ', num2str(optValue)];
            svmParamInfoStruct.probEstimates = optValue;
          else
            disp('Ignoring invalid SVM option value for probEstimates!');
          end

        end
      end
      
    end
  end
  
end