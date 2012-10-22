% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%  Low-level function for getting the SVM-command string and infos about the paramStruct.
%
%
% Parameters:
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
%
function [inputIsValid, svmParamInfoStruct, commandString] = getSVMParamInfo(kernelMode, costParam, paramStruct)
   
  commandString = '';
  svmParamInfoStruct = struct('kernelMode', '', 'degree', [], 'gamma', [], 'coef0', [], 'nu', [], 'probEstimates', []);
  inputIsValid = true;
  
  
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