% Author: Maurice Hollmann
% Date  : 10/12
%
% Description:
%  Low-level function that returns information about a learned svm model.
%
function [svmInfoStruct] = getSVMModelInfo(svmModel)

  if( ~exist('svmModel','var')) 
    error('Usage of getSVMModelInfo: [svmModel, weights] = [svmInfoStruct] = getSVMModelInfo(svmModel)');
  end

  if(svmModel.Parameters(1) < 3)
    svmInfoStruct.svmType = 'classification';
  elseif(svmModel.Parameters(1) == 3)
    svmInfoStruct.svmType = 'regression_epsilon';
  elseif(svmModel.Parameters(1) == 4)
    svmInfoStruct.svmType = 'regression_nu';
  else
    svmInfoStruct.svmType = 'unknown';
  end
  
end