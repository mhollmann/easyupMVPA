% Returns values for global properties in the toolbox.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%  Returns global properties. To set properties use : easyupMVPA_setGlobals(propertyName, propertyValue)
%
%  Properties are:
%  
%  easyupMVPA_quietMode  - TRUE/FALSE if true all outputs of routines are suppressed
%
%  
%
% Parameters:
%   propertyName       -> Name of property to be returned
%
% Returns:
%   ok                 -> returns true when setting worked, otherwise false
%                            
%
% Comments:
function [propertyValue] = easyupMVPA_getGlobals(propertyName)
  
  global propertyStruct;  
  
  if( isfield(propertyStruct, propertyName) )
    propertyValue = propertyStruct.(propertyName);
  else
    disp(['INFO: easyupMVAP_getGlobals: Variable "',propertyName,'" does not exist, returning FALSE !']);
    propertyValue = false;
  end
  
end