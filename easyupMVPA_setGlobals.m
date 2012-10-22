% Used to set values for global properties in the toolbox.
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%  Sets global properties. Properties are:
%  
%  easyupMVPA_quietMode  - TRUE/FALSE if true all outputs of routines are suppressed
%
%
% Parameters:
%   propertyName       -> Name of property to set
%   propertyValue      -> Value that should be given
%
% Returns:
%   ok                 -> returns true when setting worked, otherwise false
%                            
%
% Comments:
function [ok] = easyupMVPA_setGlobals(propertyName, propertyValue)
  
  global propertyStruct;
  propList = {'quietMode'};
  
  if( ~isempty(strmatch(propertyName, propList)))

    propertyStruct.(propertyName) = propertyValue;
    ok = 1;
    
  else
    error(['easyupMVAP_setGlobals: variable "',propertyName,'" could not be set, it is no global property !']);
  end
  
end