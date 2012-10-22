%adds leading zeros b4 number (nr = number, decimals = how many digits
%there should be)
function [numbero] = num2str_lazy(nr,decimals);


numbero=num2str(nr,['%0',num2str(decimals),'d']);
