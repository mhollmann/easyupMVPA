% Initializes the toolbox (setting number of cores for parallel processing and toolbox messages).
%
% Author: Maurice Hollmann
% Date  : 08/10
%
% Description:
%  easyupMVPA_init('nmbCores', 4, 'quietMode', true)
%
%  Sets the initial state. 
%  Via input parameter "nmbUsedCores" one may define how many cores to use for parallelization (if parallel distr. toolbox is available).
%  If no value is given the maximal number of cores is used. If you dont want to use parallel computing call with: easyupMVPA_init(0)
%  
% Parameters:
%   nmbCores    - [INT] the number of cores that should be used for parallel processing (0 means don't use parallel processing) 
%                 If nothing is given, the toolbox will use the maximum number of cores.
%   quietMode   - true/false true if messages of the toolbox should be supressed
%
% Returns:
%                            
%
% Comments:
function easyupMVPA_init(varargin)
  
  
  %set path variable
  dirEU = which('easyupMVPA_init');
  [pathstr, name, ext] = fileparts(dirEU);
  addpath(genpath(pathstr));
  
  easyupMVPA_version();
  
  easyupMVPA_setGlobals('quietMode', false);
  
  for i=1:size(varargin,2)
     
     %Read input var
     try
       varName = char(varargin(1,i));
     catch
       varName = 'X';
     end
     
     switch varName
       case 'nmbCores'  
          try 
            nmbUsedCores = cell2mat(varargin(1,i+1));
            if(~isnumeric(nmbUsedCores))
              clear('nmbUsedCores');
            end
          catch
          end
       case 'quietMode'
         try 
           quietMode = cell2mat(varargin(1,i+1));
           if(quietMode)
             easyupMVPA_setGlobals('quietMode', true);
           else
             easyupMVPA_setGlobals('quietMode', false);
           end
         catch
           easyupMVPA_setGlobals('quietMode', false);
         end
       case 'X'
         %do nothing
       otherwise
         error(['Invalid option : ''',varName,'''. Usage of easyupMVPA_init example: easyupMVPA_init() or easyupMVPA_init(''nmbCores'', 4, ''quietMode'', true).']);
     end
     
  end
  
  %check if matlabpool command is available otherwise
  %the distributed computing toolbox is not available
  try
    sPool = matlabpool('size');
  catch
    disp('INFO: Seems as the Distributed Computing Toolbox is not available. Running easyUpMVPA without parallelization.');
    return;
  end
  
  if(exist('nmbUsedCores','var') && nmbUsedCores>0)
    outP = findResource();
    if(nmbUsedCores > outP.ClusterSize || nmbUsedCores < 0) 
      nmbUsedCores = outP.ClusterSize;
      warning(['Restricted number of used cores to maximum available cores: ', num2str(nmbUsedCores),'.']);
    end
  end
  if(~exist('nmbUsedCores','var'))
    outP = findResource();
    nmbUsedCores = outP.ClusterSize;
  end
  
  sPool = matlabpool('size');
  
  %try to open a parallel session if intended
  if( sPool==0 && nmbUsedCores>0)
    matlabpool('open', nmbUsedCores);    
  elseif(sPool ~= nmbUsedCores && nmbUsedCores>0)
    matlabpool close;
    matlabpool('open', nmbUsedCores); 
  elseif(sPool > 0 && nmbUsedCores==0)
    matlabpool close;
  end
    
  if(matlabpool('size') > 0)
    disp(['INFO: running easyUpMVPA parallelized using ', num2str(matlabpool('size')),' cores. You can set the number of maximal cores to N using easyupMVPA_init(''nmbCores'', N)!']);
  else
    disp('INFO: running easyUpMVPA without parallelization.');
  end
  
  
end