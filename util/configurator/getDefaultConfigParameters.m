% Checks which Configparameters have been specified and fills the rest up
% with default parameters
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%


function [configParameters] = getDefaultConfigParameters(configParameters);
%completes the empty fields of userdata


timeVector=configParameters.timeVector;
configParameters.timeString = [num2str(timeVector(1)),'-',num2str_lazy(timeVector(2),2),'-',num2str_lazy(timeVector(3),2),'-',num2str_lazy(timeVector(4),2),num2str_lazy(timeVector(5),2)];


%load default values if unspecified above
if isfield(configParameters,'cpu_cores') == 0; configParameters.cpu_cores = 1;end
if isfield(configParameters,'dir_analysis_mvpa') == 0; configParameters.dir_analysis_mvpa   = 'analysis_mvpa'; end
if isfield(configParameters,'dataset_name') == 0; configParameters.dataset_name   = 'dataset'; end
if isfield(configParameters,'averaging') == 0; configParameters.averaging   = 0; end
if isfield(configParameters,'secondLevelSmoothingFWHM') == 0; configParameters.secondLevelSmoothingFWHM   = 0; end
if isfield(configParameters,'useSpmBetaFiles') == 0; configParameters.useSpmBetaFiles   = 0; end
if isfield(configParameters,'Comments') == 0; configParameters.Comments   = ''; end
if isfield(configParameters,'SearchLightDiameter') == 0; configParameters.SearchLightDiameter   = 5; end
if isfield(configParameters,'obtainmask') == 0; configParameters.obtainmask   = 0; end
if isfield(configParameters,'accuracymapsavename') == 0; configParameters.accuracymapsavename   = 'accmap'; end
if isfield(configParameters,'customAnalysisSubfolder') == 0; configParameters.analysisSubfolder   = configParameters.timeString; else  configParameters.analysisSubfolder = configParameters.customAnalysisSubfolder; end
if isfield(configParameters,'Searchlightcachesize') == 0; configParameters.Searchlightcachesize   = 10000; end




