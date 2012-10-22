% Averages data according to the splitting method specified in
% myDataset.dataSplitter.AveragingDataSplitScheme. This entry
% (averagingsplitsscheme) is computed by the function getSplitScheme
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:
%  This method does an averaging multiple samples of the same class into
%  one. The definition of multiple depends on the AveragingDataSplitscheme, which itself
%  depends on the configurator variable averagingSplitMethod. If it is set to "trial", then all samples of
%  ONE trial will be averaged. In case of "split-to-X", X / totalvolumes
%  samples aver averaged to one (while totalvolumes is the total number of
%  functional volumes. If it's set to "runs", then runs / totalvolumes
%  samples are averaged.
%
%
% Parameters:
%   dataset     - the dataset for averaging (types dataset2D or dataset4D)
%
% Returns:
%   dataset    - the dataset with samples averaged over chunks
%
% Comments:
%
function [myDataset] = averageOverSplitScheme(myDataset)

% --case Samples
if strfind(myDataset.configParameters.averagingSplitMethod,'ampl')
    error('Cannot average one sample over one sample!')

% --case Rest
else
    classesSpecified = max(unique(myDataset.classIDs));
    %totalvolumes = myDataset.chunks;
    if(~easyupMVPA_getGlobals('quietMode'))
        disp(['INFO: Averaging data over ', myDataset.configParameters.averagingSplitMethod]);
    end
    
    d4DIndex = 1;
    
    TotalNewSampleCount = max(myDataset.dataSplitter.AveragingDataSplitScheme) * classesSpecified;
    datadim = size(myDataset.data,1);
    newData4D   = zeros(size(myDataset.data,1),size(myDataset.data,2),size(myDataset.data,3), TotalNewSampleCount);
    newChunks   = zeros(1,TotalNewSampleCount);
    newClassIDs = zeros(1,TotalNewSampleCount);
    

    for s = 1:max(myDataset.dataSplitter.AveragingDataSplitScheme)
        for cl=1:classesSpecified
            current_class_ind = find(myDataset.classIDs == cl);
            split_X_ind = find(myDataset.dataSplitter.AveragingDataSplitScheme == s);
            avg = mean(myDataset.data(:,:,:,current_class_ind(ismember(current_class_ind,split_X_ind))),4);
            newData4D(:,:,:,d4DIndex) = avg;
            newClassIDs(d4DIndex) = cl;
            newChunks(d4DIndex) = d4DIndex;
            d4DIndex = d4DIndex + 1;
        end
    end

    
    myDataset.data = newData4D;
    myDataset.classIDs = newClassIDs;
    myDataset.chunks = newChunks;
    
        
    
end


    