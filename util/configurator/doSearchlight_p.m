% Performs a searchlight decoding for class1 versus class2. 
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:
% For every ocation (every value of the mask > 0) a spherical neighbourhood of
% voxels is cut out and a leave-one-out cross validation is performed
% (using the samples specified in the LOOCV splitmethod). The result is a
% mean accuracy map, which has values between -50 and +50 (every time wrong
% to every time correct, 0 is the chance level).
%
%
% Parameters:
%   dataset     - the dataset, class1 and class2
%
% Returns:
%   accuracymap    -a map of the same dimension of the original data. values between -50 and +50 (every time wrong
% to every time correct, 0 is the chance level).
%
% Comments:

function [accuracymap] = doSearchlight_p(myDataset,class1,class2);
% default settings
%tic

spotlightform = getSpotlightform(myDataset.configParameters.SearchLightDiameter);


%===============================

%-- get all regular (non-excluded) volumes which are in class1 and class2
ind_validsamples = find(myDataset.chunks);

ind_class1 = find(myDataset.classIDs == class1);
ind_class1 = ind_class1(ismember(ind_class1,ind_validsamples));


ind_class2 = find(myDataset.classIDs == class2);
ind_class2 = ind_class2(ismember(ind_class2,ind_validsamples));

%-- check whether sample count is identical
if numel(ind_class1) == numel(ind_class2)
    classIDsReduced = horzcat(ones(1,numel(ind_class1)),2*ones(1,numel(ind_class1)));
    ind_samples = horzcat(ind_class1,ind_class2);
else
    error('sample count is not the same for all conditions. check your data!')
end


%-- check whether split matrices are identical
if myDataset.dataSplitter.splitMatrix(:,ind_class1) == myDataset.dataSplitter.splitMatrix(:,ind_class2)
    splitMatrixReduced = horzcat(myDataset.dataSplitter.splitMatrix(:,ind_class2),myDataset.dataSplitter.splitMatrix(:,ind_class2));
else
    error('sample count is not the same for all conditions. check your data!')
end

%-- reduce the data size by copying only what is needed into a new variable
clear functionaldata
functionaldata = myDataset.data(:,:,:,ind_samples);


sz = size(spotlightform,1);
spotlightform_4d(:,:,:,1) = spotlightform;
spotlightform_4d = shiftdim(spotlightform_4d,3);
spotlightform_flat = reshape(spotlightform_4d,1,sz^3);
flatcube_killcols = find(spotlightform_flat == 0);
spotlightrad = floor(sz/2);
crossvalidations = size(myDataset.dataSplitter.splitMatrix,1);

if(~easyupMVPA_getGlobals('quietMode'))
    disp(['INFO: Starting searchlight procedure using ',num2str(crossvalidations),' CrossValidations']);
end

spotlightcount = sz^3;% - numel(find(flatcube_killcols));


%--make functionaldata double and replace NaN with zeros
functionaldata=manipulateBox(functionaldata,sz,'increase');
functionaldata = double(functionaldata);
functionaldata(isnan(functionaldata))=0;
mask = manipulateBox(myDataset.mask,sz,'increase');
hsize = size(functionaldata,1);
ksize = size(functionaldata,2);
lsize = size(functionaldata,3);
clear accuracymap
accuracymap = zeros(hsize,ksize,lsize);


%--now checkout all locations
allSearchlightLocations = find(mask > 0);

%--and split these locations into pieces of length LocSplitSize. LocSplit
%defines the splitpoints of the indexes
LocSplitSize = myDataset.configParameters.Searchlightcachesize;
LocSplitCount = ceil(numel(allSearchlightLocations)/LocSplitSize);

LocSplit = zeros(LocSplitCount,2);
LocSplit(:,1) = 1:LocSplitSize:numel(allSearchlightLocations);
LocSplit(1:end-1,2) = (LocSplitSize+1:LocSplitSize:numel(allSearchlightLocations))-1;
LocSplit(end,2) = numel(allSearchlightLocations);

%pre-fill the SlicedData-Cell
EmptySlice = zeros(sz,sz,sz);
EmptySlicedData = cell(LocSplitSize);
for i = 1:LocSplitSize
   EmptySlicedData{i} = EmptySlice;
end


for s = 1:LocSplitCount
    
    
    %prepare the data, put it into a cell array
    cellindex = 1;
    SplitElements = (LocSplit(s,2)-LocSplit(s,1)+1);
    SlicedData = EmptySlicedData;

    for i = LocSplit(s,1):LocSplit(s,2)
        [h,k,l] = ind2sub(size(mask),allSearchlightLocations(i));
        
        SlicedData{cellindex} = functionaldata(h-spotlightrad:h+spotlightrad,k-spotlightrad:k+spotlightrad,l-spotlightrad:l+spotlightrad,:);
        cellindex = cellindex + 1;
    end
    
    accuracies = zeros(1,SplitElements);
    
    clear TempDataSlice flatcube_train flatcube_test test_data train_data
    
    %PARFOR
    parfor i = 1:SplitElements
        
        TempDataSlice = SlicedData{i};
        accuracyvector = zeros(crossvalidations,1);
        
        for cv=1:crossvalidations
            splitMatrixReducedSlice = splitMatrixReduced(cv,:);
            testcount = numel(find(splitMatrixReducedSlice == 1));
            traincount = numel(find(splitMatrixReducedSlice == 2));
            
            %define train and test
            test_data = TempDataSlice(:,:,:,find(splitMatrixReducedSlice == 1));
            train_data = TempDataSlice(:,:,:,find(splitMatrixReducedSlice == 2));
            
            
            test_data = shiftdim(test_data,3);
            train_data = shiftdim(train_data,3);
            
            flatcube_test = reshape(test_data,testcount,spotlightcount);
            flatcube_test(:,flatcube_killcols) = [];
            
            flatcube_train = reshape(train_data,traincount,spotlightcount);
            flatcube_train(:,flatcube_killcols) = [];
        
            train_label = classIDsReduced(find(splitMatrixReducedSlice == 2))';
            test_label = classIDsReduced(find(splitMatrixReducedSlice == 1))';
            
            %LIBSVM
            
            model = svmtrain(train_label, double(flatcube_train),'-s 1 -t 0 -b 0');
            [predicted_label, accuracy, decision_values] = svmpredict(test_label, double(flatcube_test), model);
            
            
           
            accuracyvector(cv) = accuracy(1)-50;
        end
        
       accuracies(i)=mean(accuracyvector);
        
    end
    
    
    
    LocAccIndex = 1:SplitElements;
    LocAccIndex = LocAccIndex - 1 + LocSplit(s,1);
    accuracymap(allSearchlightLocations(LocAccIndex)) = accuracies;
    
end
accuracymap = manipulateBox(accuracymap,sz,'decrease');

clearvars -except accuracymap



    
