% Performs a searchlight decoding with cross-training: Classifier trianed
% to distinguish class1 versus class2 but tested on class3 versus class4.
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:
% For every ocation (every value of the mask > 0) a spherical neighbourhood of
% voxels is cut out and a leave-one-out cross validation is performed
% (using the samples specified in the LOOCV splitmethod). The classifier is trained on 2 conditions and the generizability is
% tested on two different conditions. The result is a
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

function [accuracymap] = doSearchlightSwap_p(myDataset,tr_c1,tr_c2,ts_c1,ts_c2);
% default settings
%tic


%values section, needs to be set later
spotlightform = mpc_spotlightformdep(myDataset.configParameters.SearchLightDiameter);


%===============================




%-- get all regular (non-excluded) volumes which are in TRAINING class1 and training class2
tr_ind_validsamples = find(myDataset.chunks);

tr_ind_class1 = find(myDataset.classIDs == tr_c1);
tr_ind_class1 = tr_ind_class1(ismember(tr_ind_class1,tr_ind_validsamples));


tr_ind_class2 = find(myDataset.classIDs == tr_c2);
tr_ind_class2 = tr_ind_class2(ismember(tr_ind_class2,tr_ind_validsamples));

%-- check whether sample count is identical
if numel(tr_ind_class1) == numel(tr_ind_class2)
    classIDsReduced = horzcat(ones(1,numel(tr_ind_class1)),2*ones(1,numel(tr_ind_class1)));
    tr_ind_samples = horzcat(tr_ind_class1,tr_ind_class2);
else
    error('sample count is not the same for all conditions. check your data!')
end


%-- check whether split matrices are identical
if myDataset.dataSplitter.splitMatrix(:,tr_ind_class1) == myDataset.dataSplitter.splitMatrix(:,tr_ind_class2)
    splitMatrixReduced = horzcat(myDataset.dataSplitter.splitMatrix(:,tr_ind_class2),myDataset.dataSplitter.splitMatrix(:,tr_ind_class2));
else
    error('sample count is not the same for all conditions. check your data!')
end

%-- reduce the data size by copying only what is needed into a new variable
tr_functionaldata = myDataset.data(:,:,:,tr_ind_samples);



%-- get all regular (non-excluded) volumes which are in TESTING class1 and training class2
ts_ind_validsamples = find(myDataset.chunks);

ts_ind_class1 = find(myDataset.classIDs == ts_c1);
ts_ind_class1 = ts_ind_class1(ismember(ts_ind_class1,ts_ind_validsamples));


ts_ind_class2 = find(myDataset.classIDs == ts_c2);
ts_ind_class2 = ts_ind_class2(ismember(ts_ind_class2,ts_ind_validsamples));

%-- check whether sample count is identical
if numel(ts_ind_class1) == numel(ts_ind_class2)
    classIDsReduced = horzcat(ones(1,numel(ts_ind_class1)),2*ones(1,numel(ts_ind_class1)));
    ts_ind_samples = horzcat(ts_ind_class1,ts_ind_class2);
else
    error('sample count is not the same for all conditions. check your data!')
end


%-- check whether split matrices are identical
if myDataset.dataSplitter.splitMatrix(:,ts_ind_class1) == myDataset.dataSplitter.splitMatrix(:,ts_ind_class2)
    splitMatrixReduced = horzcat(myDataset.dataSplitter.splitMatrix(:,ts_ind_class2),myDataset.dataSplitter.splitMatrix(:,ts_ind_class2));
else
    error('sample count is not the same for all conditions. check your data!')
end

%-- reduce the data size by copying only what is needed into a new variable
ts_functionaldata = myDataset.data(:,:,:,ts_ind_samples);



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
tr_functionaldata=manipulateBox(tr_functionaldata,sz,'increase');
tr_functionaldata = double(tr_functionaldata);
tr_functionaldata(isnan(tr_functionaldata))=0;

ts_functionaldata=manipulateBox(ts_functionaldata,sz,'increase');
ts_functionaldata = double(ts_functionaldata);
ts_functionaldata(isnan(ts_functionaldata))=0;


mask = manipulateBox(myDataset.mask,sz,'increase');
hsize = size(tr_functionaldata,1);
ksize = size(tr_functionaldata,2);
lsize = size(tr_functionaldata,3);
accuracymap = zeros(hsize,ksize,lsize);


%--now checkout all locations
allSearchlightLocations = find(mask > 0);

%--and split these locations into pieces of length LocSplitSize. LocSplit
%defines the splitpoints of the indexes
LocSplitSize = 10000;
LocSplitCount = ceil(numel(allSearchlightLocations)/LocSplitSize);

LocSplit = zeros(LocSplitCount,2);
LocSplit(:,1) = 1:LocSplitSize:numel(allSearchlightLocations);
LocSplit(1:end-1,2) = (LocSplitSize+1:LocSplitSize:numel(allSearchlightLocations))-1;
LocSplit(end,2) = numel(allSearchlightLocations);



for s = 1:LocSplitCount
    
    
    %prepare the data, put it into a cell array
    cellindex = 1;
    SplitElements = (LocSplit(s,2)-LocSplit(s,1)+1);
    tr_SlicedData = cell(SplitElements);
    ts_SlicedData = cell(SplitElements);
    for i = LocSplit(s,1):LocSplit(s,2)
        [h,k,l] = ind2sub(size(mask),allSearchlightLocations(i));
        tr_SlicedData{cellindex} = tr_functionaldata(h-spotlightrad:h+spotlightrad,k-spotlightrad:k+spotlightrad,l-spotlightrad:l+spotlightrad,:);
        ts_SlicedData{cellindex} = ts_functionaldata(h-spotlightrad:h+spotlightrad,k-spotlightrad:k+spotlightrad,l-spotlightrad:l+spotlightrad,:);
        cellindex = cellindex + 1;
    end
    
    
    accuracies = zeros(1,SplitElements);
    
    parfor i = 1:SplitElements
        tr_TempDataSlice = tr_SlicedData{i};
        ts_TempDataSlice = ts_SlicedData{i};
        accuracyvector = zeros(crossvalidations,1);
        for cv=1:crossvalidations
            splitMatrixReducedSlice = splitMatrixReduced(cv,:);
            testcount = numel(find(splitMatrixReducedSlice == 1));
            traincount = numel(find(splitMatrixReducedSlice == 2));
            
            %define train and test
            train = tr_TempDataSlice(:,:,:,find(splitMatrixReducedSlice == 2));
            test = ts_TempDataSlice(:,:,:,find(splitMatrixReducedSlice == 1));
            
            train = shiftdim(train,3);
            test = shiftdim(test,3);
            
            flatcube_train = reshape(train,traincount,spotlightcount);
            flatcube_train(:,flatcube_killcols) = [];
            
            flatcube_test = reshape(test,testcount,spotlightcount);
            flatcube_test(:,flatcube_killcols) = [];
            
            train_label = classIDsReduced(find(splitMatrixReducedSlice == 2))';
            test_label = classIDsReduced(find(splitMatrixReducedSlice == 1))';
            
            model = svmtrain(train_label, double(flatcube_train),'-s 0 -t 0 -b 0');
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




    
