function [accuracy] = doROIdecoding(myDataset,class1,class2)

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
functionaldata = myDataset.data(:,:,:,ind_samples);

%-- put data into the right format, i.e. features in x and observations in
%y
allLocations = find(myDataset.mask > 0);
FeatureCount = numel(allLocations);
flatcube_data = zeros(size(functionaldata,4),FeatureCount);
for i=1:FeatureCount
    [h,k,l] = ind2sub(size(myDataset.mask),allLocations(i));
    flatcube_data(:,i) = squeeze(functionaldata(h,k,l,:));
end

crossvalidations = size(myDataset.dataSplitter.splitMatrix,1);
accuracyvector = zeros(crossvalidations,1);

if(~easyupMVPA_getGlobals('quietMode'))
    disp(['INFO: starting decoding now using ',num2str(crossvalidations),' CrossValidations and ',num2str(FeatureCount),' Voxels']);
end

%-- LOOCV BUSINESS
for cv=1:crossvalidations
    splitMatrixReducedSlice = splitMatrixReduced(cv,:);
    testcount = numel(find(splitMatrixReducedSlice == 1));
    traincount = numel(find(splitMatrixReducedSlice == 2));
    
    %define train and test
    test_data = flatcube_data(find(splitMatrixReducedSlice == 1),:);
    train_data = flatcube_data(find(splitMatrixReducedSlice == 2),:);
   

    train_label = classIDsReduced(find(splitMatrixReducedSlice == 2))';
    test_label = classIDsReduced(find(splitMatrixReducedSlice == 1))';
    
    %LIBSVM
    model = svmtrain(train_label, double(train_data),'-s 0 -t 0 -b 0 -c 4 -e 0.1');
    [predicted_label, accuracy, decision_values] = svmpredict(test_label, double(test_data), model);
    
    
    accuracyvector(cv) = accuracy(1)-50;
end

accuracy=mean(accuracyvector);
accuracystd = std(accuracyvector);




    