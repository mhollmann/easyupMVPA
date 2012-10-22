%bla bla

%function [] = doParcellationDecoding(myDataset,class1,class2);


class1 = 1;
class2 = 2;


%get mask again (this time with all integer values
tmp = load_nii(fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.mask));
myDataset.mask = int8(tmp.img);
allLabels = unique(nonzeros(myDataset.mask));

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









for label=1:1
    ind = 






