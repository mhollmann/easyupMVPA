% Set the data field of a dataset by a given myDataset.configParameters when using SPM Beta files (as observation reduction procedure).
%
% Author: Johannes Stelzer
% Date  : 05/11
%
% Description:
%
%   [dataset] = setDataset_data_ByconfigParameters(myDataset,configParameters)

%   This methods sets the data for SPM Beta values. The names of the
%   conditions specified MUST be the same as specified in SPM.

%
% Returns:
%   dataset   - the datset with included data4D and data_3DNiftiHdr - struct
%
% Comments:
%
function myDataset = setDataset_data_bySPMbetavalues(myDataset);


%-- splitmethod needs to be Samples!
splitMethod = myDataset.configParameters.loocvSplitMethod;
if strfind(splitMethod,'ampl')
else
    disp('WARNING: loocvSplitMethod usually needs to be set to <samples> for SPM BETA decoding');
end
%-- multiple runs!
if numel(myDataset.configParameters.runs) > 1
    error('specify only ONE run (directory) for SPM BETA decoding!')
end
%-- no averaging!
if myDataset.configParameters.averaging == 1
      error('set averaging to 0 for SPM BETA decoding!')
end  
%-- no trail/tail removal!
if myDataset.configParameters.removeTrailingTransitions > 0 || myDataset.configParameters.removeTailTransitions > 0
      error('removeTrailing oder removeTail has to be set to 0 for SPM BETA decoding!')
end 




tmpdir = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.runs(1).directory);

if(~easyupMVPA_getGlobals('quietMode'))
    disp(['INFO: Loading in SPM Beta files from']);
    disp(tmpdir);
end

for c = 1:numel(myDataset.configParameters.conditions)
    [tmpdata,tmpheader,allfileslist] = loadFilteredVolumes(tmpdir,'beta',myDataset.configParameters.conditions{c});
    switch c
        case 1
            runs = size(tmpdata,4);
            data = zeros(size(tmpdata,1),size(tmpdata,2),size(tmpdata,3),runs*numel(myDataset.configParameters.conditions));
            classIDs = zeros(1,runs*numel(myDataset.configParameters.conditions));
    end
    data(:,:,:,(c-1)*runs+1:c*runs) = tmpdata;
    classIDs((c-1)*runs+1:c*runs) = c;
 
end
chunks = 1:runs*numel(myDataset.configParameters.conditions);

myDataset.data = data;
myDataset.chunks = chunks;
myDataset.classIDs = classIDs;

%--get DatasplitScheme
classesSpecified=max(classIDs);
volumes=numel(chunks);
splits = runs;
DataSplitScheme = zeros(size(chunks));
for cl=1:classesSpecified
    
    %ValidSamplesOfThisClass(find(ValidSamples == cl)) = cl;
    ValidSamplesOfThisClass = classIDs == cl;
    ValidSamplesOfThisClass = ValidSamplesOfThisClass*cl;
    
    
    sampleind = find(ValidSamplesOfThisClass);
    
    %allow split-X setting (in case there is more than one beta estimate
    %per run...
    
    if strfind(splitMethod,'split')
        splits = str2num(splitMethod(regexp(splitMethod,'\d')));
        volspersplitclass = volumes/(splits*classesSpecified);
        if mod(volumes,splits)
            error('Splitting of SPM Beta values not possible due to uneven number of classes & samples');
        end
        
        for s = 1:splits
            for v = 1:volspersplitclass
                
                %(s-1)*volspersplitclass+1+(v-1)
                DataSplitScheme(sampleind((s-1)*volspersplitclass+1+(v-1))) = s;
            end
        end
    else
        for s = 1:splits
            DataSplitScheme(sampleind(s)) = s;
            %DataSplitScheme(sampleind(samplespersplit*(s-1)+1:samplespersplit*s)) = s;
        end
    end
    
end


%--put the DataSplitScheme in the usual datasplitter format, create a
%matrix
splittings = max(DataSplitScheme);
splitMatrix = zeros(splittings,volumes);
for i = 1:splittings
    splitMatrixSlice = zeros(1,volumes);
    splitMatrixSlice(find(DataSplitScheme == i)) = 1; %test set
    splitMatrixSlice(find(DataSplitScheme > 0 & DataSplitScheme ~= i)) = 2;
    splitMatrix(i,:) = splitMatrixSlice;
end

%display split populatinos
population_size = zeros(classesSpecified,splits);
for cl=1:classesSpecified
    for s=1:splits
        population_size(cl,s) = numel(find(DataSplitScheme == s & classIDs == cl));
        population_ind{cl,s} = find(DataSplitScheme == s & classIDs == cl);
    end
end

%population_ind = population_ind


dataSplitter.splitMatrix = splitMatrix;
myDataset.dataSplitter = dataSplitter;
myDataset.dataSplitter.LOOCVDataSplitScheme = DataSplitScheme;
myDataset.data_3DNiftiHdr = tmpheader.hdr;



