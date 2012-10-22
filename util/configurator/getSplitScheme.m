% Returns a splitscheme for the leave-one-out-crossvalidation OR averaging and puts it
% into the dataset format. The purpose of the splitting is defined in the second input variable. The information is taken from
% configParameters.loocvSplitMethod or  configParameters.averagingSplitMethod
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:
%
%   [myDataset] = getSplitScheme(myDataset,splitPurpose)
%
% Parameters (value of configParameters.loocvSplitMethod or configParameters.averagingSplitMethod)
%
% 'Samples' : ONE sample will serve test set, the rest is taken for the training.
% Assuming N Samples, this procedure will be repeated for N times (N-fold
% crossvalidation). In case of averaging nothing will happen here.
%
% 'Chunks' : ONE Chunk (consisting of multiple samples) will serve as test set, the rest is taken as training.
% Assuming M Chunks, this procedure will be repeated for M times (M-fold
% crossvalidation). In case of averaging, chunks will get averaged.
%
% 'Split-Into-X' : All Chunks are split into X groups of Chunks (eventually
% Chunks have to be removed to assure a equal number of Chunks. In this case the killed Chunks will be drawn randomly and you'll be noticed). 
% Assuming X groups, this procedure will be repeated for X times (X-fold
% crossvalidation). In case of averaging, everything within the group will get averaged.
%
% 'Runs' : ONE run (consisting of multiple Chunks) will serve as test set, the rest is taken as training.
% Assuming R Runs, this procedure will be repeated for R times (R-fold
% crossvalidation). In case of averaging, everything within the run will get averaged.

 
function [myDataset] = getSplitScheme(myDataset,splitPurpose)

switch splitPurpose
    case 'LOOCV'
        splitMethod = myDataset.configParameters.loocvSplitMethod;
    case 'AVERAGING'
        splitMethod = myDataset.configParameters.averagingSplitMethod;
end


volumes = numel(myDataset.chunks);
DataSplitScheme = zeros(1,volumes);
chunks_kicked = zeros(1,volumes);


classesSpecified = max(unique(myDataset.classIDs));

ValidSamples = zeros(size(myDataset.chunks));
ValidSamples(find(myDataset.chunks)) = myDataset.classIDs(find(myDataset.chunks));

do_X_split = 0;

% --case Samples
if strfind(splitMethod,'ampl')
    dataSplitter.splitterType = 'Samples';
    
    %check if the samples are equally distributed over the classes
    
    for cl=1:classesSpecified
        samplecount(cl) = numel(find(ValidSamples==cl));
    end
    if numel(unique(samplecount)) == 1 %if they're uniformly distributed
        for cl=1:classesSpecified
            sampleind = find(ValidSamples == cl);
            for s=1:numel(sampleind)
                DataSplitScheme(sampleind(s))=s;
            end
        end
    else %if not... kill some here!
        for cl=1:classesSpecified
            sampleind = find(ValidSamples == cl);
            sampleind = sampleind(randperm(samplecount(cl)));
            kicked_ind = sampleind(min(samplecount)+1:end);
            chunks_kicked(kicked_ind) = 1;
            myDataset.chunks(kicked_ind) = 0;
            if numel(kicked_ind) > 0 
                warning(['kicked samples ',num2str(kicked_ind)])
            end
            sampleind = sort(sampleind(1:min(samplecount)));
            for s=1:numel(sampleind)
                DataSplitScheme(sampleind(s))=s;
            end
        end 
        
    end
    

% --case Chunks
elseif strfind(splitMethod,'hunk')
    dataSplitter.splitterType = 'Chunks';
    
    %find out what length the Chunks have (i.e. is the number of samples per Chunk
    %constant?)

    SamplesPerChunk = zeros(size(myDataset.chunks));
    
    
    for cl=1:classesSpecified
        ValidSamplesOfThisClass = zeros(size(ValidSamples));
        ValidSamplesOfThisClass(find(ValidSamples == cl)) = cl;
        index = 1;
        crossvaliteration = 1;
        while index <= numel(SamplesPerChunk)
            nextpos = 0;
            if ValidSamplesOfThisClass(index)
                stopinnerloop = 0;
                nextpos = 1;
                samplesinthisChunk = 1;
                
                while stopinnerloop == 0
                    if index + nextpos <= numel(SamplesPerChunk)
                        if ValidSamplesOfThisClass(index) == ValidSamplesOfThisClass(index + nextpos)
                            nextpos = nextpos + 1;
                            samplesinthisChunk = samplesinthisChunk + 1;
                        else
                            stopinnerloop = 1;
                            nextpos = nextpos-1;
                            SamplesPerChunk(index:index+nextpos) = samplesinthisChunk;
                            DataSplitScheme(index:index+nextpos) = crossvaliteration;
                            crossvaliteration = crossvaliteration + 1;
                        end
                    else
                        stopinnerloop = 1;
                    end
                end
                
            end
            index = index + nextpos + 1;
        end
    end
%     SamplesPerChunk
%     ValidSamples
%     DataSplitScheme
%     


%--case split-to-x
elseif strfind(splitMethod,'lit')
    splits = str2num(splitMethod(regexp(splitMethod,'\d')));
    dataSplitter.splitterType = ['Do ',num2str(splits),'-split'];
    %splits
    do_X_split = 1;
    %perform operation later (end)
    
    
    %--case runs
elseif strfind(splitMethod,'un')
    dataSplitter.splitterType = 'Runs';
    splits = numel(myDataset.dataFilelist);
    do_X_split = 1;
    %perform operation later (end)
    
    
    
    
    
else
    error('Please specify the DataSplitMethod in a clear way.')
end


if do_X_split
    for cl=1:classesSpecified
        samplecount(cl) = numel(find(ValidSamples==cl));
    end
    if numel(unique(samplecount)) > 1
        %error()
        disp(['samples not distributed equally across conditions. killig some.'])
    end
    if mod(samplecount(1),splits) > 0
        disp(['samples not distributed equally across conditions. killig some.'])
        %error('number of samples divided by the split must be an integer! (7/3 IS NOT an integer, 8/4 IS an integer...)')
    end
    
    
    
%     for cl=1:classesSpecified
%         samplespersplit = max(samplecount)/splits;
%         ValidSamplesOfThisClass = zeros(size(ValidSamples));
%        
%         %ValidSamplesOfThisClass(find(ValidSamples == cl)) = cl;
%         ValidSamplesOfThisClass = ValidSamples == cl;
%         ValidSamplesOfThisClass = ValidSamplesOfThisClass*cl;
%         
%        
%         sampleind = find(ValidSamplesOfThisClass);
%         for s = 1:splits
%             DataSplitScheme(sampleind(samplespersplit*(s-1)+1:samplespersplit*s)) = s;
%         end
%         
%     end

%   better split data into equal parts
    [split_inds] = getfairsplitsize(numel(myDataset.chunks),splits);
    for s=1:splits
        DataSplitScheme(split_inds(s,1):split_inds(s,2)) = s;
    end
    DataSplitScheme=(ValidSamples > 0).*DataSplitScheme;
    
    %-- CLEANUP: Here comes the samples killer: count the number of samples per
    %class and per crossvalidation. take the minimum number and randomly kill
    %exceeding samples 
    population_size = zeros(classesSpecified,splits);
    for cl=1:classesSpecified
        for s=1:splits
            population_size(cl,s) = numel(find(DataSplitScheme == s & ValidSamples == cl));
            population_ind{cl,s} = find(DataSplitScheme == s & ValidSamples == cl);
        end
    end
    
    population_ind = population_ind
    
    %-- undersampling procedure
    switch splitPurpose
        case 'LOOCV'
            uniq_population_size = unique(population_size);
            if numel(uniq_population_size) > 1
                myDataset.chunksKILLED = zeros(size(DataSplitScheme));
                killpopulations = uniq_population_size(2:end);
                killpopulations(:,2) = killpopulations(:,1) - uniq_population_size(1);
                
                for i=1:size(killpopulations,1)
                    population_oversized=find(population_size == killpopulations(i,1));
                    for j=1:numel(population_oversized)
                        tmp_ind = population_ind{population_oversized(j)};
                        tmp_ind = tmp_ind(randperm(numel(tmp_ind)));
                        %killing happens here
                        killed_ind = tmp_ind(1+end-killpopulations(i,2):end);
                        DataSplitScheme(killed_ind) = 0;
                        myDataset.chunksKILLED(killed_ind) = 1;
                        myDataset.chunks(killed_ind) = 0;
                        tmp_ind(1+end-killpopulations(i,2):end) = [];
                        tmp_ind = sort(tmp_ind);
                        population_ind{population_oversized(j)} = tmp_ind;
                    end
                end
                %allwhatwaskilled=find(myDataset.chunksKILLED)
                new_population_ind = population_ind
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

switch splitPurpose
    case 'LOOCV'
        dataSplitter.splitMatrix = splitMatrix;
        myDataset.dataSplitter = dataSplitter;
        myDataset.dataSplitter.LOOCVDataSplitScheme = DataSplitScheme;
        myDataset.configParameters.dataSplitter = dataSplitter;
        myDataset.configParameters.dataSplitter.LOOCVDataSplitScheme = DataSplitScheme;
    case 'AVERAGING'
        myDataset.dataSplitter.AveragingDataSplitScheme = DataSplitScheme;
        myDataset.configParameters.dataSplitter.AveragingDataSplitScheme = DataSplitScheme;
end
        
        
        









