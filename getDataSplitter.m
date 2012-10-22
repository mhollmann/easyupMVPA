% Returns a defined data splitter used in LOOCV or RFE.
%
% Author: Maurice Hollmann
% Date  : 09/10
%
% Description:
%
%   [dataSplitter] = getDataSplitter(dataset, splitterType [one of: "oneSampleOut", "oneChunkOut", "oddEvenSamples", "oddEvenChunks", "repeatPattern", "custom"], [pattern])
%
%   Using this method the data can be split into "not in use", "training" and "test" parts for a repeated analysis (like LOOCV or RFE).
%   The result is a struct, that describes the whole splitting in a 2D matrix (nmbSplits x nmbSamples).
%   Coding in splitMatrix: 0 = exclude, 1 = test, 2 = train.
%
%   Please consider that the number of splits can be lower than the number of samples.
% 
%   Splitter Types = "oneSampleOut" 
%                    "oneChunkOut" 
%                    "oddEvenSamples"
%                    "oddEvenChunks"
%                    "repeatPattern" 
%                    "custom"
% 
%   There are following possibilities to split the data:
%
%   1. splitterType == 'oneSampleOut'
%
%      This is the typical split for a one sample out cross validation. It creates nmbSamples splits where
%      in every split all samples are training-samples except one. 
%      First split : first sample is test - rest is training
%      Second split: just second sample is test - rest (including first sample) is training   
%      and so on ...
%
%      Example: chunks:       [1 2 1 3 2 3]
%               splitMatrix:   1 2 2 2 2 2
%                              2 1 2 2 2 2
%                              2 2 1 2 2 2
%                              2 2 2 1 2 2
%                              2 2 2 2 1 2
%                              2 2 2 2 2 1
%
%   2. splitterType == 'oneChunkOut'
%
%      This extract a split-pattern that selects all chunks with the same id
%      and iterates the test dataset inside this set. The other samples in the
%      same chunk are not included in the training data.
%
%      Example: chunks:       [1 2 1 3 2 3]
%               splitMatrix:   1 2 1 2 2 2
%                              2 1 2 2 1 2
%                              2 2 2 1 2 1
%                              
%
%   3. splitterType == 'oddEvenSamples'
% 
%      Splits the samples in two parts. 
%      odd sample numbers  = test
%      even sample numbers = training
%
%      Example: Two adjacent samples belong to one run. 
%               Call: getDataSplitter(dataset, 'oddEven')
%
%               samples (nmb):  [1 2 3 4 5 6]
%               splitMatrix:     2 1 2 1 2 1
%
%
%   4. splitterType == 'oddEvenChunks'
% 
%      Splits the samples in two parts. 
%      odd chunk numbers  = test
%      even chunk numbers = training
%
%      Example: Two adjacent samples belong to one run. 
%               Call: getDataSplitter(dataset, 'oddEven')
%
%               chunks:         [1 1 2 3 2 3 3 4 4]
%               splitMatrix:     1 1 2 1 2 1 1 2 2
%
%
%   5. splitterType == 'repeatPattern'
% 
%      May be used for excluding samples from training if they are in the same
%      run as the test sample. In this condition it is possible to exclude samples 
%      e.g. if you dont want to include the third sample of every run use 
%      pattern: [1 1 0 1] (Here one run consists of 4 samples).
%
%      Example: Two adjacent samples belong to one run. 
%               Call: getDataSplitter(dataset, 'repeatPattern', [1 1])
%               chunks:       [1 2 1 2 1 2]
%               splitMatrix:   1 0 2 2 2 2
%                              0 1 2 2 2 2
%                              2 2 1 0 2 2
%                              2 2 0 1 2 2
%                              2 2 2 2 1 0
%                              2 2 2 2 0 1
%
%   6. splitterType == 'custom'
%
%      This option may be used for arbitrary splitting. The given pattern must be 2D with the size
%      N x nmbSamples, where N must be > 0. The pattern has to contain just 0, 1 or 2 (0 = exclude, 1 = test, 2 = train).
%
% Parameters:
%   dataset       - the dataset to print info 
%   splitterType  - string one of: "oneSampleOut" "oneChunkOut" "repeatPattern" "custom"
%   pattern       - the split pattern
%                    in splitterType "repeatPattern" pattern must be in dimensions [1 n] and just contain elements 0 (non-use) and 1 (use)
%                    in splitterType "custom"        pattern must be in dimensions [nmbSamples nmbSamples] and just contain elements 0 (non-used), 1 (test), or 2 (training) 
%
% Returns:
%   dataSplitter  - a struct with the fields: splitterType and splitMatrix
%
% Comments:
%
function [dataSplitter] = getDataSplitter(dataset, splitterType, pattern, balanced)

   dataSplitter = struct('splitterType', '', 'splitMatrix', []);
   
   if( ~exist('dataset','var') || ~exist('splitterType','var') ) 
    error('Usage of getDataSplitter: [dataSplitter] = getDataSplitter(dataset, splitterType [one of: "oneSampleOut", "oneChunkOut", "oddEvenSamples", "oddEvenChunks", "repeatPattern", "custom"], [pattern])');
   end
      
   lengthSet = length(dataset.classIDs);
   
   
   if(~exist('balanced', 'var'))
     if(~exist('pattern', 'var'))
       balanced = false;
     elseif(isscalar(pattern))
       balanced = pattern;
     end
   end
   
   
   %check for splitter type
   if(strcmp(splitterType, 'oneSampleOut'))
     
     dataSplitter.splitterType = 'oneSampleOut';
     dataSplitter.splitMatrix = ones(lengthSet, lengthSet, 'uint8');
     dataSplitter.splitMatrix(:) = 2;
     for i=1:lengthSet
       dataSplitter.splitMatrix(i,i) = 1;
     end
     
   elseif(strcmp(splitterType, 'oneChunkOut'))
     
     dataSplitter.splitterType = 'oneChunkOut';
     
     %get unique elements in chunks
     uniqueElements = unique(dataset.chunks);
     
     dataSplitter.splitMatrix    = zeros(length(uniqueElements), lengthSet, 'uint8');
     dataSplitter.splitMatrix(:) = 2;
     smIndexY = 1;
     
     %loop over unique chunk elements
     for i = 1:length(uniqueElements)
       sel = dataset.chunks == uniqueElements(i);
       dataSplitter.splitMatrix(i, sel) = 1;
     end
     
   elseif(strcmp(splitterType, 'oddEvenSamples'))  
     
     dataSplitter.splitterType  = 'oddEvenSamples';
     dataSplitter.splitMatrix   = ones(1,size(dataset.chunks,2),'uint8');
     
     nmbArr = 1:size(dataset.chunks,2);
     dataSplitter.splitMatrix(~mod(nmbArr,2)) = 2;
     
   elseif(strcmp(splitterType, 'oddEvenChunks'))
     
     dataSplitter.splitterType  = 'oddEvenChunks';
     dataSplitter.splitMatrix   = ones(1,size(dataset.chunks,2),'uint8');
     
     dataSplitter.splitMatrix(~mod(dataset.chunks,2)) = 2;
     
   elseif(strcmp(splitterType, 'repeatPattern'))
     
     if(~exist('pattern','var') || size(pattern,1)>1 || size(find(pattern>1),2)>0 || size(find(pattern<0),2)>0 )
       error('getDataSplitter: for splitter Type "repeatPattern" the argument "pattern" must be defined and in dimensions [1 n] and just contain elements 0 and 1!');
     end
     
     dataSplitter.splitterType = 'repeatPattern';
    
     dataSplitter.splitMatrix    = zeros(lengthSet, lengthSet, 'uint8');
     dataSplitter.splitMatrix(:) = 2;
     smIndexY  = 1;
     lengthPat = length(pattern);
     
     nmbRuns = floor(lengthSet/lengthPat);
     
     %all zero entriess in patternVector will be set to zero in 
     %the resulting splitMatrix
     for i=1:nmbRuns
       indexX = (i-1)*lengthPat+1;
       for j=1:lengthPat
         dataSplitter.splitMatrix(smIndexY,indexX:indexX+lengthPat-1) = 0;
         if(pattern(j) == 1)
           dataSplitter.splitMatrix(smIndexY,indexX+j-1) = 1;
         end
         smIndexY = smIndexY+1;
       end
     end
     
   elseif(strcmp(splitterType, 'custom'))
     
     if(~exist('pattern','var') || size(pattern,2)~=lengthSet || ~isempty(find(pattern>2)) || ~isempty(find(pattern<0)) )
       error('getDataSplitter: for splitter Type "custom" the argument "pattern" must be defined, and in dimensions [N nmbSamples] and just contain elements 0, 1, or 2!)');
     end
     
     dataSplitter.splitterType = 'custom';
     dataSplitter.splitMatrix  = pattern;
     
   else
     error('Usage of getDataSplitter: [dataSplitter] = getDataSplitter(splitterType [one of: "oneSampleOut", "oneChunkOut", "oddEvenSamples", "oddEvenChunks", "repeatPattern", "custom"], [pattern])');
   end
   
   
   %remove all splits that are containing no 1 and 2 in parallel
   smIndexY = 1;
   for i=1:size(dataSplitter.splitMatrix, 1)
     line = dataSplitter.splitMatrix(i,:);
     nmbOnes = size(line(line==1),2);
     nmbTwos = size(line(line==2),2);
     if(nmbOnes>=1 && nmbTwos >= 1)
       resMat(smIndexY, :) = line;
       smIndexY = smIndexY+1;
     end
   end
   dataSplitter.splitMatrix = resMat;
   
   
   
   %now balance the training-samples if intended
   %that means the training set should always have 
   %the same number of samples in each class
   if(balanced)
     
     uniqueClassIDs  = unique(dataset.classIDs);
     if(length(uniqueClassIDs)~=2)
       error('getDataSplitter: Sorry, at the moment the creation of balanced sets is possible just for 2-class problems!');
     end
     c1ID = uniqueClassIDs(1);
     c2ID = uniqueClassIDs(2);
     
     for i=1:size(dataSplitter.splitMatrix, 1)
        smline =  dataSplitter.splitMatrix(i,:);
        %collect all necessary indices
        trainSamples = find(smline==2);
        C1Samples    = find(dataset.classIDs==c1ID);
        C2Samples    = find(dataset.classIDs==c2ID);
        
        intersectTsC1 = intersect(trainSamples,C1Samples);
        intersectTsC2 = intersect(trainSamples,C2Samples);
        
        if(length(intersectTsC1)>length(intersectTsC2))
          %set the randomly picked number of C1's to zero
          randIndices = randperm(length(intersectTsC1));
          smline(intersectTsC1(randIndices(1:(length(intersectTsC1)-length(intersectTsC2))))) = 0;
        elseif(length(intersectTsC2)>length(intersectTsC1))
          %set the randomly picked number of C2's to zero
          randIndices = randperm(length(intersectTsC2));
          smline(intersectTsC2(randIndices(1:(length(intersectTsC2)-length(intersectTsC1))))) = 0;
        end
        
        %check if still enough training samples are available
        if(find(smline==2)<2)
          error('getDataSplitter: Sorry, balancing leads to training sample number lower than 2!');
        end
        
        dataSplitter.splitMatrix(i,:) = smline;
     end %endfor
     
     
   end %end if balanced
     
end