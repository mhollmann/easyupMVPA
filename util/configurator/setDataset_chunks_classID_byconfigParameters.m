% Returns the chunks and classIDs for a dataset specified by the input
% configuratorparameters
%
% Author: Johannes Stelzer
% Date  : 05/11
% Description:

% myDataset.configParameters is assumed to be a structure containing fields
% specifiying the baseDirectory, subjectlist, conditions, unitofdesign,
% removeTrailingTransitions, removeTailTransitions and
% directories (runs) of the functional data. 
% Furthermore it is assumed that a 'condition1.ons' (a vector containing the onsets),
% 'condition1.dur' (a vector containing the durations), 'condition2.ons', 'condition2.dur' ... file in each
% specified directory and obtains the chunks and classID vectors from this.


function [myDataset] = setDataset_chunks_classID_byconfigParameters(myDataset)

% repeat this for each specified run and put it together in the end

chunkindex = 1;
chunks = [];
classIDs = [];

% if unit of design was seconds, round the removetail and removetrailing to
% units of TR
if strfind(myDataset.configParameters.unitofdesign,'econd')
    myDataset.configParameters.removeTrailingTransitions=round(myDataset.configParameters.removeTrailingTransitions/myDataset.configParameters.TR);
    myDataset.configParameters.removeTailTransitions=round(myDataset.configParameters.removeTailTransitions/myDataset.configParameters.TR);
    if isfield(myDataset.configParameters,'sameTrialDurations')
        myDataset.configParameters.sameTrialDurations=round(myDataset.configParameters.sameTrialDurations/myDataset.configParameters.TR);
    end
end
    

for r=1:numel(myDataset.configParameters.runs)
    % get how many volumes are present in this run
    dir_funct = fullfile(myDataset.configParameters.baseDirectory,myDataset.configParameters.subjectname,myDataset.configParameters.runs(r).directory);
    temp_header = load_nii_hdr(fullfile(dir_funct,myDataset.configParameters.dataFileFormat));
    temp_header_dimensions = temp_header.dime.dim;
    volumes = temp_header_dimensions(5);
    
    current_chunks =    zeros(1,volumes);
    current_classIDs =  zeros(1,volumes);
    
    % load the .ons and .dur files
    for c=1:numel(myDataset.configParameters.conditions)
        onsets = load(fullfile(dir_funct,[myDataset.configParameters.conditions{c},'.ons']));
        
        %--load durations if it is not overridden in configurator
        if isfield(myDataset.configParameters,'sameTrialDurations')
            durations = myDataset.configParameters.sameTrialDurations * ones(numel(onsets),1);
        else
            durations = load(fullfile(dir_funct,[myDataset.configParameters.conditions{c},'.dur']));
        end
        
        %convert time to scan number of needed (unit of design)
        if strfind(myDataset.configParameters.unitofdesign,'ec')
            onsets = round(onsets/myDataset.configParameters.TR);
            durations = floor(durations/myDataset.configParameters.TR);
        end
        
        %check whether specifications are integers now
        if sum(rem(onsets,1).^2) > 0
            potentialerrormessage = ['Problem: You need to specify your onsets either in integer numbers (=scans) or in seconds. In the latter case change the variable <unitofdesign> to <seconds>! Problem occuring in ',dir_funct,' at condition ',myDataset.configParameters.conditions{c}];
            error(potentialerrormessage);
        end
        if sum(rem(durations,1).^2) > 0
            potentialerrormessage = ['Problem: You need to specify your durations either in integer numbers (=scans) or in seconds. In the latter case change the variable <unitofdesign> to <seconds>! Problem occuring in ',dir_funct,' at condition ',myDataset.configParameters.conditions{c}];
            error(potentialerrormessage);
        end
        
            
        
        
        if myDataset.configParameters.spm_onsets_specification; onsets = onsets + 1; end
        if myDataset.configParameters.spm_onsets_specification && isfield(myDataset.configParameters,'sameTrialDurations') ==0; durations + 1; end
        % convert them into chunks and classIDs 
        for k=1:numel(onsets)
            % determine a durationscheme for this trial and check if
            % everything's OK
            durationscheme = ones(durations(k),1);
            potentialerrormessage = ['Problem: You are cutting off too much, change your removeTrailingTransitions variable. Problem occuring in ',dir_funct,' at condition ',myDataset.configParameters.conditions{c}];
            
              
            if durations(k) <= myDataset.configParameters.removeTrailingTransitions error(potentialerrormessage); end
            if durations(k) <= myDataset.configParameters.removeTailTransitions error(potentialerrormessage); end
            
            for d=1:myDataset.configParameters.removeTrailingTransitions durationscheme(d) = 0; end
            
            for d=durations(k):-1:(durations(k)-myDataset.configParameters.removeTailTransitions+1) durationscheme(d) = 0; end
            
            %durationscheme
            if sum(durationscheme(:))==0 error(potentialerrormessage); end
            
            % set the classIDs   
            %check that nothing is overwritten and overlapping (also to the
            %end)
            if onsets(k)+durations(k)-1 > numel(current_classIDs) 
                potentialerrormessage = ['Problem: You have defined more scans (onset+duration) than there are in ',dir_funct,' at condition ',myDataset.configParameters.conditions{c}, ', entry nr. ', num2str(k)];
                error(potentialerrormessage);
            end
            
            if sum(current_classIDs(onsets(k):onsets(k)+durations(k)-1)) > 0
                 potentialerrormessage = ['Problem: You have defined an overlapping scan-scheme (at least one scan is now assigned to TWO conditions at once) in ',dir_funct,' at condition ',myDataset.configParameters.conditions{c}, ', entry nr. ', num2str(k)];
                error(potentialerrormessage);
            end
                
            current_classIDs(onsets(k):onsets(k)+durations(k)-1) = ones(durations(k),1)*c; %durationscheme*c
            current_chunks(onsets(k):onsets(k)+durations(k)-1) = durationscheme*chunkindex; %this indexing is temporary to make chunks distinguishable
            chunkindex = chunkindex + 1;
            
        end
    end
    chunks = horzcat(chunks,current_chunks);
    classIDs = horzcat(classIDs,current_classIDs);
end

% reorder chunk naming
unordered_chunks = chunks;
index = 1;
reordered_chunks = zeros(1,numel(unordered_chunks));
newchunknr = 1;
while index <= numel(unordered_chunks)
    nextpos = 1;
    if chunks(index) > 0
        stopinnerloop = 0;
        while stopinnerloop == 0
            oldchunknr = unordered_chunks(index);
            reordered_chunks(index) = newchunknr;
            if index + nextpos <= numel(unordered_chunks)  
                if unordered_chunks(index + nextpos) == oldchunknr
                    reordered_chunks(index + nextpos) = newchunknr;
                    nextpos = nextpos + 1;
                else
                    stopinnerloop = 1;
                    newchunknr = newchunknr + 1;
                end
            else
                stopinnerloop = 1;
            end
        end
        
    end
    index = index + nextpos;
end


myDataset.classIDs = classIDs;
myDataset.chunks = reordered_chunks;

        
    
    
    




