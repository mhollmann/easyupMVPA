function error = dicom2Nifti111 (input_pathname, output_pathname , handles)
%
%
% Quick and dirty dicom2analyse converter
%
% This has the following capabilities:
%   Conversion to analyze of Siemens standard diffusion data
%   Conversion to analyze of MDR non-mosaic diffusion data
%   Conversion to analyze of MDR Mosaic diffusion data
%   Comversion to analyze of EPI data (including Mosaic)
%   Conversion to analyze of pilot images
%   Conversion to analyze of FLASH 3D images
%
% For the diffusion data .bval and .bvec files are automatically created.
%
% It cannot handle
%   Spectroscopy
%   Argus
%   Archive folders from the scanner
% anything complicated
%
% The code is faster (x4) than the Dicom2Ana1xx (1-series) code, but it is not as reliable as
% yet.
%
% Good luck,
%
% Matt (2nd December 2003)

%
% Revision history:
% version 2.00e:  
% added additional diagnostics in the output (including number of files
% found, number of series created, time for headers, time for writing,
% number of errors
%
% version 2.01
%   Changes to correct the problems with rotated and not rotated Mosaic
%   data, I don't really understand these fixes
%
% version 2.02
%   Changes to sort out the rotation stuff again.
%
% version 2.03
%   and again (changes to fix rotation of Mosaic data)
%   added code to enable it to just dump the sequences of interest
% 
% version 2.04
%   updated to output diffusion information (i.e. bvec and bval files) for
%   diffusion acquisitions
%
% version 2.05
%   updated to handle Mosaic Diffusion!!
%
% version 2.06
%   corrected to deal with non-square Mosaic Diffusion (not input_dim and
%   output_dim are separated which makes flipping easier).
%
% version 2.07
%   Fixed to work with new and old versios of the diffusion input data
%
% version 2.08
%   added compatibility with auto_d2a script
%
% version 2.09
% Added extra carriage return to enable compatibility with bedpost
% Corrected bug associated with "handles" not being present (oops)
%
% CONVERTED TO Dicom2Nifti100
%
% version 2.11
% additional "end" removed.  Dopey error.
%
% DICOM2NIFTI 102
%
% 1.03
% Diffusion conversion changed slightly so it doesn't fail with b=0 images
% for diffusion under 25a
%
% 1.03b
% Fix for MATLAB 7 (thanks Dan)
%
% 1.05
% Fix for the qform for the niftii output
%
% 1.06
% Added functionality to handle multiple contrast datasets
%
% 1.07
% Added functionality to CORRECTLY order the images in the output when the
% input order is not the same as the file order
%
% 1.08
% Fixed a bug for empty series
%
% 1.09
% Hacked bug for large number of diffusion directions (i.e. 60)
%
% 1.11
% Fixed bug which failed to clear an important variable (i.e. clear images_in_this_EchoNumbers;)
%
% Structure of the code is as follows:
%
% 1) Use dicom_scan_singlefile on all the files in the directory that don't
% have a '.' in their name
% 2) For each study amd for each series 
%   read in the info using dicom_info ( the standard MATLAB method) for JUST THE FIRST IMAGE
%   read in all the images
%   concatenate an appropriate analyse file
%   Perform mosaic correction as necessary
%
% Finish
%

% Suppress warnings which appear in Matlab V7
warning off MATLAB:intConvertOverflow
warning off MATLAB:intMathOverflow

global delimiter single_file_nifti single_file_nifti_compress double_file_nifti output_mdr_aux_file

single_file_nifti = 0; %This will create a single uncompress .nii file
single_file_nifti_compress = 0; %This will take .nii file and convert to .nii.gz
%single_file_nifti = 0; %turn off the above output
%double_file_nifti = 1; %This will create two a .hdr and .img which are backwardly compatible with analyze
double_file_nifti = 1; %turn off the above output

%output_mdr_aux_file = 1;
output_mdr_aux_file = 0; %Don't output the mdr_aux_file

if (single_file_nifti + double_file_nifti < 1)
    disp('No ouput requested, aborting');
end

%
% This string results in only those series that contain this text being
% dumped
% To dump everything make this string = ''
only_dump_these = ''; %Dump everything
%only_dump_these = 'EPI';
%
error = 0;
error_string = '';

if (ispc)
    %This means that the machine is a PC and presumably the delimiter is a '\' for filenames
    delimiter = '\';
    wrong_delimiter = '/';
else
    if (isunix == 0)
        % this machine is neither Windows or UNIX, RUN AND HIDE
        error_string = 'Machine is neither Windows or UNIX';
        error = 99;
        return
    else
        delimiter = '/';
        wrong_delimiter = '\';
    end
end
delimiter;
nargin;

if (nargin < 1)
    [input_filename, input_pathname] = uigetfile ('*','Pick one of the input files');
end
if ( isequal(exist(input_pathname,'dir'),0) )
    disp('Path does not exist');
    error_string = 'Input path does not exist';
    error = 1;
    return
end
if (isempty(findstr(input_pathname(length(input_pathname)),delimiter)))
    input_pathname = [input_pathname delimiter];
end

if (nargin < 2)
    %previously this used uigetdir, but unfortunately this sometimes hangs
    %on certain disks (we don't know why)
    [output_filename_NOT_USED , output_pathname] = uigetfile ('*','Pick a directory for Analyze output (by selecting a file in that directory)');
    
    if (output_pathname == 0  | exist(output_pathname,'dir')  == 0)
        disp('You must select an output directory. Quitting!');
        return
    end
end

tic
version = 'nifti 111'; 
disp(['Version: ' version]);

filenames_list = dir (input_pathname);


h = waitbar(0,'DICOM2ANA: Reading headers');

%
% Read in the information about each of the files
%
max_series = 0;
actual_file_counter = 0;
for counter = 1 : size(filenames_list,1)
    waitbar(counter/size(filenames_list,1),h);
    occurences_of_IMA = strfind(filenames_list(counter).name , '.IMA');
    occurences_of_ima = strfind(filenames_list(counter).name , '.ima'); %Added to handle data that has its sufix mangled! 210
    occurences_of_dot = strfind(filenames_list(counter).name , '.');
    file_directory = [input_pathname filenames_list(counter).name];
    
    if ((length(occurences_of_dot)==0  | length(occurences_of_IMA)~=0 | length(occurences_of_ima)~=0) & (exist(file_directory,'file') ~= 0 & exist(file_directory,'dir') == 0)) %only consider files that have no '.' in the  
        actual_file_counter = actual_file_counter + 1;
        % the filename ends in IMA
        
        %file_directory
        info_tmp = dicom_scan_singlefile(file_directory);
        filename_data(actual_file_counter).info = info_tmp;
        filename_data_series(actual_file_counter) = info_tmp.series;
        filename_data_image(actual_file_counter) = info_tmp.image;
        filename_data_EchoNumbers(actual_file_counter) = info_tmp.EchoNumbers;
        
        filename_data(actual_file_counter).filename = filenames_list(counter).name;
        
        max_series = max(max_series,info_tmp.series);
        
        store_file_directory = file_directory;   
    end 
    if (exist('handles' ,'var'))    
        pause(get(handles.speed_slider,'Value')); 
    end
end

disp(['Number of image files found  ' num2str(actual_file_counter)]);
disp(['Number of series files found ' num2str(max_series)]);

%Create a directory with the patient_name
% 'Patient_' PatientName
try
    info = dicominfo(store_file_directory);
    corrected_name = legitimise_filename (['Patient_' info.PatientID ]);
    if (exist([output_pathname delimiter corrected_name ],'dir') == 0)      
        mkdir(output_pathname,corrected_name);
    end
    output_pathname = [output_pathname delimiter corrected_name ];
catch
    disp('Failed to create specific Patient output directory');
end

close (h)
h = waitbar(0,'DICOM2ANA: writing series');

error_counter = 0;

for series = 1 : max_series
    waitbar((series-0.5)/max_series,h);
    
    images_in_this_series = find(filename_data_series == series);
    
    maximum_EchoNumbers = max(filename_data_EchoNumbers(images_in_this_series));
    
    [values,index] = sort(filename_data_image(images_in_this_series));
    
    for EchoNumbers_counter = 0 : maximum_EchoNumbers
        index_counter = 0;
        
        for index_dump_loop = 1:length(index)
            %Is this the right echo?
            if (filename_data_EchoNumbers(index(index_dump_loop)) == EchoNumbers_counter)
                index_counter = index_counter + 1;
                images_in_this_EchoNumbers(index_counter) = images_in_this_series(index(index_dump_loop));
            end
        end
        %images_in_this_EchoNumbers = images_in_this_series(find(filename_data_EchoNumbers(images_in_this_series(index)) == EchoNumbers_counter));
        
        if (index_counter < 1)
            %Empty series
            
        else
            images_to_dump = images_in_this_EchoNumbers; %images_in_this_series(index);
            clear images_in_this_EchoNumbers;
            clear images_in_this_EchoNumbers;
            
            %Use dicominfo to find out the Study description
            info = dicominfo([input_pathname filename_data(images_to_dump(1)).filename]);
            if (exist('handles' ,'var'))    
                pause(get(handles.speed_slider,'Value')); 
            end
            %if the text appears in the ProtocolName or there is no text to
            %match
            dump = 0;
            if ( isempty(only_dump_these) )
                dump = 1;
            else
                if ( ~isempty(strfind( info.ProtocolName, only_dump_these )) ) 
                    dump = 1;
                end
            end
            if (dump == 1)
                %try
                if EchoNumbers_counter > 0
                     output_filename = [output_pathname delimiter 'images_' num2str(series) '_echo_' num2str(EchoNumbers_counter)];
                else
                    output_filename = [output_pathname delimiter 'images_' num2str(series)];
            end
                if (exist('handles' ,'var')) 
                    M2A (filename_data , output_filename , images_to_dump, input_pathname , handles);
                else
                    M2A (filename_data , output_filename , images_to_dump, input_pathname );
                end
                
                %    dumpoutthe indexed filenames
                %catch
                %error_counter = error_counter + 1;
                %disp(['Failed for series ' num2str(series)]);
                %end
            end
        end
    end
end

if (error_counter > 0)
    error = 2;
    error_string = 'Conversion failed for one or more series';
end
disp(['Total writing errors         ' num2str(error_counter)]);
close(h)
toc_time = toc;
disp(['Total time for conversion    ' num2str(toc_time)]);


function nothing = M2A (filename_data , output_filename , image , input_pathname , handles)
%
% Write a set of images to an analyze file
%
% function VO = M2A (filename_data , output_filename , image )
%
% This takes a range of data (i.e.   M2A(filename_data, 1:2 , 1:128) ) and
% bungs it into an Analyze file format.  If the data are not of a consistent size then
% this probably will work but the output will be ugly!
%
% dicom_str = this is not the full dicom_structure but just the
% filename_data as created previously
% output_filename = filename for the output data (note that there will be a 'output_filename'.hdr & an 'output_filename.img'
%
% image     = as above (optional argument)
%
%    flip_flag = 0; % Don't flip            (dimensions remain the same)
%    flip_flag = 1; % flip diagonally       (dimensions flip)
%    flip_flag = 2; % rotate by 90 degrees  (dimensions flip)
%    flip_flag = 3; % rotate by 180 degrees (dimensions remain the same)
%    flip_flag = 4; % rotate by 270 degrees (dimensions flip)
%    flip_flag = 5; % flip up down          (dimensions remain the same)
%    flip_flag = 6; % flip left right       (dimensions remain the same)
%
% Note: if the image argument is not included, all the images in that series will be used
%
% This code can also handle MOSAIC data
%
% MDR 1/04/2002
%
global delimiter single_file_nifti single_file_nifti_compress double_file_nifti output_mdr_aux_file DIFFUSION_FLAG
nothing = 0;
flip_flag = 4;%4; %THis seems to be the default required for conversion from DICOM to ANALYZE

%Read in the first file completely using the MATLAB DicomInfo
info_tmp = dicominfo([input_pathname delimiter filename_data(image(1)).filename]); %Read all the information as then we can have an easy ride for pixel dimension etc.

%Append info_tmp.ProtcolName to the output_filename
try
    output_filename = [output_filename '_' info_tmp.ProtocolName];
catch
    %Sometimes the ProtocolName doesn't exist.  In this case don't include
    %it!!!
end

% Open the image file
fname = strcat(output_filename , '.img');
fp    = fopen(fname,'r+');
if fp == -1,
    fp = fopen(fname,'w');
    if fp == -1,
        VO = [];
        disp('Can''t open image file');
        return;
    end;
end;

% Note this main section works with MOSAIC, the special case are non-MOSAIC where the Mosaic factor is 1x1

% Need to determine whether the data are Mosaic, and if so what the
% dimensions are
% This is a MOSAIC dataset this will be handled later

%For the time-being lets assume non-mosaic
% Here we will add code that finds out about the MOSAIC format
Mos_factor = 1;
MOSAIC_FLAG = 0;
DIFFUSION_FLAG = 0;

if (strcmp(info_tmp.ImageType , 'ORIGINAL\PRIMARY\M\MOSAIC'))
    MOSAIC_FLAG = 1;
end

if (strcmp(info_tmp.ImageType , 'ORIGINAL\PRIMARY\M\ND\MOSAIC')) 
    MOSAIC_FLAG = 1;
end
if (~isempty(strmatch('Private_0029_1020' , fieldnames(info_tmp))))
    if (~isempty(findstr('<ParamFunctor."MosaicUnwrapper">', char(info_tmp.Private_0029_1020)')))
        MOSAIC_FLAG = 1;
    end
end

if (~isempty(strmatch('Private_0029_1020' , fieldnames(info_tmp))))
    if (~isempty(findstr('sDiffusion.lDiffWeightings', char(info_tmp.Private_0029_1020)')))
        DIFFUSION_FLAG = 1;
    end
end

if (MOSAIC_FLAG == 1)
    %The following code is required to handle non-axial scans
    if (info_tmp.AcquisitionMatrix(1)  ~= 0)
        Mos_factor = info_tmp.Columns   / info_tmp.AcquisitionMatrix(1);
    else
        if (info_tmp.AcquisitionMatrix(2)  ~= 0)
            Mos_factor = info_tmp.Columns   / info_tmp.AcquisitionMatrix(2);
        else
            Mos_factor = info_tmp.Columns   / info_tmp.AcquisitionMatrix(3);
        end
    end
end

% So, the first actual data is at "tmp_counter"
% The format of the MOSAIC data is such that:
% if there are N images in the mosaic
% the data will be at N-1, (2N)-1, (3N)-1, etc.
%
% The Mosaic file seems also always to be square
%Mos_factor = sqrt(tmp_counter+1);
%fprintf(['Mosaic factor ' num2str(Mos_factor) 'x' num2str(Mos_factor) '\n']);
flag = 0;
counter = 1;

for L = image
    
    if (flag == 0)
        
        %This is the first time through here so set up some things that
        %don't need changing
        dim = zeros(8,1);
        
        input_dim(1) = 4; % number of dimensions
        input_dim(2) = info_tmp.Rows / Mos_factor;
        input_dim(3) = info_tmp.Columns / Mos_factor;
        
        if (DIFFUSION_FLAG == 1)
            
            hmm = char(info_tmp.Private_0029_1020');
            
            [slice_array_size, success] = values_in_strings(hmm,'sSliceArray.lSize');
            if (success >= 1)
                input_dim(4) = eval(slice_array_size);
            else
                disp(['problem with finding slice array size in diffusion dicom header']);
            end
            
            [number_of_diffusion_weightings, success] = values_in_strings(hmm,'sDiffusion.lDiffWeightings');
            
            if (success >= 1)
                [number_of_diffusion_directions, success] = values_in_strings(hmm,'sDiffusion.lDiffDirections');
                if (success >= 1)
                    total_volumes = 1 + (eval(number_of_diffusion_directions)*(eval(number_of_diffusion_weightings)-1));
                    diffusion.weightings = eval(number_of_diffusion_weightings);
                    diffusion.directions = eval(number_of_diffusion_directions);
                else
                    number_of_diffusion_directions = '1';
                    total_volumes = 1 + (eval(number_of_diffusion_directions)*(eval(number_of_diffusion_weightings)-1));
                    diffusion.weightings = eval(number_of_diffusion_weightings);
                    diffusion.directions = eval(number_of_diffusion_directions);
                    %disp(['problem with finding number of diffusion encoding directions in diffusion dicom header']);
                    %total_volumes = 0;
                end
            else
                disp(['problem with finding number of diffusion weightings (i.e. 2 for b=0,1000 in diffusion dicom header']);
                total_volumes = 0;
            end  
            
            if (Mos_factor == 1)
                
                % This is the Siemens diffusion sequence
                input_dim(5) = length(image)/input_dim(4);
            else
                %This is my diffusion sequence
                %disp(['This is my diffusion sequence']);
                
                input_dim(5) = length(image);
            end
            if (input_dim(5) ~= total_volumes)
                disp(['WARNING: Total_volumes according to diffusion information ' num2str(total_volumes) ' , does not match the number of images divided by the number of slices' num2str(input_dim(4))]);
            end              
            try
                %
                % We can create the b-val data now.
                b_val_counter = 1;
                diffusion.b_val(b_val_counter) = 0; %This is always the first value (I think), andhas only one direction
                for b_weights = 2:eval(number_of_diffusion_weightings)
                    %
                    % What is the n'th diffusion weighting
                    [value, success] = values_in_strings(hmm,['sDiffusion.alBValue[' num2str(b_weights-1) ']']);
                    if (success ~= 1)
                        disp(['Problem reading the diffusion weightings']);
                        crash_out;    
                    end
                    for b_dirns = 1 : eval(number_of_diffusion_directions)
                        b_val_counter = b_val_counter + 1;
                        
                        diffusion.b_val(b_val_counter) = eval(value);
                    end
                end                   
            catch
                disp(['No diffusion b-val results']);
            end            
            try    
                [value, success] = values_in_strings(hmm,['sDiffusion.lNoiseLevel']);
                mdr_lBerhens_flag = eval(value);
            catch
                mdr_lBerhens_flag = 1;
            end
            
            
            %Make the b-vec and b-val files
            Make_bvecs103(output_filename, info_tmp, mdr_lBerhens_flag, diffusion);               
            
        else
            % NOT diffusion == 1
            
            if (Mos_factor == 1)
                input_dim(4) = length(image); %Ironically this is correct!
                input_dim(5) = 1; %This has been assumed but we probably need to be cleverer about input_dim(4) and input_dim(5)
            else   
                
                %input_dim(4) = Mos_factor*Mos_factor; %Ironically this is correct! But doesn't take account for the extra blank images added by Siemens
                input_dim(5) = length(image);%)/(Mos_factor*Mos_factor); %This has been assumed but we probably need to be cleverer about input_dim(4) and input_dim(5)
            end
        end % End of Diffusion == 1 if/else/end
        
        input_dim(6) = 1; % only 4D therefore irrelevant
        input_dim(7) = 1; % only 4D therefore irrelevant
        input_dim(8) = 1; % only 4D therefore irrelevant
        
        pixdim = zeros(8,1);
        pixdim(1) = 0;
        pixdim(2) = info_tmp.PixelSpacing(1); % these may be the wrong way around
        pixdim(3) = info_tmp.PixelSpacing(2); %these may be the wrong way around
        %I can only get the z_pixel spacing from the interslice distance (how do I do this for Mosaic?)
        pixdim(4) = 1;% I have no idea what this should be, but .....dicom_str{i}.Study{j}.Series{k}.Image{l}.Whole_Structure.SpacingBetweenSlices;
        if (Mos_factor >= 2  |  strcmp(info_tmp.ScanningSequence,'SE\EP'))
            %This EPI, so set the time scaling to the TR
            pixdim(5) = info_tmp.RepetitionTime; %This is arbitrary but saves working outthe real numbers which are tricky
            if (~isempty(strmatch('Private_0029_1020' , fieldnames(info_tmp))))
                if (~isempty(findstr('lDelayTimeInTR', char(info_tmp.Private_0029_1020)')))
                    %See if there is an additional time delay after the images,
                    %look in 
                    hmm = char(info_tmp.Private_0029_1020');
                    xtra_tr_delay = eval(values_in_strings(hmm,'lDelayTimeInTR'));
                    pixdim(5) = pixdim(5) + xtra_tr_delay/1000; %Convert from us to ms
                end
            end
        else
            %This isn't EPI, so calculating the time between images is
            %quite tricky
            if ( strmatch(info_tmp.MRAcquisitionType, '2D') )
                pixdim(5) = info_tmp.RepetitionTime * info_tmp.NumberOfPhaseEncodingSteps; %assume slice selective, and not 3D.  If anyone is doint a series of 3D's then this will be wrong. MDR
            end
            if ( strmatch(info_tmp.MRAcquisitionType, '3D') )
                pixdim(5) = info_tmp.RepetitionTime * info_tmp.NumberOfPhaseEncodingSteps * input_dim(4);
            end
        end
        
        pixdim(6) = 1; % only 4D therefore irrelevant
        pixdim(7) = 1; % only 4D therefore irrelevant
        pixdim(8) = 1; % only 4D therefore irrelevant
        
        A = zeros(input_dim(2)*input_dim(3),1);
        %flip_flag = 0; % Don't flip            (dimensions remain the same)
        %flip_flag = 1; % flip diagonally       (dimensions flip)
        %flip_flag = 2; % rotate by 90 degrees  (dimensions flip)
        %flip_flag = 3; % rotate by 180 degrees (dimensions remain the same)
        %flip_flag = 4; % rotate by 270 degrees (dimensions flip)
        %flip_flag = 5; % flip up down       (dimensions remain the same)
        %flip_flag = 6; % flip left right    (dimensions remain the same)
        
        output_dim = input_dim; % This is the default condition
        if (flip_flag == 1 | flip_flag == 2 | flip_flag == 4)
            %Swap the dimensions 
            output_dim(2) = input_dim(3);
            output_dim(3) = input_dim(2);
            % And the voxel sizes
            temp = pixdim(2);
            pixdim(2) = pixdim(3);
            pixdim(3) = temp;
        end
    end %flag == 0
    
    
    
    % Open the dicomfile
    
    tmp_image = dicomread([input_pathname delimiter filename_data(L).filename]);
    if (exist('handles' ,'var'))    
        pause(get(handles.speed_slider,'Value')); 
    end
    
    %Split the data into blocks before putting it into the A matrix
    
    
    
    if (flag == 0) 
        if (Mos_factor ~= 1)
            %If it is mosaic then this variable may be stored as
            %SpacingBetweenSlices (in the info for these data)
            try
                pixdim(4) = info_tmp.SpacingBetweenSlices;
            catch
                pixdim(4) = 0;
            end%time_of_first_slice = str2num(info_tmp.AcquisitionTime);
        else
            %We can determine this much better by comparing the centre of the
            %slices for the first and second images in the series
            %So, this is the first iteration around the loop, so store the
            %centre of the slice, for the next iteration
            z_slice_posn = info_tmp.ImagePositionPatient;
        end    
    end
    if (flag == 1)
        if (Mos_factor == 1)
            %This is the second time around this loop so work out the slice
            %spacing by differencing the positions of the 1st and second slices
            
            %This requires a second call to dicominfo (which is slow) but
            %thats life
            info_tmp2 = dicominfo([input_pathname delimiter filename_data(L).filename]);
            pixdim(4) = magnitude(info_tmp2.ImagePositionPatient - info_tmp.ImagePositionPatient);
        end
    end    
    
    %if (flag == 2)
    %info = dicominfo(  char(dicom_str{i}.Study{j}.Series{k}.Image{l}.Filename));
    %         pixdim(5) = str2num(info.AcquisitionTime) - time_of_first_slice;
    %         
    %         if (pixdim(5) < 0)
    %             %Data acquired at mid-night or midday, may get buggered
    %             pixdim(5) = pixdim(5) +  120000.0;
    %         end
    %         if (pixdim(5) > 4000.0)
    %             %it is likely that this has gone over a "hour" boundary
    %             % Time is stored as a character array of
    %             pixdim(5) = pixdim(5) - 4000.0;
    %         end
    %         if (pixdim(5) > 40.0)
    %             %it is likely that this has gone over a "minute" boundary
    %             % Time is stored as a character array of
    %             pixdim(5) = pixdim(5) - 40.0;
    %         end
    %     end
    counter = 1;
    image_output_counter = 0;
    for x = 1 : Mos_factor
        for y = 1 : Mos_factor
            
            if (Mos_factor ~= 1)
                try
                    sub_image = tmp_image((1+(input_dim(2)*(x-1))):(input_dim(2)*x),(1+(input_dim(3)*(y-1))):(input_dim(3)*y));           
                catch
                    disp(['Wierd orientation for Mosaic data']);
                    sub_image = tmp_image((1+(input_dim(3)*(y-1))):(input_dim(3)*y),(1+(input_dim(2)*(x-1))):(input_dim(2)*x));
                end
            else           
                sub_image = tmp_image;%((1:input_dim(3)),(1:input_dim(2)));
            end    
            %flip_flag = 0; % Don't flip            (dimensions remain the same)
            %flip_flag = 1; % flip diagonally       (dimensions flip)
            if (flip_flag == 1)
                sub_image = (rot90(flipdim(sub_image,2),1)); 
            end
            %flip_flag = 2; % rotate by 90 degrees  (dimensions flip)
            if (flip_flag == 2)
                sub_image = rot90(sub_image,1);  
            end
            %flip_flag = 3; % rotate by 180 degrees (dimensions remain the same)
            if (flip_flag == 3)
                sub_image = rot90(sub_image,2);    
            end
            %flip_flag = 4; % rotate by 270 degrees (dimensions flip)
            if (flip_flag == 4)
                sub_image = rot90(sub_image,3);  
            end
            %flip_flag = 5; % flip up down       (dimensions remain the same)
            if (flip_flag == 5)
                sub_image = flipdim(sub_image,1); 
            end
            %flip_flag = 6; % flip left right    (dimensions remain the same)
            if (flip_flag == 6)
                sub_image = flipdim(sub_image,2); 
            end
            
            if (Mos_factor ~= 1)
                maximum_pixel = max(max(sub_image));
                if (maximum_pixel == 0)
                    % this image is one of the padded images at the end of the mosaic set, so therefore don't output it
                    break;  
                else
                    try
                    A(counter: counter+(output_dim(2)*output_dim(3))-1) = reshape(sub_image,output_dim(2)*output_dim(3),1);
                    counter = counter + output_dim(2)*output_dim(3);
                catch
                    disp('hmm');
                end
                    image_output_counter = image_output_counter + 1;
                end
            else
                try
                A(counter: counter+(output_dim(2)*output_dim(3))-1) = reshape(sub_image,output_dim(2)*output_dim(3),1);
                counter = counter + output_dim(2)*output_dim(3);
                catch
                    disp('hmm');
                end
                image_output_counter = image_output_counter + 1;
            end
            
        end
    end
    if (Mos_factor == 1)
        %In this case dim(4) already contains the correct number
    else
        output_dim(4) = image_output_counter;
    end
    
    analyze_data_type = 'short';
    pause(0.00001);
    fwrite(fp,A,analyze_data_type);
    flag = flag + 1;
end

if (flag == 1 & Mos_factor == 1) 
    % Only one image in the dataset
    output_dim(4) = 1;
    pixdim(4) = 1;
end

fclose(fp);



%Now output a header for this data
header_filename                 = [output_filename '.hdr'];

% For byte swapped data-types, also swap the bytes around in the headers.
mach = 'native';
%if spm_type(TYPE,'swapped'),
%       if spm_platform('bigend'),
%               mach = 'ieee-le';
%       else,
%               mach = 'ieee-be';
%       end;
%       TYPE = spm_type(spm_type(TYPE));
%end;
fid             = fopen(header_filename,'w',mach);

if (fid == -1),
    error(['Error opening ' header_filename '. Check that you have write permission.']);
end;
%---------------------------------------------------------------------------
data_type       = ['dsr      ' 0]; %Leave it as this!


% set header variables
%---------------------------------------------------------------------------
%DIM            = DIM(:)'; if size(DIM,2) < 4; DIM = [DIM 1]; end  %This is the dimensions of the voxel data from SPM
%VOX            = VOX(:)'; if size(VOX,2) < 4; VOX = [VOX 0]; end  %This is the Voxel data in SPM
%output_dim            = [4 DIM(1:4) 0 0 0];   
%pixdim         = [0 VOX(1:4) 0 0 0];
vox_offset      = 0; %This is overwritten if we are doing a .nii file
funused1        = 1;  % spm uses SCALE;
glmax           = 1;
glmin           = 0;
bitpix          = 0;
descrip         = zeros(1,80);
my_text = 'Siemens Dicom to nifti-1 version 100. MDR';
descrip(1:length(my_text)) = my_text;

%The following removes the pathname from the filename, hopefully bringing it down to 24characters!
delims = strfind(output_filename,delimiter);


mdr_aux_file = [ output_filename( delims(length(delims))+1 : length(output_filename))  '.mdr'  ];
aux_output_pathname = [output_filename(1:delims(length(delims)))];
if (length(mdr_aux_file) > 24)
    mdr_aux_file = [mdr_aux_file(1:19) '.mdr' 0];
else
    mdr_aux_file = [mdr_aux_file '                            '];
    if (output_mdr_aux_file == 1)
        mdr_aux_file        = [mdr_aux_file(1:23) 0];
    end
end

diff_aux_file = ['                                            '];
diff_aux_file        = [diff_aux_file(1:23) 0];

if (DIFFUSION_FLAG == 1)
    diff_aux_file = [ output_filename( delims(length(delims))+1 : length(output_filename))  '.bmat'  ];
    aux_output_pathname = [output_filename(1:delims(length(delims)))];
    if (length(diff_aux_file) > 24)
        diff_aux_file = [diff_aux_file(1:18) '.bmat' 0];
    else
        diff_aux_file = [diff_aux_file '                            '];
        if (output_diff_aux_file == 1)
            diff_aux_file        = [diff_aux_file(1:23) 0];
        end
    end
    
end

aux_file = diff_aux_file;

origin          = [0 0 0 0 0];

TYPE = 4; %This seems to be the format of Siemens Data
%---------------------------------------------------------------------------
if TYPE == 1;   bitpix = 1;  glmax = 1;        glmin = 0;       end
if TYPE == 2;   bitpix = 8;  glmax = 255;      glmin = 0;       end
if TYPE == 4;   bitpix = 16; glmax = 32767;    glmin = 0;       end
if TYPE == 8;   bitpix = 32; glmax = (2^31-1); glmin = 0;       end
if TYPE == 16;  bitpix = 32; glmax = 1;        glmin = 0;       end
if TYPE == 64;  bitpix = 64; glmax = 1;        glmin = 0;       end

%---------------------------------------------------------------------------

data_type = 'OCMR DCMn1';

tmp_db_name = [filename_data(image(1)).info.idStr '_'  filename_data(image(1)).info.patientStr '                        '];
db_name   = tmp_db_name(1:18);


fseek(fid,0,'bof');

% write (struct) header_key
%---------------------------------------------------------------------------
fwrite(fid,348,         'int32');       %size of header
fwrite(fid,data_type,   'char' );       %data_type[10]
fwrite(fid,db_name,     'char' );       %db_name[18]
fwrite(fid,0,           'int32');       %extents
fwrite(fid,0,           'int16');       %session_error
fwrite(fid,'r',         'char' );       %regular


%DIM_INFO
if isfield(info_tmp,'InPlanePhaseEncodingDirection') % Matlab V.7 renamed field...
    info_tmp.PhaseEncodingDirection = info_tmp.InPlanePhaseEncodingDirection;
end
if (strcmp(info_tmp.PhaseEncodingDirection, 'COL'))
    read_direction = 1;
    phase_direction = 2;
else % 'ROW'
    read_direction = 2;
    phase_direction = 1; 
end


if (output_dim(read_direction+1) < output_dim(phase_direction+1))
    disp(['Possible miss-labelling of phase-encoding direction in Nifti file']);
end
slice_direction = 3;

%Calculate the quaternion
%
first_vector = info_tmp.ImageOrientationPatient(1:3);% 3 floats for first direction
second_vector  = info_tmp.ImageOrientationPatient(4:6);       % 3 floats for second direction
slice_vector = cross(first_vector,second_vector); %I probably need to sort this out!!! MDR

orientation_matrix = [first_vector'; second_vector'; slice_vector'];
pixel_scaling = zeros(3,3); 
pixel_scaling(1,1) = pixdim(1);
pixel_scaling(2,2) = pixdim(2); 
pixel_scaling(3,3) = pixdim(3);
orientation_matrix_scaled = orientation_matrix * pixel_scaling;

pixel_origin = info_tmp.ImagePositionPatient;       % 3 floats, location in x,y,z
offset = - orientation_matrix_scaled * pixel_origin;

R = zeros(4,4);

R(1:3,1:3) = orientation_matrix;

R(1,4) = offset(1);
R(2,4) = offset(2);
R(3,4) = offset(3);

R(4,4) = 1;
[qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac] =  mat44_to_quatern( R );

pixdim(1) =  qfac;%this is a dreadful bodge


dim_info = read_direction + phase_direction*4 + slice_direction*16;  %This is the wierd format chosen by nifti
fwrite(fid,dim_info,         'char' );       %hkey_un0

% write (struct) image_dimension
%---------------------------------------------------------------------------

%fseek(fid,40,'bof');

%disp(pixdim');

fwrite(fid,output_dim,         'int16');       %output_dim[8]
fwrite(fid,0,           'float');            %intent_p1
fwrite(fid,0,           'float');            %intent_p2
fwrite(fid,0,           'float');            %intent_p3
fwrite(fid,0,           'int16');           %intent_code
%fwrite(fid,'mm',        'char' );       %unused8
%fwrite(fid,0,           'char' );       %unused9
%fwrite(fid,0,           'char' );       %unused9
%
%fwrite(fid,zeros(1,8),  'char' );   %unused 10 , 11 , 12 , 13 
%fwrite(fid,0,           'int16');       %unused14
fwrite(fid,TYPE,        'int16');       %datatype
fwrite(fid,bitpix,      'int16');       %bitpix
fwrite(fid,0,           'int16');       %dim_un0
fwrite(fid,pixdim,      'float');       %pix_dim[8], width, height, thickness etc.
fwrite(fid,vox_offset,  'float');       %vox_offset
scl_slope = 1.0;                        %Effectively a unit transform
fwrite(fid,scl_slope,    'float');      %funused1 SPM uses this for SCALE, we can copy this
scl_inter = 0.0;                         %Effectively a unit transform
fwrite(fid,scl_inter,    'float');       %funused2

fwrite(fid,0,           'short');       %slice_end
fwrite(fid,0,           'char');         %slice_code
fwrite(fid, 18 ,        'char');        %xyz_t units (18 = 2 + 16) 2 = mm, 16=ms

fwrite(fid,0,           'float');       %cal_max (float)
fwrite(fid,0,           'float');       %cal_min   (float)
fwrite(fid,0,           'int32');       %compressed (float)
fwrite(fid,0,           'int32');       %verified (float)
fwrite(fid,glmax,       'int32');       %glmax
fwrite(fid,glmin,       'int32');       %glmin

% write (struct) image_dimension
%---------------------------------------------------------------------------
fwrite(fid,descrip,     'char');
fwrite(fid,aux_file,    'char');


% NIFTI INFORMATION FROM http://nifti.nimh.nih.gov/dfwg/src/nifti1.h
if (flip_flag == 4  )
    NIFTI_XFORM_SCANNER_ANAT = -1;  %Rotation through 270 degrees means we need a q-form of -1
else
    NIFTI_XFORM_SCANNER_ANAT = 1;
end
qform_code=NIFTI_XFORM_SCANNER_ANAT;
fwrite(fid,qform_code, 'int16');

sform_code=0;
fwrite(fid,sform_code, 'int16');

fwrite(fid,qb,           'float'); 
fwrite(fid,qc,           'float'); 
fwrite(fid,qd,           'float'); 
fwrite(fid,qx,           'float'); 
fwrite(fid,qy,           'float'); 
fwrite(fid,qz,           'float'); 

fwrite(fid,zeros(1,12), 'float');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');

fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');

fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');

fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');
fwrite(fid,0, 'char');

%This is the magic information, if it is 2 files then have
fwrite(fid,'ni1');
fwrite(fid,0,'char');
%Otherwise have
%fwrite(fid,'n+1\0');

fclose(fid);
%spm_unlink(P);
%error(['Error writing ' P '. Check your disk space.']);
%end

%s   = ftell(fid);
%fclose(fid);

if (single_file_nifti == 1)
    
    
    header_filename                 = [output_filename '.hdr'];
    fid_header       = fopen(header_filename ,'r' , mach);
    header_read = fread(fid_header);
    fclose(fid_header);
    
    
    %Write the nifi file
    single_file_nifti_filename                 = [output_filename '.nii'];    
    fid             = fopen(single_file_nifti_filename,'w',mach);
    
    %Overrigtht the vox-offset values
    fwrite(fid , header_read(1:108));
    fwrite(fid , 352 , 'float');  %The data starts at 352 now, we add 4 bytes for some reason!
    fwrite(fid, header_read(113:344));
    
    %The last 4 bytes are different
    fwrite(fid,'n+1');
    fwrite(fid,0,'char');
    
    %And we now need an extra 4 bytes (for some reason)
    fwrite(fid , 0, 'char');
    fwrite(fid , 0, 'char');
    fwrite(fid , 0, 'char');
    fwrite(fid , 0, 'char');
    
    
    %Now append the data
    fid_data = fopen(fname,'r');
    data_read = fread(fid_data);
    fclose(fid_data);
    
    fwrite(fid,data_read);
    fclose(fid);
end

if (double_file_nifti == 1)
    %Retain both the img and hdr files
else
    disp('Deleting .hdr and .img files');
    delete(fname);
    delete(header_filename);
end

if single_file_nifti_compress
    tmp = ['! gzip ' single_file_nifti_filename];
    eval(tmp)
end

%
%THIS IS THE OLD AUX file that includes orientation information (that
%wasn't included in analyze but is now included in nifti)
if (output_mdr_aux_file == 1)
    % %Now write out the auxilliary file.  This should contain information on slice locations and orientations
    fid             = fopen([aux_output_pathname mdr_aux_file],'w',mach);
    
    if (fid == -1),
        error(['Error opening ' [aux_output_pathname mdr_aux_file] '. Check that you have write permission.']);
    end;
    % 
    % 
    pixel_origin = info_tmp.ImagePositionPatient;       % 3 floats
    % %Next the vector describing how to get from one slice to the next
    slice_vector = [0.0 , 0.0 , 0.0 ]; %this will be created when flag == 1 (when we get to the next slice)
    % %Next the in-plane vector for the data as it is in the Dicom file
    first_vector = info_tmp.ImageOrientationPatient(1:3);% 3 floats for read direction
    second_vector  = info_tmp.ImageOrientationPatient(4:6);       % 3 floats for phase encode direction
    % %This describes the orientation of the patient within the magnet
    Patient_Position = info_tmp.PatientPosition;
    % %First the slice offset information, for the first slice
    fwrite(fid,pixel_origin,                'float');       % 3 floats
    % %Next the vector describing how to get from one slice to the next
    fwrite(fid,slice_vector,                'float');       % 3 floats
    % %Next the in-plane vector for the data as it is in the Dicom file
    fwrite(fid,first_vector,         'float');       % 3 floats for read direction
    fwrite(fid,second_vector,                'float');       % 3 floats for phase encode direction
    % %Next include the patient orientation information i.e.
    % % HFS (head first suppine), HFP (head first prone), FFS (feet first supine) etc.
    fwrite(fid,Patient_Position, 'char'); %3 char
    % %Finally the information concerning whether the images have been flipped in the output file
    fwrite(fid,flip_flag,           'int16');       % 1 int (as described above)
    fclose(fid);
end

if (DIFFUSION_FLAG == 1)
    %Read in the bval and bvec files
    %fid_bvecs = fopen([ output_filename '.bvecs'],'r');
    bvecs = load ([ output_filename '.bvecs']);
    
    %bvecs = eval(output_filename);  %fread(fid_bvecs,'float');
    %fclose(fid_bvecs);
    
    
    %create the b r r(transpose)  matrix
    %fid_bvals = fopen([ output_filename '.bvals'],'r');
    bvals = load ([output_filename '.bvals']);
    %bvals = eval(output_filename); %fread(fid_bvals,'float');
    %fclose(fid_bvals);
    %output the xx xy xz yy yz zz elements of the matrix
    % %Now write out the auxilliary file.  This should contain information on slice locations and orientations
    fid             = fopen([aux_output_pathname aux_file],'w',mach);
    if (fid == -1),
        error(['Error opening ' [aux_output_pathname aux_file] '. Check that you have write permission.']);
    end
    
    for counter  = 1 : length(bvals)
        %
        % Determine the bvec vector
        %Determine bval * bvec bvec'
        normalisation = bvecs(:,counter)'*bvecs(:,counter);
        if (normalisation <= 0.00000001)
            normalisation = 1;
        end
        bmat = bvals(counter)*bvecs(:,counter)*bvecs(:,counter)'/normalisation;
        %Choose the unique components.
        fprintf(fid,'%f %f %f %f %f %f \n', bmat(1,1),bmat(1,2), bmat(1,3), bmat(2,2), bmat(2,3), bmat(3,3));
    end
    fclose(fid);
    
    
    
end




return;
%_______________________________________________________________________

% Scan data base of DICOM file 
% See dicom_read_singlefile() for details. 
%
%       infoStruct = dicom_scan_singlefile(fileName,byteOrderFlag)
%       Does not read the data but extracts only essential information
%       from the DICOM file. Returns also the position parameter necessary
%       for cardiac data evaluation.
%       Example:
%       infoStruct = dicom_scan_singlefile('c:\tmp\myDicom.acr');
%
%               Harald Fischer
%               1/01
%
%     PC
%
% Comments and Changes
%     - 17.09.01 Harald Fischer: now reads also the series number
%     - 17.09.01 Harald Fischer: skipps elements of length 0
%     - 02.11.01 Harald Fischer: reads also ac-date and Patient-ID
%     - 12.12.01 Harald Fischer: some warnings now send to gui handle, not command line


function infoStruct = dicom_scan_singlefile(fileName,byteOrderFlag,messageHandle)


%%% set default values, error check
infoStruct    = [];
success       = 1;
byteorder     ='ieee-le'; % 'ieee-be'
explicit      = 0;
msgHandle     = []; 

if nargin>3 | nargin<1
    warning('Wrong number of input arguments');
    return;
end
if nargin>=2
    if byteOrderFlag>0
        byteorder = 'ieee-be';
    end
end
if nargin==3
    msgHandle = messageHandle;
end
%%% End of: set default values, error check



%%% open file
if nargin==0 | isempty(exist(fileName)) | exist(fileName)~=2
    [fileStr,dirStr]=uigetfile('*.*','Load a DICOM file');
  %  fileName = full_filename(dirStr,fileStr);
  fileName = fullfile(dirStr,fileStr);
end
[fileID,message]=fopen(fileName,'r',byteorder);
if fileID<0,
    disp('file could not be opened');
    return;
end;
%%% End of: open file



%%% get length of file
fseek(fileID,0,'eof');
fileLength = ftell(fileID);
fseek(fileID,0,'bof');
%%% End of: get length of file



%%% detect EXPLICIT dicom 
epilog = fread(fileID,132,'uchar');
dicmStr = sprintf('%s',char(epilog(129:132)'));
if strcmp(dicmStr,'DICM')
    explicit = 1;
else
    fseek(fileID,0,'bof');
end
%%% End of: detect EXPLICIT dicom



%%% read file
readStruct = local_read_dicom_file(fileID, fileLength, explicit);
fclose(fileID);
%%% End of: read file


%%% use other byte order if necessary
if (readStruct.errorFlag)
    
    messageStr  = sprintf('Warning: byteorder %s probably not correct, try to use other byteorder',byteorder);
    set(msgHandle,'String',messageStr);
    
    if byteorder=='ieee-le'
        byteorder='ieee-be';
    else
        byteorder='ieee-le';
    end
    
    [fileID,message]=fopen(fileName,'r',byteorder);
    if fileID<0,
        disp('file could not be opened');
        return;
    end;
    readStruct = local_read_dicom_file(fileID, fileLength, explicit);
    fclose(fileID);
    
end
%%% End of: use other byte order if necessary


%%% error check
if (readStruct.errorFlag)
    messageStr  = sprintf('Warning: cannot load file: %s',fileName);
    set(msgHandle,'String',messageStr);
    return;
end   
%%% End  of: error check


%%% return
infoStruct.size_x     = readStruct.size_x;
infoStruct.size_y     = readStruct.size_y;
infoStruct.patientStr = readStruct.patientStr;
infoStruct.idStr      = readStruct.idStr;
infoStruct.voxVc      = readStruct.voxVc;
infoStruct.tr         = readStruct.tr;
infoStruct.teVc       = readStruct.teVc;
infoStruct.tiVc       = readStruct.tiVc;
infoStruct.pos        = readStruct.position;
infoStruct.image      = readStruct.image;
infoStruct.series     = readStruct.series;
infoStruct.EchoNumbers     = readStruct.EchoNumbers;
infoStruct.date       = readStruct.dateStr;
%%% End of: return


%%%%%%%%%%%%%%%%%%%%%%%%%% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%

function readStruct = local_read_dicom_file(fileID, fileLength, explicit)

success = 1;
changeTransferSyntax  = 0;
transferSyntax        = '';
readStruct.size_x     = -1;
readStruct.size_y     = -1;
readStruct.patientStr = '';
readStruct.idStr      = '';
readStruct.voxVc      = [];
readStruct.tr         = [];
readStruct.teVc       = [];
readStruct.tiVc       = [];
readStruct.errorFlag  = 0;
readStruct.position   = [];
readStruct.image      = -1;
readStruct.series     = -1;
readStruct.EchoNumbers = -1;
readStruct.dateStr    = '';


while success>0
    
    
    %%% reading group and element of cluster
    if (ftell(fileID)==fileLength | ftell(fileID)+8==fileLength)
        break;
    elseif (ftell(fileID)+8>fileLength)
        readStruct.errorFlag=1; break;
    end
    
    [group,success]   = fread(fileID,2,'uchar');
    [element,success] = fread(fileID,2,'uchar');
    group   = group(2)*256+group(1);
    element = element(2)*256+element(1);
    %%% End of: reading group and element of cluster
    
    
    
    
    %%% eventually change transfer syntax
    if group>2 & explicit & changeTransferSyntax
        explicit = 0;   
    end
    %%% End of: eventually change transfer syntax
    
    
    %%% reading length of cluster
    if explicit
        [typVc ,success] = fread(fileID,2,'uchar');
        typStr = sprintf('%s',char(typVc'));
        if strcmp(typStr,'OB') | strcmp(typStr,'OW') | strcmp(typStr,'SQ') | strcmp(typStr,'UN')
            [dummy,success]  = fread(fileID,2,'uchar');
            [len,success]    = fread(fileID,4,'uchar');
        else
            [len,success]    = fread(fileID,2,'uchar');
        end
    else
        [len,success]    = fread(fileID,4,'uchar');
    end
    
    if length(len(:))==2
        length2  = len(1)+len(2)*256;
    else
        length2  = len(1)+len(2)*256+len(3)*256*256+len(4)*256*256*256;
    end
    %%% End of: reading length of cluster
    
    
    %%% error check
    if group<0 | element<0 | length2<0 | length2>fileLength
        readStruct.errorFlag=1; break;
    end
    %%% End of: error check
    
    
    
    %%% skipp group of no interesst 
    if element==0 & group~=2 & group~=16 & group~=24 & group~=32 &group~=40 & group~=32736
        [forward,success]=fread(fileID,1,'int');
        fseek(fileID,forward,'cof');
        %%% End of: skipp a group of no interesst
        
        
        %%% searching and reading clusters of interest
    elseif group==2 & element==16 & length2>0,
        [something,success]=fread(fileID,length2,'uchar');
        transferSyntax = sprintf('%s',char(something)');
        if strcmp(transferSyntax,'1.2.840.10008.1.2')
            changeTransferSyntax=1;
        end
        
    elseif group==8 & element==34 & length2>0,
        [readStruct.dateStr,success]=fread(fileID,length2,'uchar');
        readStruct.dateStr = sprintf('%s',char(readStruct.dateStr)');
        readStruct.dateStr = deblank(readStruct.dateStr);
        
    elseif group==16 & element==16 & length2>0,
        [readStruct.patientStr,success]=fread(fileID,length2,'uchar');
        readStruct.patientStr = sprintf('%s',char(readStruct.patientStr)');
        readStruct.patientStr = deblank(readStruct.patientStr);
        
    elseif group==16 & element==32 & length2>0,
        [readStruct.idStr,success]=fread(fileID,length2,'uchar');
        readStruct.idStr = sprintf('%s',char(readStruct.idStr)');
        readStruct.idStr = deblank(readStruct.idStr);
        
    elseif group==24 & element==80 & length2>0,
        [sliceThickness,success]=fread(fileID,length2,'uchar');
        sliceThickness      = sprintf('%s',char(sliceThickness)');
        sliceThickness      = str2num(sliceThickness);
        readStruct.voxVc(3) = sliceThickness;
        
    elseif group==24 & element==128 & length2>0,
        [readStruct.tr,success]=fread(fileID,length2,'uchar');
        readStruct.tr = sprintf('%s',char(readStruct.tr)');
        readStruct.tr = str2num(readStruct.tr);
        
    elseif group==24 & element==129 & length2>0,
        [readStruct.teVc,success]=fread(fileID,length2,'uchar');
        readStruct.teVc = sprintf('%s',char(readStruct.teVc)');
        readStruct.teVc = str2num(readStruct.teVc);
        
    elseif group==24 & element==130 & length2>0,
        [readStruct.tiVc,success]=fread(fileID,length2,'uchar');
        readStruct.tiVc = sprintf('%s',char(readStruct.tiVc)');
        readStruct.tiVc = str2num(readStruct.tiVc);     
        
    elseif group==24 & element==134 & length2>0,
        [readStruct.EchoNumbers,success]=fread(fileID,length2,'uchar');
        readStruct.EchoNumbers = sprintf('%s',char(readStruct.EchoNumbers)');
        readStruct.EchoNumbers = str2num(readStruct.EchoNumbers);
        
    elseif group==24 & element==136 & length2>0,
        [gap,success]=fread(fileID,length2,'uchar');
        gap          = sprintf('%s',char(gap)');
        gap          = str2num(gap);
        readStruct.voxVc(4) = gap-readStruct.voxVc(3);
           
    elseif group==32 & element==4161 & length2>0,
        [readStruct.position,success]=fread(fileID,length2,'uchar');
        readStruct.position = sprintf('%s',char(readStruct.position)');
        readStruct.position = str2num(readStruct.position);
        
    elseif group==32 & element==17 & length2>0,
        [readStruct.series,success]=fread(fileID,length2,'uchar');
        readStruct.series = sprintf('%s',char(readStruct.series)');
        readStruct.series = str2num(readStruct.series);
        
   
        
    elseif group==32 & element==19 & length2>0,
        [readStruct.image,success]=fread(fileID,length2,'uchar');
        readStruct.image= sprintf('%s',char(readStruct.image)');
        readStruct.image= str2num(readStruct.image);
        
    elseif group==40 & element==48 & length2>0,
        [voxStr,success]=fread(fileID,length2,'uchar');
        voxStr = sprintf('%s',char(voxStr)');
        voxStr = strrep(voxStr,'\',' ');
        [vox_xStr, vox_yStr] = strtok(voxStr);
        readStruct.voxVc(1)  = str2num(vox_xStr);
        readStruct.voxVc(2)  = str2num(vox_yStr);
        
    elseif group==40 & element==16 & length2>0,
        [readStruct.size_y,success]=fread(fileID,1,'ushort');
        
    elseif group==40 & element==17 & length2>0,
        [readStruct.size_x,success]=fread(fileID,1,'ushort');
        
        %elseif group==40 & element==256,
        %%%   [dataType,success]=fread(fileID,length2,'uchar');
        
    elseif group==32736 & element==16 & length2>0,
        if (ftell(fileID)+readStruct.size_x*2*readStruct.size_y > fileLength)
            readStruct.errorFlag=1; break;
        end
        return;
        
    else
        fseek(fileID,length2,'cof');
    end
    %%% searching and reading clusters of interrest
    
end %endwhile

function magn = magnitude (vector)
total = 0.0;
for i = 1 : length(vector)
    total = total + vector(i)*vector(i);
end
magn = sqrt(total);
return

function parsed_name = legitimise_filename (input_filename)
% Check the Patient ID for dodgy characters and replace them, also change the PatientID in the Patient structure as the same time
parsed_name = input_filename;        
colons = strfind(parsed_name, ':');
if (isempty(colons))
    % Patient_ID does not contain any ':'
else
    parsed_name(colons) = '-';
    disp(['Colons found while running legitimise_filename, therefore these have been swapped for ;']);
end
spaces = strfind(parsed_name, ' ');
if (isempty(spaces))
    % Patient_ID does not contain any spaces
else
    parsed_name(spaces) = '_';
    disp(['Spaces found while running legitimise_filename, therefore these have been swapped for _']);
end

return

function [value, success] = values_in_strings(input_string, variable_to_match)

%
% Take in a string containing a load of variables i.e.
% sGroupArray.anMember[37]                 = 37
% sGroupArray.anMember[38]                 = 38
% sGroupArray.anMember[39]                 = 39
% sGroupArray.anMember[40]                 = 40
% sGroupArray.anMember[41]                 = 41
%
% and a variable_to_match
% e.g. sGroupArray.anMember[39]
% 
% and return the values
% e.g. 39
%
% If a single version of this variable is found with no problems then
% success == 1
% Multiple occurences success == number of occurences, value will give the
% result of the first occurence
% No occurences or failure from other cause
% Success == 0
%

success = 0;
value = '';

indices  = strfind(input_string, variable_to_match);

if (isempty(indices))
    success = 0;
else
    success = length(indices);
    chunk = input_string(indices(1):indices(1)+1000);
    next_equals = strfind(chunk,'=');
    if (isempty(next_equals))
        success = 0;
    else
        chunk = chunk(next_equals(1)+1:end);
        next_return = strfind(chunk,char(10));
        if (isempty(next_return))
            success = 0;
        else
            value = chunk(1:next_return(1)-1);
        end
    end
end

return 

function [error] = Make_bvecs103(output_filename, info_tmp, mdr_lBerhens_flag, diffusion)

version = 1.06;

% disp(['version = ' num2str(version)])

% Read_mdr
%
% Read the orientation file
% Tim, Matt (and Clare)
% 5/06/03
%
%
% Updated to allow it to be incorporated into the Dicom2Ana code
% This includes input that provides the slice orientation information (via
% info_tmp), 
% the number of diffusion encoding directions, and 
% information on the filenames for the output
%structure describing the diffusion weightings (diffusion structure)

%aux_output_pathname = 'M:\CLARE\diffusion\Patient_N01369\Study_2\';
%aux_file = 'images_10_ep2d_diff.mdr';
global error vector_b


if (nargin < 4)
    
    %Read the input information from the user
    [aux_file,aux_output_pathname] = uigetfile( {'*.mdr', 'MDR orientation files (*.mdr)'},'Select an mdr file');
    mach = 'native';
    
    fid             =  fopen([aux_output_pathname aux_file ],'r',mach);
    
    pixel_origin = fread(fid,  3,   'float')       % 3 floats
    %Next the vector describing how to get from one slice to the next, this
    %gives weird values for non-multi slice data sets, and for multi-angle
    %obliquing!
    slice_vector = fread(fid,  3,   'float')       % 3 floats
    %Next the in-plane vector for the data as it is in the Dicom file
    first_vector = fread(fid,   3,   'float')       % 3 floats for read direction
    second_vector = fread(fid,  3,   'float')       % 3 floats for phase encode direction
    %Next include the patient orientation information i.e.
    % HFS (head first suppine), HFP (head first prone), FFS (feet first supine) etc.
    Patient_Position = char(fread(fid, 3,'char'))'        %3 char
    %Finally the information concerning whether the images have been flipped in the output file
    flip_flag = fread(fid,     1,   'int16')       % 1 int (as described above)
    fclose(fid);
    
    directions = input('HOW MANY DIRECTIONS (12,60, 124, 252 or 512)?')
    
else
    pixel_origin = info_tmp.ImagePositionPatient;       % 3 floats
    %Next the vector describing how to get from one slice to the next
    slice_vector = [0.0 , 0.0 , 0.0 ]; %this will be created when flag == 1 (when we get to the next slice)
    %Next the in-plane vector for the data as it is in the Dicom file
    first_vector = info_tmp.ImageOrientationPatient(1:3);% 3 floats for read direction
    second_vector  = info_tmp.ImageOrientationPatient(4:6);       % 3 floats for phase encode direction
    %This describes the orientation of the patient within the magnet
    slice_vector = cross(first_vector, second_vector);
    Patient_Position = info_tmp.PatientPosition;
    
end

if (strcmp(Patient_Position,'HFS') == 1)
    %This is fine
else
    disp('Patient put into magnet upside down!!!')
end
if (sum(slice_vector) == 0)
    %y = sin((128/pi)*([0:3:8192 , 0:4096])).^4;
    %sound(y,12345)
    disp('Slice vector is zero, this may be single slice, fixing');
    slice_vector = cross (first_vector, second_vector);
end

if ( abs(dot(slice_vector, first_vector)) <= 0.00001 & abs(dot(slice_vector,second_vector)) <= 0.000001 ) 
    % Leave the slice_vector as it is
else
    %y = sin((128/pi)*([0:3:8192 , 0:4096])).^4;
    %sound(y,12345)
    disp('Slice vector is not perpendicular to read and phase! this may be MAO, fixing');
    slice_vector = cross (first_vector, second_vector);
end

if (nargin < 4)
    % This is run real-time by the user, these special cases have been
    % programmed in
    
    switch diffusion.directions
        case 12
            b_val = 1000;
            Num_b_zeros = 1;
            vector_b = [];
            diffusion_code(12)        
            zero_vectors = zeros(3,Num_b_zeros); % b = 0 scans
            bvecs_old = [ zero_vectors vector_b' ];
            
        case 60
            b_val = 1000;
            Num_b_zeros = 1;
            vector_b = [];
            diffusion_code(29)
            first_29 = vector_b;
            vector_b = [];
            diffusion_code(31);
            second_31 = vector_b;
            
            zero_vectors = zeros(3,Num_b_zeros); % b = 0 scans
            bvecs_old = [ zero_vectors first_29' second_31' ];
            
        case 124
            b_val = 10000;
            Num_b_zeros = 2;
            zero_vectors = zeros(3,Num_b_zeros); % b = 0 scans
            bvecs_old = [zero_vectors];
            for mdr_lBerhens_flag = 1 : 2
                vector_b = [];
                diffusion_code(62,mdr_lBerhens_flag);
                bvecs_old = [bvecs_old vector_b'];
            end
        case 252
            b_val = 10000;
            Num_b_zeros = 4;
            zero_vectors = zeros(3,Num_b_zeros); % b = 0 scans
            bvecs_old = [zero_vectors];
            for mdr_lBerhens_flag = 1 : 4
                vector_b = [];
                diffusion_code(63,mdr_lBerhens_flag);
                bvecs_old = [bvecs_old vector_b'];
            end
        case 512
            b_val = 10000;
            Num_b_zeros = 9;
            zero_vectors = zeros(3,Num_b_zeros); % b = 0 scans
            bvecs_old = [zero_vectors];
            for mdr_lBerhens_flag = 1 : 8
                vector_b = [];
                diffusion_code(64,mdr_lBerhens_flag);
                bvecs_old = [bvecs_old vector_b'];
            end
        otherwise
            disp('Only 12, 60, 124, 252 or 512 directions are catered for in the code at present');
    end
    
    
else
    
    % This is baing called as a function (probably from Dicom2Ana)
    % Presently restrict this to 2 diffusion weighting values (0 and
    % something)
    Num_b_zeros = 1;
    zero_vectors = zeros(3, Num_b_zeros);
    
    bvecs_old = [zero_vectors];
    
    vector_b = [];
    diffusion_code(diffusion.directions , mdr_lBerhens_flag);
    if (diffusion.weightings > 1)
        
        for weightings_counter = 2 : diffusion.weightings
            bvecs_old = [bvecs_old vector_b']; % Keep appending the information for each direction
        end
    end
    
end
% The analyse format uses foot to head as the positive coordinate (for z)
% The Siemens uses foot to heaed as the NEGATIVE coordinate (for z)
% SO, we need this very ugly flipping.  Be afraid!
bvecs_old(3,:) = -bvecs_old(3,:);

bvecs_new=correct_bvecs(bvecs_old,slice_vector,first_vector,second_vector);
% plot3(bvecs_new(1,:)',bvecs_new(2,:)',bvecs_new(3,:)','r-');

if (nargin < 4)
    fid2 = fopen([aux_output_pathname '\bvecs'],'w');
else
    fid2 = fopen([output_filename '.bvecs'],'w');
end



for counter2 = 1:3
    for counter = 1: size(bvecs_new,2)
        fprintf(fid2,'%d ',bvecs_new(counter2,counter));
    end
    fprintf(fid2,'\n');
end



fclose(fid2);

% write the b-vals file too

if(nargin < 4)
    fid2 = fopen([aux_output_pathname '\bvals'],'w');
else
    fid2 = fopen([output_filename '.bvals'],'w'); 
end

if (nargin < 4)
    for counter = 1: size(bvecs_new,2)
        if (counter <= Num_b_zeros)
            fprintf(fid2,'%d ',0);
        else
            fprintf(fid2,'%d ',b_val);
        end
    end
else
    
    for counter = 1 : size(bvecs_new,2)
        fprintf(fid2,'%d ',diffusion.b_val(counter));
    end
    
    fprintf(fid2,'\n');
end

fclose(fid2);





function bvecs_new=correct_bvecs(bvecs_old,slice_dir,read_dir,phase_dir);

%CORRECT_BVECS - corrects bvecs when slices are acquired off-axis
%
% bvecs_new=correct_bvecs(bvecs_old,slice_dir,read_dir,phase_dir);
%
% bvecs_old is 3xN
slice_dir = slice_dir ./ (sqrt(sum(slice_dir.^2)));
read_dir = read_dir ./ (sqrt(sum(read_dir.^2)));
phase_dir = phase_dir ./ (sqrt(sum(phase_dir.^2)));

n=size(bvecs_old,2);
sd=repmat(slice_dir,1,n);
rd=repmat(read_dir,1,n);
pd=repmat(phase_dir,1,n);

bvecs_new=[sum(bvecs_old.*rd);sum(bvecs_old.*pd);sum(bvecs_old.*sd)];

function  setVector(index, x,y,z,hmm)
global error vector_b
if (hmm == 0)
    error = 1;
end
vector_b(index+1,1) = x;
vector_b(index+1,2) = y;
vector_b(index+1,3) = z;

function nothing = diffusion_code(m_Directions, mdr_lBerhens_flag)

% This is the code in the Siemens code.
% If the Siemens code is changed THIS HAS TO BE CHANGED
%
XYZ = 1;

switch (m_Directions) 
    case 6
        setVector (0,  1.0, 0.0,  1.0, XYZ) ;
        setVector (1, -1.0, 0.0,  1.0, XYZ) ;
        setVector (2, 0.0,  1.0,  1.0, XYZ) ;
        setVector (3, 0.0,  1.0, -1.0, XYZ) ;
        setVector (4,  1.0,  1.0, 0.0, XYZ) ;
        setVector (5, -1.0,  1.0, 0.0, XYZ) ;
        
    case 12
        setVector (0, 1.0, 0., +0.5, XYZ) ;
        setVector (1, 0., +0.5, 1.0, XYZ) ;
        setVector (2, +0.5, 1.0, 0., XYZ) ;
        setVector (3, 1.0, +0.5, 0., XYZ) ;
        setVector (4, 0., 1.0, +0.5, XYZ) ;
        setVector (5, +0.5, 0., 1.0, XYZ) ;
        setVector (6, 1.0, 0., -0.5, XYZ) ;
        setVector (7, 0., -0.5, 1.0, XYZ) ;
        setVector (8, -0.5, 1.0, 0., XYZ) ;
        setVector (9, 1.0, -0.5, 0., XYZ) ;
        setVector (10, 0., 1.0, -0.5, XYZ) ;
        setVector (11, -0.5, 0., 1.0, XYZ) ;
        
    case 7  
        setVector (0 , 0.99998142 , 0.43098525 , 0.30088211, XYZ) ; 
        setVector (1 , -0.30129064 , 0.99999636 , -0.43066507, XYZ) ; 
        setVector (2 , -0.15482230 , -0.99991262 , -0.50244002, XYZ) ; 
        setVector (3 , -1.00000000 , 0.50227889 , 0.15478077, XYZ) ; 
        setVector (4 , -0.50248532 , 0.15470045 , -0.99990872, XYZ) ; 
        setVector (5 , -0.65194289 , -0.65232102 , 0.65244841, XYZ) ; 
        setVector (6 , 0.43091770 , -0.30117948 , -0.99992101, XYZ) ; 
        
    case 8  
        setVector (0 , 0.92099681 , -0.45502853 , 0.45455671, XYZ) ; 
        setVector (1 , 0.98266142 , 0.46242991 , 0.28712888, XYZ) ; 
        setVector (2 , -0.47387523 , -0.99992077 , 0.19367141, XYZ) ; 
        setVector (3 , -0.19919359 , 0.78176980 , -0.78170697, XYZ) ; 
        setVector (4 , 0.47374712 , 0.19359248 , -0.99999675, XYZ) ; 
        setVector (5 , -0.30207865 , -0.41311646 , -0.99999558, XYZ) ; 
        setVector (6 , -0.30161411 , 1.00000000 , 0.41344506, XYZ) ; 
        setVector (7 , 0.98243070 , -0.28740150 , -0.46275067, XYZ) ; 
        
    case 9  
        setVector (0 , 0.99998927 , 0.53837331 , -0.20408203, XYZ) ; 
        setVector (1 , 0.99743701 , -0.41514757 , -0.40527243, XYZ) ; 
        setVector (2 , -0.08170540 , 0.58366806 , 0.99203309, XYZ) ; 
        setVector (3 , -0.99456798 , -0.07770708 , -0.57988789, XYZ) ; 
        setVector (4 , -0.58197340 , 0.99298976 , -0.08217132, XYZ) ; 
        setVector (5 , -0.20283431 , -1.00000000 , 0.53882471, XYZ) ; 
        setVector (6 , -0.40617235 , -0.99747775 , -0.41416905, XYZ) ; 
        setVector (7 , -0.41327285 , 0.40087724 , -0.99998841, XYZ) ; 
        setVector (8 , 0.53855100 , 0.20364413 , -0.99998286, XYZ) ; 
        
    case 10  
        setVector (0 , 0.43165406 , -0.28687434 , -0.99999971, XYZ) ; 
        setVector (1 , -0.44714056 , -0.26208323 , -0.99999962, XYZ) ; 
        setVector (2 , 0.85973707 , -0.65145758 , -0.32415540, XYZ) ; 
        setVector (3 , -0.01439669 , -0.51808713 , 1.00000000, XYZ) ; 
        setVector (4 , 0.98340042 , 0.22931730 , -0.49895763, XYZ) ; 
        setVector (5 , 0.96717035 , -0.28818444 , 0.50015275, XYZ) ; 
        setVector (6 , -0.40984299 , 0.98852664 , -0.35137636, XYZ) ; 
        setVector (7 , -0.02979540 , -0.99961216 , -0.51817884, XYZ) ; 
        setVector (8 , -0.89688332 , -0.59912778 , -0.32444992, XYZ) ; 
        setVector (9 , -0.46903176 , -0.96125005 , 0.35302845, XYZ) ; 
        
    case 11  
        setVector (0 , -0.99889874 , -0.38175995 , 0.31085451, XYZ) ; 
        setVector (1 , -0.33891423 , 0.84599477 , -0.63999990, XYZ) ; 
        setVector (2 , -0.84863867 , 0.36839447 , 0.61989338, XYZ) ; 
        setVector (3 , -0.99999919 , 0.42621226 , -0.24189793, XYZ) ; 
        setVector (4 , 0.44823584 , -0.99998999 , -0.19817801, XYZ) ; 
        setVector (5 , -0.41439564 , -0.98417303 , 0.31598984, XYZ) ; 
        setVector (6 , -0.91463376 , -0.38227902 , -0.50742263, XYZ) ; 
        setVector (7 , -0.39072656 , -0.37400055 , 0.97346095, XYZ) ; 
        setVector (8 , -0.12049865 , 0.47502627 , 1.00000000, XYZ) ; 
        setVector (9 , -0.47673847 , 0.11366899 , -0.99998484, XYZ) ; 
        setVector (10 , 0.33182972 , 0.93838319 , 0.49949565, XYZ) ; 
        
    case 13  
        setVector (0 , -0.79299033 , -0.69134262 , -0.37538795, XYZ) ; 
        setVector (1 , 0.98360719 , 0.38308820 , -0.36532824, XYZ) ; 
        setVector (2 , 0.43852157 , 0.40626874 , -0.94358303, XYZ) ; 
        setVector (3 , -0.11387966 , -0.97589016 , -0.53138895, XYZ) ; 
        setVector (4 , 0.99999604 , -0.33591398 , -0.36725211, XYZ) ; 
        setVector (5 , -0.99989752 , -0.02496075 , -0.49727891, XYZ) ; 
        setVector (6 , -0.32819412 , -0.37415640 , -1.00000000, XYZ) ; 
        setVector (7 , 0.45614001 , -0.38442336 , -0.94438306, XYZ) ; 
        setVector (8 , -0.82677691 , 0.65032791 , -0.37578976, XYZ) ; 
        setVector (9 , 0.48112610 , 0.97437396 , -0.25849075, XYZ) ; 
        setVector (10 , -0.52778774 , 0.95023387 , 0.25729376, XYZ) ; 
        setVector (11 , 0.34408874 , -0.35959377 , 0.99999982, XYZ) ; 
        setVector (12 , -0.16233358 , 0.97009637 , -0.52940082, XYZ) ; 
        
    case 14  
        setVector (0 , -0.39430490 , 0.73327016 , -0.70463335, XYZ) ; 
        setVector (1 , -0.42576291 , 0.09191992 , -0.99997315, XYZ) ; 
        setVector (2 , -0.43085715 , 0.99997372 , -0.06390866, XYZ) ; 
        setVector (3 , -0.32273871 , -0.59439366 , -0.85569009, XYZ) ; 
        setVector (4 , -0.91625516 , -0.29161748 , -0.51488383, XYZ) ; 
        setVector (5 , -0.30932815 , -0.95108098 , 0.43523637, XYZ) ; 
        setVector (6 , -0.93446325 , 0.45529698 , 0.33038269, XYZ) ; 
        setVector (7 , 0.27584399 , 0.46152540 , -0.94898579, XYZ) ; 
        setVector (8 , -0.35121327 , 0.83496175 , 0.60758354, XYZ) ; 
        setVector (9 , 0.95837396 , -0.38189611 , 0.35404004, XYZ) ; 
        setVector (10 , 0.90271407 , 0.60568914 , -0.08898103, XYZ) ; 
        setVector (11 , -0.36197972 , 0.24217477 , 0.99999585, XYZ) ; 
        setVector (12 , 0.34133660 , 1.00000000 , 0.27047913, XYZ) ; 
        setVector (13 , 0.88601485 , 0.17536065 , -0.61147033, XYZ) ; 
        
    case 15  
        setVector (0 , -0.64443738 , 0.56844502 , -0.64853185, XYZ) ; 
        setVector (1 , -0.61625596 , -0.63485840 , -0.61335654, XYZ) ; 
        setVector (2 , -0.44801040 , -0.02491320 , -0.97861579, XYZ) ; 
        setVector (3 , -0.99999632 , 0.04226632 , 0.39653973, XYZ) ; 
        setVector (4 , 0.16185960 , -0.99996617 , 0.36454352, XYZ) ; 
        setVector (5 , 0.00183795 , 0.55183640 , 0.92438956, XYZ) ; 
        setVector (6 , -0.00146677 , -0.55304613 , 0.92366697, XYZ) ; 
        setVector (7 , -1.00000000 , 0.04119509 , -0.39664317, XYZ) ; 
        setVector (8 , 0.61731977 , 0.63225065 , -0.61497825, XYZ) ; 
        setVector (9 , -0.42876452 , -0.98750979 , 0.00287118, XYZ) ; 
        setVector (10 , 0.64373913 , -0.56966577 , -0.64815406, XYZ) ; 
        setVector (11 , 0.44715627 , 0.02430073 , -0.97902175, XYZ) ; 
        setVector (12 , 0.94647341 , 0.51304086 , -0.00002688, XYZ) ; 
        setVector (13 , -0.81454044 , 0.70395078 , 0.00005999, XYZ) ; 
        setVector (14 , -0.15830811 , 0.99997965 , 0.36606296, XYZ) ; 
        
    case 16  
        setVector (0 , -0.00001553 , 0.38201430 , 0.99999988, XYZ) ; 
        setVector (1 , 0.61798434 , -0.61791468 , -0.61823256, XYZ) ; 
        setVector (2 , 0.61781860 , 0.61804283 , 0.61827011, XYZ) ; 
        setVector (3 , 0.56273600 , 0.00022431 , -0.91063870, XYZ) ; 
        setVector (4 , -0.00009234 , 0.91073155 , 0.56258575, XYZ) ; 
        setVector (5 , 0.91065055 , -0.56271672 , -0.00039701, XYZ) ; 
        setVector (6 , 0.38201954 , -0.99999788 , 0.00012249, XYZ) ; 
        setVector (7 , 0.99997300 , 0.00028532 , -0.38208455, XYZ) ; 
        setVector (8 , 0.56275987 , -0.00009859 , 0.91062397, XYZ) ; 
        setVector (9 , 0.61821521 , 0.61816372 , -0.61775259, XYZ) ; 
        setVector (10 , -0.38201398 , -1.00000000 , 0.00012463, XYZ) ; 
        setVector (11 , -0.00004025 , 0.38202584 , -0.99999547, XYZ) ; 
        setVector (12 , 0.61818899 , -0.61809020 , 0.61785239, XYZ) ; 
        setVector (13 , 0.00005357 , -0.91052656 , 0.56291746, XYZ) ; 
        setVector (14 , 0.99999867 , -0.00023739 , 0.38201740, XYZ) ; 
        setVector (15 , 0.91068264 , 0.56266457 , 0.00063772, XYZ) ; 
        
    case 17  
        setVector (0 , 0.89373461 , -0.05747407 , -0.56589914, XYZ) ; 
        setVector (1 , -0.58617154 , 0.60969194 , 0.63795402, XYZ) ; 
        setVector (2 , -0.99999840 , 0.34954119 , 0.01143814, XYZ) ; 
        setVector (3 , -0.13878729 , 0.99999857 , -0.32101030, XYZ) ; 
        setVector (4 , -0.62146151 , -0.76154930 , -0.39513912, XYZ) ; 
        setVector (5 , 0.06637616 , 0.66579182 , -0.82135384, XYZ) ; 
        setVector (6 , 0.55136996 , 0.34431610 , 0.83650718, XYZ) ; 
        setVector (7 , 0.01301419 , -0.57362614 , -0.89055620, XYZ) ; 
        setVector (8 , 0.34972091 , -0.00139951 , -1.00000000, XYZ) ; 
        setVector (9 , -0.93224533 , -0.00186983 , -0.50321151, XYZ) ; 
        setVector (10 , -0.98198105 , -0.39549871 , 0.04000829, XYZ) ; 
        setVector (11 , 0.64118148 , -0.59765061 , 0.59498465, XYZ) ; 
        setVector (12 , 0.47975574 , 0.92447044 , -0.19363756, XYZ) ; 
        setVector (13 , -0.65960291 , -0.51722316 , 0.64785097, XYZ) ; 
        setVector (14 , 0.00930960 , -0.99433994 , -0.36538758, XYZ) ; 
        setVector (15 , 0.29595558 , -0.18632837 , 0.99999935, XYZ) ; 
        setVector (16 , -0.62008225 , 0.85732593 , 0.05288603, XYZ) ; 
        
    case 18  
        setVector (0 , -0.81295112 , 0.62227054 , 0.29180983, XYZ) ; 
        setVector (1 , 0.31535921 , 0.18388779 , -0.99999849, XYZ) ; 
        setVector (2 , 0.06284082 , 0.42146740 , 0.97554056, XYZ) ; 
        setVector (3 , -0.57263159 , -0.62036467 , -0.64846270, XYZ) ; 
        setVector (4 , 0.29077077 , -0.81363483 , -0.62186321, XYZ) ; 
        setVector (5 , -0.33051688 , 0.15498941 , -1.00000000, XYZ) ; 
        setVector (6 , 0.03183174 , -0.69439955 , 0.80626245, XYZ) ; 
        setVector (7 , 0.64796519 , 0.57164812 , -0.62178989, XYZ) ; 
        setVector (8 , -0.97561281 , -0.06287621 , 0.42129485, XYZ) ; 
        setVector (9 , 0.18560876 , 0.99954494 , 0.31578856, XYZ) ; 
        setVector (10 , -0.42062673 , 0.97595542 , -0.06202665, XYZ) ; 
        setVector (11 , 0.69488425 , 0.80582133 , -0.03241905, XYZ) ; 
        setVector (12 , -0.80647820 , -0.03030406 , -0.69421733, XYZ) ; 
        setVector (13 , 0.62190088 , -0.29168193 , -0.81327983, XYZ) ; 
        setVector (14 , 0.15520852 , 0.99999896 , -0.33041719, XYZ) ; 
        setVector (15 , -0.99988517 , 0.31593801 , -0.18351007, XYZ) ; 
        setVector (16 , -0.62112283 , 0.64888197 , -0.57133331, XYZ) ; 
        setVector (17 , 0.99996725 , 0.33030462 , 0.15565177, XYZ) ; 
        
    case 19  
        setVector (0 , -0.33295126 , -0.33236353 , -0.95044126, XYZ) ; 
        setVector (1 , 0.36030716 , -0.26054477 , 0.96278545, XYZ) ; 
        setVector (2 , 0.67483393 , 0.67692185 , -0.45938724, XYZ) ; 
        setVector (3 , -0.75063284 , 0.74913943 , 0.00105219, XYZ) ; 
        setVector (4 , -0.58720617 , 0.64047154 , -0.60798501, XYZ) ; 
        setVector (5 , 0.64079708 , -0.58511980 , -0.60965119, XYZ) ; 
        setVector (6 , -0.11582429 , -0.76219765 , 0.72821709, XYZ) ; 
        setVector (7 , -0.99920697 , -0.35333520 , 0.03742097, XYZ) ; 
        setVector (8 , -0.24881386 , -0.25050623 , 0.99999947, XYZ) ; 
        setVector (9 , 0.20523830 , -0.99996255 , 0.28742442, XYZ) ; 
        setVector (10 , -0.17879048 , -0.81690884 , -0.65219213, XYZ) ; 
        setVector (11 , -1.00000000 , 0.20285743 , 0.28898013, XYZ) ; 
        setVector (12 , 0.35120234 , 0.99999083 , -0.03654996, XYZ) ; 
        setVector (13 , 0.70724612 , 0.70827261 , 0.35044757, XYZ) ; 
        setVector (14 , 0.98034506 , -0.24905749 , 0.31867629, XYZ) ; 
        setVector (15 , -0.25064179 , 0.97938661 , 0.32037668, XYZ) ; 
        setVector (16 , 0.76108433 , 0.11844220 , -0.72896004, XYZ) ; 
        setVector (17 , 0.26001294 , -0.35850303 , -0.96360235, XYZ) ; 
        setVector (18 , -0.81743245 , -0.17964934 , -0.65129944, XYZ) ; 
        
    case 20  
        setVector (0 , -0.42237219 , 0.87193592 , -0.42236450, XYZ) ; 
        setVector (1 , -0.82387533 , 0.65659671 , -0.08469060, XYZ) ; 
        setVector (2 , -0.55973177 , -0.61935036 , 0.64820349, XYZ) ; 
        setVector (3 , 0.33761585 , 0.94241131 , 0.33902618, XYZ) ; 
        setVector (4 , 0.99637052 , -0.05208659 , 0.34870483, XYZ) ; 
        setVector (5 , -0.43186221 , -0.48361558 , -0.83466957, XYZ) ; 
        setVector (6 , 0.15866862 , 1.00000000 , -0.30312796, XYZ) ; 
        setVector (7 , 0.69753177 , -0.38000305 , 0.69721547, XYZ) ; 
        setVector (8 , -0.65025696 , 0.61727867 , 0.55963847, XYZ) ; 
        setVector (9 , -0.08516271 , 0.65700848 , -0.82349831, XYZ) ; 
        setVector (10 , 0.99999540 , -0.23326202 , -0.25032045, XYZ) ; 
        setVector (11 , 0.30379675 , -0.99999542 , -0.15741340, XYZ) ; 
        setVector (12 , 0.35065279 , -0.05090425 , 0.99574779, XYZ) ; 
        setVector (13 , 0.68716360 , 0.79878669 , -0.08251246, XYZ) ; 
        setVector (14 , -0.21647163 , 0.33094552 , 0.98014151, XYZ) ; 
        setVector (15 , -0.74729761 , -0.00210110 , 0.74739826, XYZ) ; 
        setVector (16 , -0.24877246 , -0.23495166 , 0.99998614, XYZ) ; 
        setVector (17 , -0.98014750 , -0.33251949 , 0.21401856, XYZ) ; 
        setVector (18 , -0.08772471 , 0.79925257 , 0.68597521, XYZ) ; 
        setVector (19 , 0.83317147 , 0.47994535 , 0.43879386, XYZ) ; 
        
    case 21  
        setVector (0 , -0.21383014 , 0.83170127 , 0.61819362, XYZ) ; 
        setVector (1 , -0.61802407 , -0.21404862 , -0.83177108, XYZ) ; 
        setVector (2 , -0.61782826 , 0.34603206 , 0.78623390, XYZ) ; 
        setVector (3 , 0.21407743 , -0.27241485 , 0.99978732, XYZ) ; 
        setVector (4 , -1.00000000 , -0.00034631 , -0.34585194, XYZ) ; 
        setVector (5 , 0.78576457 , 0.27235669 , -0.65422439, XYZ) ; 
        setVector (6 , -0.99353062 , -0.34394676 , 0.11921085, XYZ) ; 
        setVector (7 , -0.34591441 , 0.99997834 , 0.00046970, XYZ) ; 
        setVector (8 , -0.83237879 , 0.61740821 , -0.21346272, XYZ) ; 
        setVector (9 , 0.00039737 , 0.34588445 , 0.99998874, XYZ) ; 
        setVector (10 , -0.70777241 , 0.34352364 , -0.70757573, XYZ) ; 
        setVector (11 , -0.65374586 , -0.78612900 , 0.27245409, XYZ) ; 
        setVector (12 , -0.78609506 , -0.61815876 , -0.34575712, XYZ) ; 
        setVector (13 , -0.34359558 , -0.11886104 , 0.99369402, XYZ) ; 
        setVector (14 , -0.99988003 , 0.21349205 , 0.27253397, XYZ) ; 
        setVector (15 , 0.34384765 , 0.70782462 , 0.70736609, XYZ) ; 
        setVector (16 , -0.11876178 , -0.99353168 , 0.34409900, XYZ) ; 
        setVector (17 , 0.27114748 , 0.65396056 , -0.78640214, XYZ) ; 
        setVector (18 , 0.27207095 , 0.99999945 , 0.21352327, XYZ) ; 
        setVector (19 , 0.70768179 , -0.70751268 , -0.34384005, XYZ) ; 
        setVector (20 , 0.34717058 , -0.78646197 , 0.61689857, XYZ) ; 
        
    case 22  
        setVector (0 , 0.76370272 , -0.76444935 , 0.00085437, XYZ) ; 
        setVector (1 , -1.00000000 , 0.29017427 , 0.28883263, XYZ) ; 
        setVector (2 , -0.62507853 , 0.62497202 , 0.62154019, XYZ) ; 
        setVector (3 , -0.76445467 , -0.76369788 , -0.00001732, XYZ) ; 
        setVector (4 , 0.76369779 , 0.00115879 , 0.76445388, XYZ) ; 
        setVector (5 , -0.29075120 , 0.28997165 , -0.99950267, XYZ) ; 
        setVector (6 , 0.29035056 , -0.99963727 , -0.28990910, XYZ) ; 
        setVector (7 , -0.00010358 , -0.76511552 , -0.76303579, XYZ) ; 
        setVector (8 , -0.29056661 , 0.99882560 , -0.29247881, XYZ) ; 
        setVector (9 , -0.62497971 , -0.62318885 , -0.62342715, XYZ) ; 
        setVector (10 , -0.29168269 , -0.99881045 , 0.29141772, XYZ) ; 
        setVector (11 , -0.76449178 , 0.00188046 , 0.76365842, XYZ) ; 
        setVector (12 , 0.29233949 , -0.29396788 , -0.99787069, XYZ) ; 
        setVector (13 , -0.62285514 , -0.62349150 , 0.62524813, XYZ) ; 
        setVector (14 , 0.29020872 , 0.29235739 , 0.99896519, XYZ) ; 
        setVector (15 , -0.00101588 , -0.76302239 , 0.76512822, XYZ) ; 
        setVector (16 , 0.99790277 , -0.29273283 , 0.29346710, XYZ) ; 
        setVector (17 , -0.62241799 , 0.62381750 , -0.62535827, XYZ) ; 
        setVector (18 , 0.29155135 , 0.28851479 , -0.99969116, XYZ) ; 
        setVector (19 , -0.29223558 , -0.99874512 , -0.29108752, XYZ) ; 
        setVector (20 , 0.99998534 , 0.28958005 , 0.28947902, XYZ) ; 
        setVector (21 , 0.99798524 , 0.29268877 , -0.29323053, XYZ) ; 
        
    case 23  
        setVector (0 , -0.99208446 , -0.33288968 , -0.13066926, XYZ) ; 
        setVector (1 , 0.28942892 , 0.76912786 , -0.66090453, XYZ) ; 
        setVector (2 , 0.99518533 , -0.16839710 , 0.30556530, XYZ) ; 
        setVector (3 , 0.32403864 , 0.61931717 , 0.78966245, XYZ) ; 
        setVector (4 , 0.79638871 , 0.64325906 , -0.25318841, XYZ) ; 
        setVector (5 , -0.33465970 , -0.01119964 , -0.99999951, XYZ) ; 
        setVector (6 , -0.64725520 , -0.75989219 , -0.34021486, XYZ) ; 
        setVector (7 , 0.68876551 , -0.30639456 , 0.73745917, XYZ) ; 
        setVector (8 , -0.17619324 , -0.28474114 , 1.00000000, XYZ) ; 
        setVector (9 , -0.23654096 , 0.75930768 , 0.69254733, XYZ) ; 
        setVector (10 , -0.12246971 , 0.31164611 , 0.99999972, XYZ) ; 
        setVector (11 , 0.37715274 , 0.97526904 , -0.13684913, XYZ) ; 
        setVector (12 , 0.76388254 , 0.27472267 , 0.67315116, XYZ) ; 
        setVector (13 , 0.90767930 , -0.53274730 , -0.06648446, XYZ) ; 
        setVector (14 , 0.69468776 , 0.37701701 , -0.69813227, XYZ) ; 
        setVector (15 , 0.23852769 , -0.63119697 , 0.81044216, XYZ) ; 
        setVector (16 , -0.61141733 , 0.11898649 , 0.85095983, XYZ) ; 
        setVector (17 , 0.09702081 , 0.99756117 , 0.32799427, XYZ) ; 
        setVector (18 , -0.64490878 , 0.74242382 , -0.38081638, XYZ) ; 
        setVector (19 , 0.46119039 , -0.93440777 , -0.16219470, XYZ) ; 
        setVector (20 , -0.71702774 , 0.53908780 , 0.55441603, XYZ) ; 
        setVector (21 , -0.99490049 , 0.00602820 , 0.34965449, XYZ) ; 
        setVector (22 , -0.14802589 , 0.99538276 , -0.31531425, XYZ) ; 
        
    case 24  
        setVector (0 , -0.78948442 , 0.61667435 , -0.33390887, XYZ) ; 
        setVector (1 , -0.73929581 , 0.73188269 , 0.18126635, XYZ) ; 
        setVector (2 , 0.90033068 , 0.45734275 , 0.30872370, XYZ) ; 
        setVector (3 , -0.02542745 , 0.65158805 , 0.83057481, XYZ) ; 
        setVector (4 , 0.18549377 , 0.99807745 , 0.29069143, XYZ) ; 
        setVector (5 , -0.27068117 , -0.20522010 , -0.99984222, XYZ) ; 
        setVector (6 , -0.48878951 , -0.68916893 , -0.63340276, XYZ) ; 
        setVector (7 , 0.30851674 , -0.90023990 , -0.45766100, XYZ) ; 
        setVector (8 , -0.32207305 , 1.00000000 , -0.10647534, XYZ) ; 
        setVector (9 , -0.45722117 , -0.30879471 , 0.90036808, XYZ) ; 
        setVector (10 , 0.99999670 , -0.10665125 , 0.32202510, XYZ) ; 
        setVector (11 , -0.33538426 , 0.78924585 , -0.61617889, XYZ) ; 
        setVector (12 , 0.63362317 , -0.48859157 , -0.68910667, XYZ) ; 
        setVector (13 , 0.61715530 , -0.33480464 , 0.78872887, XYZ) ; 
        setVector (14 , 0.68828844 , 0.63402884 , -0.48921826, XYZ) ; 
        setVector (15 , 0.73196222 , 0.18068693 , 0.73935890, XYZ) ; 
        setVector (16 , 0.83051078 , 0.02577683 , -0.65165593, XYZ) ; 
        setVector (17 , 0.65164356 , 0.83053220 , 0.02539665, XYZ) ; 
        setVector (18 , -0.99999823 , 0.27023447 , 0.20504859, XYZ) ; 
        setVector (19 , 0.20501012 , 0.99972020 , -0.27129028, XYZ) ; 
        setVector (20 , -0.29111135 , 0.18539408 , 0.99797358, XYZ) ; 
        setVector (21 , -0.99771567 , -0.29176316 , 0.18575722, XYZ) ; 
        setVector (22 , -0.10746987 , 0.32175648 , -0.99999552, XYZ) ; 
        setVector (23 , -0.17934432 , -0.73888877 , 0.73276671, XYZ) ; 
        
    case 25  
        setVector (0 , -0.31749379 , -0.72811714 , -0.70217819, XYZ) ; 
        setVector (1 , -0.29928754 , 0.18564641 , -0.99998673, XYZ) ; 
        setVector (2 , 0.73101598 , -0.73860658 , -0.20996917, XYZ) ; 
        setVector (3 , -0.61445026 , -0.60760988 , 0.61422488, XYZ) ; 
        setVector (4 , 0.27164498 , 0.99420282 , 0.24855752, XYZ) ; 
        setVector (5 , 0.98322221 , 0.32141863 , -0.23232572, XYZ) ; 
        setVector (6 , 0.72778140 , 0.76509728 , -0.09471786, XYZ) ; 
        setVector (7 , 0.09481549 , 0.76517177 , -0.72769037, XYZ) ; 
        setVector (8 , 0.70225471 , -0.72835403 , 0.31678041, XYZ) ; 
        setVector (9 , 0.23244622 , 0.32166573 , -0.98311291, XYZ) ; 
        setVector (10 , 0.26333714 , -0.99277653 , -0.26279157, XYZ) ; 
        setVector (11 , -0.40000535 , 0.63799945 , -0.74629988, XYZ) ; 
        setVector (12 , -0.74542079 , -0.11157958 , 0.74559301, XYZ) ; 
        setVector (13 , 0.24838866 , -0.99414406 , 0.27201426, XYZ) ; 
        setVector (14 , 0.66658381 , -0.48399477 , -0.66740256, XYZ) ; 
        setVector (15 , -0.74651052 , -0.63781780 , -0.39990197, XYZ) ; 
        setVector (16 , -0.66033148 , -0.24808273 , -0.79147229, XYZ) ; 
        setVector (17 , -0.98175708 , 0.19560994 , 0.34914306, XYZ) ; 
        setVector (18 , -0.14690511 , -0.32019653 , -0.99995208, XYZ) ; 
        setVector (19 , 0.27899397 , 0.98410674 , -0.27876043, XYZ) ; 
        setVector (20 , -1.00000000 , 0.32019518 , -0.14658150, XYZ) ; 
        setVector (21 , 0.34858878 , -0.19547485 , -0.98198093, XYZ) ; 
        setVector (22 , 0.99999376 , 0.18537256 , 0.29943376, XYZ) ; 
        setVector (23 , -0.79149075 , 0.24858440 , -0.66012065, XYZ) ; 
        setVector (24 , -0.20905419 , 0.73914756 , 0.73073136, XYZ) ; 
        
    case 26  
        setVector (0 , 0.67107328 , 0.69162324 , -0.42257061, XYZ) ; 
        setVector (1 , -0.71337850 , -0.06865400 , -0.77047110, XYZ) ; 
        setVector (2 , 0.25683989 , -0.72624391 , 0.71683403, XYZ) ; 
        setVector (3 , -0.36915155 , -0.89919277 , -0.40302292, XYZ) ; 
        setVector (4 , 0.99132311 , 0.27602585 , -0.21985496, XYZ) ; 
        setVector (5 , 0.31796291 , 0.07840637 , -1.00000000, XYZ) ; 
        setVector (6 , 0.47795683 , -0.73190714 , -0.58576205, XYZ) ; 
        setVector (7 , 0.25027875 , 0.25295424 , 0.99026394, XYZ) ; 
        setVector (8 , -0.97180553 , 0.28497994 , -0.28570688, XYZ) ; 
        setVector (9 , -0.26796066 , 0.99998589 , -0.18834352, XYZ) ; 
        setVector (10 , 0.63049809 , 0.54242383 , 0.64459020, XYZ) ; 
        setVector (11 , -0.32539080 , -0.99842869 , 0.06714875, XYZ) ; 
        setVector (12 , -0.97851822 , -0.22814859 , -0.31256725, XYZ) ; 
        setVector (13 , 0.19421096 , -0.26374385 , 0.99998463, XYZ) ; 
        setVector (14 , -0.67997369 , 0.73688125 , -0.31920178, XYZ) ; 
        setVector (15 , -0.28471080 , 0.39407471 , 0.93321640, XYZ) ; 
        setVector (16 , 0.78874322 , 0.68869575 , 0.10406856, XYZ) ; 
        setVector (17 , -0.16462864 , -0.90429332 , 0.51224894, XYZ) ; 
        setVector (18 , 0.07409987 , 0.71817772 , 0.76549196, XYZ) ; 
        setVector (19 , 0.65367680 , -0.40314071 , 0.71932759, XYZ) ; 
        setVector (20 , -0.16001770 , 0.99733816 , 0.29488796, XYZ) ; 
        setVector (21 , 0.75313138 , 0.23369086 , -0.69672783, XYZ) ; 
        setVector (22 , 0.75732349 , -0.27971501 , -0.67488415, XYZ) ; 
        setVector (23 , 0.29763440 , 0.55061255 , -0.84586497, XYZ) ; 
        setVector (24 , -0.99997723 , 0.23844002 , 0.22458822, XYZ) ; 
        setVector (25 , 0.75291097 , -0.71469088 , -0.17201741, XYZ) ; 
        
    case 27  
        setVector (0 , -0.72696804 , 0.25060480 , -0.70482408, XYZ) ; 
        setVector (1 , 0.10967734 , 0.32636336 , 0.98464212, XYZ) ; 
        setVector (2 , 0.40471230 , -0.86835126 , -0.41259704, XYZ) ; 
        setVector (3 , 0.99999447 , -0.22201429 , -0.19693399, XYZ) ; 
        setVector (4 , -0.16730403 , -0.24509522 , 0.99999998, XYZ) ; 
        setVector (5 , -0.41748038 , 0.18659441 , 0.93752597, XYZ) ; 
        setVector (6 , 0.85208957 , 0.53518742 , 0.27491830, XYZ) ; 
        setVector (7 , 0.16271958 , -0.66649853 , -0.78572535, XYZ) ; 
        setVector (8 , 0.83601014 , -0.02097365 , -0.62346566, XYZ) ; 
        setVector (9 , -0.68835717 , -0.28212564 , -0.73118521, XYZ) ; 
        setVector (10 , -0.99479505 , -0.07855673 , -0.30376623, XYZ) ; 
        setVector (11 , -0.59231211 , 0.73049934 , -0.45121986, XYZ) ; 
        setVector (12 , 0.39124034 , 0.68531794 , 0.68215290, XYZ) ; 
        setVector (13 , 0.54321678 , 0.87025777 , 0.18875706, XYZ) ; 
        setVector (14 , 0.63011658 , 0.37744908 , -0.74063997, XYZ) ; 
        setVector (15 , 0.69114268 , -0.49864075 , -0.60144948, XYZ) ; 
        setVector (16 , -0.12398814 , 0.89641952 , -0.51876898, XYZ) ; 
        setVector (17 , 0.25053438 , 0.99329938 , -0.19659893, XYZ) ; 
        setVector (18 , 0.97319724 , 0.31154695 , -0.20949441, XYZ) ; 
        setVector (19 , 0.69185745 , 0.71452198 , -0.31441039, XYZ) ; 
        setVector (20 , -0.29115835 , 1.00000000 , -0.05735061, XYZ) ; 
        setVector (21 , 0.25231128 , -0.54327622 , 0.85396267, XYZ) ; 
        setVector (22 , -0.08425717 , -0.97678566 , -0.35616399, XYZ) ; 
        setVector (23 , -0.28745004 , -0.71118806 , 0.70685663, XYZ) ; 
        setVector (24 , -0.35063634 , 0.06997215 , -0.97990833, XYZ) ; 
        setVector (25 , 0.92406818 , -0.40364609 , 0.26688971, XYZ) ; 
        setVector (26 , -0.75611793 , 0.71283782 , 0.09061016, XYZ) ; 
        
    case 28  
        setVector (0 , 0.09944788 , -1.00000000 , -0.27162349, XYZ) ; 
        setVector (1 , 0.06107997 , 0.68092435 , -0.78503533, XYZ) ; 
        setVector (2 , 0.53333951 , -0.81114082 , -0.37585735, XYZ) ; 
        setVector (3 , -0.98494172 , -0.32182888 , -0.09992587, XYZ) ; 
        setVector (4 , -0.50285931 , -0.27706828 , -0.86835182, XYZ) ; 
        setVector (5 , 0.40850157 , 0.87612976 , -0.38625421, XYZ) ; 
        setVector (6 , 0.23974087 , -0.24938200 , 0.98183610, XYZ) ; 
        setVector (7 , -0.99990919 , 0.06369847 , 0.28247710, XYZ) ; 
        setVector (8 , 0.72997689 , 0.74157252 , -0.02954899, XYZ) ; 
        setVector (9 , -0.37573146 , 0.72605907 , -0.64446358, XYZ) ; 
        setVector (10 , -0.75655932 , 0.41917629 , 0.57929132, XYZ) ; 
        setVector (11 , 0.73253962 , 0.07144199 , -0.73617318, XYZ) ; 
        setVector (12 , -0.04010688 , 0.99982603 , -0.28706890, XYZ) ; 
        setVector (13 , 0.86048519 , 0.15821127 , 0.56409542, XYZ) ; 
        setVector (14 , -0.87342290 , 0.55812587 , 0.09642176, XYZ) ; 
        setVector (15 , -0.99997967 , 0.17768509 , -0.22833720, XYZ) ; 
        setVector (16 , 0.86417958 , 0.43004043 , -0.38977954, XYZ) ; 
        setVector (17 , -0.76697734 , 0.57960410 , -0.39934202, XYZ) ; 
        setVector (18 , 0.66177965 , -0.25407466 , 0.76234045, XYZ) ; 
        setVector (19 , 0.26182662 , 0.74636063 , 0.67680265, XYZ) ; 
        setVector (20 , 0.51861527 , 0.53147626 , -0.72954807, XYZ) ; 
        setVector (21 , 0.41540063 , -0.28390680 , -0.91132236, XYZ) ; 
        setVector (22 , -0.48936239 , 0.91033498 , -0.12443425, XYZ) ; 
        setVector (23 , -0.34859885 , -0.96650146 , -0.16740060, XYZ) ; 
        setVector (24 , 0.06097322 , 0.28275751 , 0.99999983, XYZ) ; 
        setVector (25 , -0.20815614 , 0.70344146 , 0.73858672, XYZ) ; 
        setVector (26 , 0.24279892 , 0.19092247 , -0.99411594, XYZ) ; 
        setVector (27 , -0.68497145 , -0.62915507 , -0.46759727, XYZ) ; 
        
    case 30  
        setVector (0 , -0.37813910 , 0.39432155 , -0.88358565, XYZ) ; 
        setVector (1 , -0.14754600 , -0.85148967 , -0.57653950, XYZ) ; 
        setVector (2 , -0.99972185 , 0.08487479 , -0.26935992, XYZ) ; 
        setVector (3 , 0.85130361 , -0.57636940 , 0.14927410, XYZ) ; 
        setVector (4 , -0.74830149 , -0.71924849 , 0.04391769, XYZ) ; 
        setVector (5 , 0.28353074 , -0.20653661 , -0.97783190, XYZ) ; 
        setVector (6 , 0.69128845 , -0.34340784 , -0.69526514, XYZ) ; 
        setVector (7 , 0.69625402 , -0.69076339 , -0.34245956, XYZ) ; 
        setVector (8 , 0.07528089 , -0.99999399 , 0.27119567, XYZ) ; 
        setVector (9 , -0.38571644 , -0.58215335 , 0.76910505, XYZ) ; 
        setVector (10 , 0.99995502 , 0.27151028 , -0.07466196, XYZ) ; 
        setVector (11 , 0.76828535 , -0.38487577 , 0.58378982, XYZ) ; 
        setVector (12 , 0.20668878 , -0.97802672 , -0.28274682, XYZ) ; 
        setVector (13 , 0.47647922 , -0.76307707 , 0.51950285, XYZ) ; 
        setVector (14 , -0.26827714 , -1.00000000 , -0.08502736, XYZ) ; 
        setVector (15 , 0.08567737 , 0.26808752 , -0.99999537, XYZ) ; 
        setVector (16 , -0.51843280 , -0.47735027 , -0.76326039, XYZ) ; 
        setVector (17 , -0.76348279 , -0.51868950 , 0.47671534, XYZ) ; 
        setVector (18 , -0.57726463 , -0.14759357 , 0.85098999, XYZ) ; 
        setVector (19 , 0.27118843 , 0.07537287 , 0.99998902, XYZ) ; 
        setVector (20 , -0.34368870 , 0.69592276 , 0.69048669, XYZ) ; 
        setVector (21 , -0.50917086 , 0.90353853 , -0.05971129, XYZ) ; 
        setVector (22 , -0.97824310 , 0.28305186 , 0.20524221, XYZ) ; 
        setVector (23 , -0.90416493 , -0.05889779 , 0.50815264, XYZ) ; 
        setVector (24 , 0.05868784 , 0.50898262 , 0.90371163, XYZ) ; 
        setVector (25 , 0.39431019 , 0.88307679 , -0.37933777, XYZ) ; 
        setVector (26 , -0.88292213 , -0.37965226 , -0.39435384, XYZ) ; 
        setVector (27 , -0.71869086 , -0.04442942 , -0.74880688, XYZ) ; 
        setVector (28 , -0.58233607 , -0.76920803 , -0.38523498, XYZ) ; 
        setVector (29 , -0.04294536 , 0.74899650 , -0.71858347, XYZ) ; 
        
    case 32  
        setVector (0 , -0.44213760 , -0.51075569 , 0.79115308, XYZ) ; 
        setVector (1 , -0.48226067 , 0.81771156 , 0.42550285, XYZ) ; 
        setVector (2 , -0.33701523 , -0.51032270 , -0.84158880, XYZ) ; 
        setVector (3 , -0.83154789 , -0.53095924 , 0.32998579, XYZ) ; 
        setVector (4 , -0.01499202 , -0.97940815 , 0.35044992, XYZ) ; 
        setVector (5 , -0.05274517 , 0.73249636 , 0.73684955, XYZ) ; 
        setVector (6 , -0.29290233 , -0.99806365 , -0.01890473, XYZ) ; 
        setVector (7 , -0.41282381 , -0.82377350 , -0.48296372, XYZ) ; 
        setVector (8 , 0.78397376 , 0.53076824 , 0.43121971, XYZ) ; 
        setVector (9 , -0.68116174 , -0.23608921 , -0.75004052, XYZ) ; 
        setVector (10 , -0.44963617 , 0.50432305 , 0.79105360, XYZ) ; 
        setVector (11 , -0.06686161 , 0.27894399 , 1.00000000, XYZ) ; 
        setVector (12 , 0.68524341 , -0.22093975 , 0.75093763, XYZ) ; 
        setVector (13 , -0.77156437 , -0.23403834 , 0.65741517, XYZ) ; 
        setVector (14 , -0.99015922 , -0.31312859 , -0.06176915, XYZ) ; 
        setVector (15 , 0.34724857 , -0.50417362 , 0.84113472, XYZ) ; 
        setVector (16 , 0.99415822 , -0.29972901 , 0.06397016, XYZ) ; 
        setVector (17 , 0.42420225 , -0.81796427 , 0.48297737, XYZ) ; 
        setVector (18 , -0.31622437 , -0.00331821 , -0.99109604, XYZ) ; 
        setVector (19 , -0.06079133 , -0.28034030 , 0.99999698, XYZ) ; 
        setVector (20 , 0.30787941 , -0.99354009 , 0.01919865, XYZ) ; 
        setVector (21 , -0.77410491 , 0.22327824 , 0.65817068, XYZ) ; 
        setVector (22 , 0.99999089 , 0.00592504 , -0.28681587, XYZ) ; 
        setVector (23 , 0.02908737 , -0.97930217 , -0.34985913, XYZ) ; 
        setVector (24 , 0.46970641 , 0.82445768 , -0.42652743, XYZ) ; 
        setVector (25 , -0.95569432 , -0.00773421 , -0.41093646, XYZ) ; 
        setVector (26 , -0.04072099 , -0.73221845 , 0.73788763, XYZ) ; 
        setVector (27 , -0.83931658 , 0.51890180 , 0.32949786, XYZ) ; 
        setVector (28 , -0.79042377 , 0.51838543 , -0.43449630, XYZ) ; 
        setVector (29 , -0.70744739 , 0.76139454 , -0.04557167, XYZ) ; 
        setVector (30 , 0.69563265 , 0.77230605 , 0.04380415, XYZ) ; 
        setVector (31 , 0.43990030 , 0.00470587 , -0.94273316, XYZ) ; 
        
        
        %//Made of the first part of the 60 case (combine 31 with 29)
    case 31  
        setVector (0 , -0.51711309 , -0.47983043 , 0.73648370, XYZ) ; 
        setVector (1 , -0.75808487 , 0.49330514 , -0.47117809, XYZ) ; 
        setVector (2 , 0.80234240 , -0.60886478 , 0.15994303, XYZ) ; 
        setVector (3 , -0.01973501 , 0.48934133 , -0.89454291, XYZ) ; 
        setVector (4 , 0.33516758 , 0.93013284 , -0.25013402, XYZ) ; 
        setVector (5 , -0.63448454 , 0.31901719 , -0.73192133, XYZ) ; 
        setVector (6 , 0.08116896 , -0.18293329 , 0.99999922, XYZ) ; 
        setVector (7 , -0.10445113 , 0.66242333 , 0.76833373, XYZ) ; 
        setVector (8 , 0.26766319 , -0.88671109 , 0.42679184, XYZ) ; 
        setVector (9 , 0.67016636 , -0.42004243 , -0.64381118, XYZ) ; 
        setVector (10 , 0.15391958 , -0.72580995 , 0.69968571, XYZ) ; 
        setVector (11 , 0.41367573 , -0.09455360 , 0.92735292, XYZ) ; 
        setVector (12 , 0.88398617 , 0.46713278 , -0.20101452, XYZ) ; 
        setVector (13 , -0.21753901 , -0.66282962 , 0.74389858, XYZ) ; 
        setVector (14 , 0.93271228 , 0.37641775 , 0.16854941, XYZ) ; 
        setVector (15 , 0.69488384 , -0.64060372 , -0.38316412, XYZ) ; 
        setVector (16 , 0.88211118 , -0.49143385 , -0.14291280, XYZ) ; 
        setVector (17 , -0.69509295 , -0.01487792 , -0.74610714, XYZ) ; 
        setVector (18 , 0.87347998 , -0.16054820 , 0.50130672, XYZ) ; 
        setVector (19 , -0.75283496 , -0.54879866 , -0.41486259, XYZ) ; 
        setVector (20 , 0.88371946 , -0.28861432 , -0.41927689, XYZ) ; 
        setVector (21 , -0.88840244 , -0.18916095 , -0.46369242, XYZ) ; 
        setVector (22 , -0.55232725 , -0.15049253 , 0.84400121, XYZ) ; 
        setVector (23 , 0.93679074 , 0.06730489 , -0.39742244, XYZ) ; 
        setVector (24 , 0.94961311 , -0.31584415 , 0.19628767, XYZ) ; 
        setVector (25 , 0.64277964 , 0.76457665 , -0.20568982, XYZ) ; 
        setVector (26 , -0.34381087 , 0.43058758 , -0.85816073, XYZ) ; 
        setVector (27 , 0.99971116 , -0.16450303 , -0.11648077, XYZ) ; 
        setVector (28 , 0.99821503 , 0.18210187 , -0.10225989, XYZ) ; 
        setVector (29 , 0.46606173 , 0.58200226 , 0.69578104, XYZ) ; 
        setVector (30 , -0.13874523 , -0.14422617 , -1.00000000, XYZ) ; 
        
        
        %//Made of the first part of the 60 case (combine 31 with 29)
    case 29 
        setVector (0 , 0.76937618 , 0.66306062 , 0.09199093, XYZ) ; 
        setVector (1 , 0.99999924 , 0.01670781 , 0.19943366, XYZ) ; 
        setVector (2 , -0.04021492 , 0.99903292 , -0.20091641, XYZ) ; 
        setVector (3 , -0.19526217 , -0.76591725 , -0.64443377, XYZ) ; 
        setVector (4 , 0.51096926 , 0.77929365 , 0.41432263, XYZ) ; 
        setVector (5 , 0.64449110 , -0.78747432 , -0.06757843, XYZ) ; 
        setVector (6 , -0.11639260 , 0.99961185 , 0.16516763, XYZ) ; 
        setVector (7 , -0.41472436 , 0.88677651 , 0.28580159, XYZ) ; 
        setVector (8 , -0.68663556 , -0.34283814 , -0.67159887, XYZ) ; 
        setVector (9 , 0.38949705 , -0.48693386 , -0.80699374, XYZ) ; 
        setVector (10 , 0.50524772 , 0.88074278 , 0.09522775, XYZ) ; 
        setVector (11 , 0.11717751 , -0.88352002 , -0.49569470, XYZ) ; 
        setVector (12 , 0.76094012 , -0.07899991 , -0.67437421, XYZ) ; 
        setVector (13 , -0.78334360 , -0.26876947 , 0.59513629, XYZ) ; 
        setVector (14 , -0.44028061 , -0.25290552 , -0.88444514, XYZ) ; 
        setVector (15 , -0.20589133 , -0.93262713 , -0.35758470, XYZ) ; 
        setVector (16 , 0.13436626 , -0.33491810 , -0.95384852, XYZ) ; 
        setVector (17 , -0.27000378 , -0.32593888 , 0.92785410, XYZ) ; 
        setVector (18 , 0.44543086 , 0.75894353 , -0.51541002, XYZ) ; 
        setVector (19 , -0.21816264 , 0.00266380 , 0.99621754, XYZ) ; 
        setVector (20 , -0.36644717 , 0.94866592 , -0.07616346, XYZ) ; 
        setVector (21 , -0.19651463 , -0.99998147 , -0.03834686, XYZ) ; 
        setVector (22 , -0.56797954 , 0.78726452 , -0.31251438, XYZ) ; 
        setVector (23 , -0.42364378 , 0.72459950 , 0.57925205, XYZ) ; 
        setVector (24 , -0.48425159 , 0.64709189 , -0.62195169, XYZ) ; 
        setVector (25 , -0.48837442 , 0.17825689 , 0.87736329, XYZ) ; 
        setVector (26 , -0.71876309 , -0.55582022 , 0.46313598, XYZ) ; 
        setVector (27 , -0.09360757 , -0.88285078 , 0.50186008, XYZ) ; 
        setVector (28 , -0.17722160 , -0.48544349 , -0.87919767, XYZ) ; 
        
    case 62  
        switch (mdr_lBerhens_flag) 
            case 1 
                setVector (0 , 0.11728079 , 0.35366914 , -0.93839532, XYZ) ; 
                setVector (1 , 0.31607521 , -0.94671839 , 0.15245711, XYZ) ; 
                setVector (2 , 0.36933501 , 0.93730292 , 0.06691264, XYZ) ; 
                setVector (3 , 0.29153602 , 0.88119900 , -0.39738833, XYZ) ; 
                setVector (4 , -0.67523831 , 0.10640151 , 0.74307090, XYZ) ; 
                setVector (5 , 0.54827642 , -0.07114892 , 0.84483917, XYZ) ; 
                setVector (6 , -0.89481670 , 0.45706603 , 0.09907637, XYZ) ; 
                setVector (7 , 0.83615565 , 0.56532838 , 0.02588377, XYZ) ; 
                setVector (8 , 0.48290356 , -0.88555408 , -0.04495049, XYZ) ; 
                setVector (9 , -0.15517455 , 0.82728763 , -0.55761855, XYZ) ; 
                setVector (10 , -0.67122559 , 0.75416560 , -0.01062440, XYZ) ; 
                setVector (11 , 0.10228347 , 0.96535389 , -0.27758312, XYZ) ; 
                setVector (12 , 0.60055740 , -0.77800415 , -0.23122019, XYZ) ; 
                setVector (13 , -0.93674534 , -0.06077828 , 0.37180182, XYZ) ; 
                setVector (14 , -0.23854879 , 0.50589374 , 0.84058814, XYZ) ; 
                setVector (15 , -0.96156476 , 0.04086344 , -0.30519799, XYZ) ; 
                setVector (16 , 0.55131399 , 0.84527219 , -0.03146783, XYZ) ; 
                setVector (17 , -0.51286493 , 0.73632785 , 0.46283182, XYZ) ; 
                setVector (18 , -0.48011075 , 0.21134064 , 0.86269996, XYZ) ; 
                setVector (19 , -0.76943823 , -0.46845260 , 0.45600373, XYZ) ; 
                setVector (20 , -0.52137600 , 0.84365941 , -0.18928363, XYZ) ; 
                setVector (21 , 0.05512599 , -0.99999845 , -0.12801033, XYZ) ; 
                setVector (22 , -0.00579203 , -0.72430344 , 0.70340131, XYZ) ; 
                setVector (23 , -0.99999412 , 0.10954384 , -0.08622259, XYZ) ; 
                setVector (24 , 0.45597497 , 0.77594892 , 0.45761613, XYZ) ; 
                setVector (25 , 0.92660801 , -0.17527261 , -0.36069313, XYZ) ; 
                setVector (26 , 0.36866697 , 0.05618976 , 0.93826957, XYZ) ; 
                setVector (27 , -0.59404488 , -0.38232515 , -0.72136023, XYZ) ; 
                setVector (28 , 0.57958003 , -0.55616789 , -0.61170801, XYZ) ; 
                setVector (29 , -0.40703104 , -0.49432543 , -0.78063469, XYZ) ; 
                setVector (30 , 0.73694541 , 0.03195869 , 0.68942913, XYZ) ; 
                setVector (31 , 0.04188097 , -0.34754252 , -0.94703886, XYZ) ; 
                setVector (32 , 0.52509642 , 0.83280373 , 0.22390649, XYZ) ; 
                setVector (33 , -0.32737549 , 0.17606154 , -0.93874919, XYZ) ; 
                setVector (34 , -0.95549538 , -0.17789172 , -0.27350600, XYZ) ; 
                setVector (35 , -0.70714956 , -0.70976142 , 0.12490256, XYZ) ; 
                setVector (36 , -0.77166739 , 0.62802140 , 0.17187490, XYZ) ; 
                setVector (37 , -0.09708335 , 0.27190798 , -0.96750364, XYZ) ; 
                setVector (38 , 0.16512818 , 0.58631795 , -0.80522442, XYZ) ; 
                setVector (39 , 0.16917715 , 0.98105027 , 0.16835047, XYZ) ; 
                setVector (40 , 0.27559126 , -0.96742099 , -0.08699719, XYZ) ; 
                setVector (41 , -0.30114896 , 0.81582530 , 0.51299200, XYZ) ; 
                setVector (42 , 0.83332841 , 0.22426646 , -0.52410948, XYZ) ; 
                setVector (43 , -0.37260530 , -0.62908013 , 0.69630877, XYZ) ; 
                setVector (44 , -0.29530963 , 0.05248661 , 0.96408495, XYZ) ; 
                setVector (45 , 0.71977694 , 0.48479913 , 0.51605559, XYZ) ; 
                setVector (46 , -0.54365207 , 0.61005360 , -0.59304254, XYZ) ; 
                setVector (47 , -0.69865442 , -0.71707579 , -0.13079274, XYZ) ; 
                setVector (48 , -0.45501518 , -0.75246073 , 0.49617180, XYZ) ; 
                setVector (49 , 0.68833696 , 0.35081904 , -0.65003127, XYZ) ; 
                setVector (50 , -0.01836275 , 0.55779468 , 0.84139784, XYZ) ; 
                setVector (51 , -0.49259456 , 0.30276002 , -0.82771335, XYZ) ; 
                setVector (52 , 0.99947733 , 0.13141500 , 0.05654714, XYZ) ; 
                setVector (53 , 0.97890374 , 0.18044735 , -0.16914088, XYZ) ; 
                setVector (54 , 0.38100374 , -0.64802674 , -0.67403258, XYZ) ; 
                setVector (55 , -0.31459045 , -0.18524336 , 0.94135018, XYZ) ; 
                setVector (56 , -0.84563576 , 0.47282842 , -0.28417577, XYZ) ; 
                setVector (57 , -0.44804642 , 0.43752139 , 0.79199234, XYZ) ; 
                setVector (58 , -0.85888244 , -0.34254696 , -0.40546882, XYZ) ; 
                setVector (59 , 0.94544041 , -0.33896507 , 0.10328377, XYZ) ; 
                setVector (60 , 0.82365210 , -0.57871726 , 0.07814067, XYZ) ; 
                setVector (61 , -0.58146572 , -0.16087883 , -0.80959128, XYZ) ; 
                
            case 2 
                setVector (0 , 0.09391635 , -0.99953888 , 0.10735069, XYZ) ; 
                setVector (1 , 0.63224273 , -0.41381821 , 0.66966114, XYZ) ; 
                setVector (2 , -0.20718755 , 0.64074087 , -0.75229441, XYZ) ; 
                setVector (3 , -0.26764933 , 0.28637466 , 0.93047075, XYZ) ; 
                setVector (4 , -0.51674907 , -0.24552082 , 0.83193290, XYZ) ; 
                setVector (5 , -0.59589397 , -0.58351484 , 0.56907226, XYZ) ; 
                setVector (6 , 0.83106571 , -0.00307204 , -0.57336094, XYZ) ; 
                setVector (7 , 0.53908412 , -0.75009975 , 0.40762863, XYZ) ; 
                setVector (8 , -0.96844911 , 0.25255969 , 0.13320040, XYZ) ; 
                setVector (9 , 0.92095295 , -0.25717229 , 0.32423833, XYZ) ; 
                setVector (10 , -0.40228632 , 0.88003300 , 0.28832284, XYZ) ; 
                setVector (11 , -0.54408209 , -0.59916884 , -0.60365040, XYZ) ; 
                setVector (12 , 0.66595058 , 0.66886148 , -0.35854785, XYZ) ; 
                setVector (13 , 0.82749683 , 0.53293284 , -0.22506445, XYZ) ; 
                setVector (14 , -0.05852276 , -0.87897404 , 0.49335803, XYZ) ; 
                setVector (15 , 0.69969028 , -0.68508080 , 0.24600860, XYZ) ; 
                setVector (16 , -0.18879143 , -0.39830936 , -0.90836660, XYZ) ; 
                setVector (17 , 0.89925032 , 0.32437553 , -0.32488735, XYZ) ; 
                setVector (18 , 0.13929643 , 0.99847821 , -0.05531896, XYZ) ; 
                setVector (19 , 0.51630296 , 0.82410991 , -0.27147105, XYZ) ; 
                setVector (20 , 0.31897279 , 0.90268957 , 0.32067169, XYZ) ; 
                setVector (21 , 0.75151930 , -0.43194278 , -0.51775148, XYZ) ; 
                setVector (22 , 0.09780702 , 0.11464404 , -0.99835512, XYZ) ; 
                setVector (23 , -0.20163222 , -0.59210016 , -0.79258077, XYZ) ; 
                setVector (24 , -0.75351596 , -0.25969023 , -0.61983635, XYZ) ; 
                setVector (25 , 0.69770139 , -0.19626905 , 0.70293218, XYZ) ; 
                setVector (26 , 0.71140724 , -0.55705784 , 0.45056490, XYZ) ; 
                setVector (27 , 0.04852382 , -0.51248977 , 0.86857475, XYZ) ; 
                setVector (28 , -0.87547118 , -0.10732945 , -0.49137870, XYZ) ; 
                setVector (29 , 0.50637402 , 0.01941685 , -0.87328733, XYZ) ; 
                setVector (30 , -0.93420222 , -0.37609197 , 0.07241180, XYZ) ; 
                setVector (31 , -0.32062054 , -0.41646481 , 0.86208002, XYZ) ; 
                setVector (32 , 0.79781596 , -0.34643020 , 0.51273604, XYZ) ; 
                setVector (33 , 0.70219919 , -0.61275373 , -0.38842191, XYZ) ; 
                setVector (34 , -0.86536577 , 0.40726294 , 0.32357595, XYZ) ; 
                setVector (35 , -0.51006994 , -0.45584968 , 0.74259824, XYZ) ; 
                setVector (36 , -0.33940353 , -0.93566998 , 0.16955634, XYZ) ; 
                setVector (37 , 0.07214131 , -0.85417716 , -0.53347862, XYZ) ; 
                setVector (38 , -0.18429366 , 0.93773164 , 0.32575703, XYZ) ; 
                setVector (39 , 0.80778210 , 0.54038503 , 0.27366865, XYZ) ; 
                setVector (40 , -0.35610409 , 0.73583062 , -0.59259227, XYZ) ; 
                setVector (41 , -0.41289294 , 0.52358866 , -0.75815352, XYZ) ; 
                setVector (42 , -0.63899435 , 0.33204867 , 0.70770921, XYZ) ; 
                setVector (43 , 0.81170384 , -0.22290634 , -0.55755902, XYZ) ; 
                setVector (44 , 0.39626523 , 0.28227214 , 0.88471393, XYZ) ; 
                setVector (45 , -0.26415889 , 0.40902543 , -0.88450025, XYZ) ; 
                setVector (46 , 0.85342947 , -0.12018939 , 0.52596109, XYZ) ; 
                setVector (47 , 0.16301647 , -0.70416973 , -0.70497734, XYZ) ; 
                setVector (48 , 0.92609080 , 0.37096823 , 0.15543752, XYZ) ; 
                setVector (49 , -0.19927354 , -0.84588837 , -0.51398964, XYZ) ; 
                setVector (50 , -0.07055140 , -0.74246457 , -0.68058158, XYZ) ; 
                setVector (51 , -0.16600573 , -0.18475263 , -0.97863731, XYZ) ; 
                setVector (52 , 0.05171295 , 0.94218440 , 0.35921685, XYZ) ; 
                setVector (53 , 0.12843041 , -0.94296299 , 0.33726673, XYZ) ; 
                setVector (54 , -0.64968755 , -0.67488824 , -0.37663560, XYZ) ; 
                setVector (55 , -0.06914268 , 0.12104217 , 0.99999525, XYZ) ; 
                setVector (56 , 0.69219009 , 0.12604078 , -0.72416091, XYZ) ; 
                setVector (57 , -0.34213262 , 0.86954638 , -0.38243534, XYZ) ; 
                setVector (58 , 0.22282423 , 0.77181741 , -0.61161236, XYZ) ; 
                setVector (59 , 0.99815794 , -0.03796024 , -0.14718074, XYZ) ; 
                setVector (60 , -0.33399858 , -0.70328565 , -0.64285043, XYZ) ; 
                setVector (61 , 0.13320384 , -0.04097755 , 1.00000000, XYZ) ; 
                
                
        end
        
    case 63  
        switch (mdr_lBerhens_flag) 
            case 1 
                setVector (0 , 0.63180206 , 0.04702343 , -0.77962996, XYZ) ; 
                setVector (1 , -0.96970569 , 0.21965191 , -0.14363789, XYZ) ; 
                setVector (2 , -0.87518846 , -0.44988812 , -0.20212314, XYZ) ; 
                setVector (3 , 0.56148168 , -0.76377527 , -0.33255616, XYZ) ; 
                setVector (4 , -0.33916401 , 0.65961872 , 0.67755357, XYZ) ; 
                setVector (5 , 0.39735926 , -0.10358838 , -0.91683313, XYZ) ; 
                setVector (6 , 0.78889796 , 0.04168436 , 0.62057260, XYZ) ; 
                setVector (7 , -0.87463801 , -0.09346983 , -0.48526247, XYZ) ; 
                setVector (8 , 0.36850141 , 0.64411982 , 0.67714422, XYZ) ; 
                setVector (9 , 0.31454323 , -0.37960604 , -0.87531122, XYZ) ; 
                setVector (10 , 0.04618650 , -0.91186358 , -0.41902212, XYZ) ; 
                setVector (11 , 0.97235467 , 0.17940615 , 0.17761686, XYZ) ; 
                setVector (12 , 0.01120389 , -0.50978273 , 0.86556569, XYZ) ; 
                setVector (13 , 0.55954637 , -0.33712870 , 0.76319068, XYZ) ; 
                setVector (14 , -0.01465823 , -0.10068692 , -0.99942744, XYZ) ; 
                setVector (15 , 0.47302424 , -0.55029823 , -0.69471423, XYZ) ; 
                setVector (16 , 0.30518013 , 0.13355733 , 0.94775284, XYZ) ; 
                setVector (17 , 0.69754538 , 0.69107585 , 0.21225581, XYZ) ; 
                setVector (18 , -0.21700394 , 0.97388872 , -0.11686736, XYZ) ; 
                setVector (19 , 0.10200966 , 0.01715218 , -0.99925360, XYZ) ; 
                setVector (20 , -0.89161958 , 0.26097957 , 0.38224615, XYZ) ; 
                setVector (21 , -0.44276728 , -0.02660419 , -0.90136412, XYZ) ; 
                setVector (22 , -0.45158451 , 0.46592251 , -0.76693909, XYZ) ; 
                setVector (23 , -0.35504820 , 0.69723469 , -0.63008927, XYZ) ; 
                setVector (24 , -0.33164979 , 0.77651108 , 0.54428567, XYZ) ; 
                setVector (25 , -0.53596230 , -0.53940323 , 0.65650323, XYZ) ; 
                setVector (26 , -0.83087927 , -0.46997447 , 0.31300411, XYZ) ; 
                setVector (27 , -0.21923411 , -0.13229047 , 0.97141318, XYZ) ; 
                setVector (28 , -0.20414320 , 0.85876278 , 0.47964567, XYZ) ; 
                setVector (29 , -0.35123224 , 0.94046479 , -0.03701111, XYZ) ; 
                setVector (30 , -0.17629931 , 0.93311594 , 0.32775162, XYZ) ; 
                setVector (31 , -0.73694320 , 0.52977323 , -0.43065413, XYZ) ; 
                setVector (32 , -0.78442113 , -0.31157116 , 0.54480715, XYZ) ; 
                setVector (33 , -0.57551323 , -0.10175011 , -0.81709202, XYZ) ; 
                setVector (34 , 0.61934223 , -0.57042979 , 0.54793520, XYZ) ; 
                setVector (35 , -0.31077336 , -0.75381535 , 0.58684774, XYZ) ; 
                setVector (36 , -0.30745787 , 0.02791585 , -0.95598028, XYZ) ; 
                setVector (37 , 0.78217189 , 0.55935794 , 0.29074690, XYZ) ; 
                setVector (38 , -0.46782591 , -0.88048253 , 0.12287133, XYZ) ; 
                setVector (39 , -0.25130912 , 0.04621240 , 0.97155343, XYZ) ; 
                setVector (40 , 0.59557968 , 0.80858492 , -0.02613816, XYZ) ; 
                setVector (41 , -0.87147384 , -0.25373559 , -0.43053399, XYZ) ; 
                setVector (42 , -0.18373698 , 0.48373840 , 0.86107248, XYZ) ; 
                setVector (43 , 0.83988000 , 0.41206129 , 0.36608061, XYZ) ; 
                setVector (44 , 0.04881442 , 1.00000000 , -0.08261404, XYZ) ; 
                setVector (45 , 0.54530933 , -0.13521599 , -0.83280387, XYZ) ; 
                setVector (46 , -0.11215521 , -0.51140539 , -0.85737604, XYZ) ; 
                setVector (47 , -0.92804454 , 0.02712110 , -0.38367396, XYZ) ; 
                setVector (48 , -0.92613632 , 0.36976637 , 0.12145895, XYZ) ; 
                setVector (49 , 0.24422431 , -0.97175230 , -0.07252509, XYZ) ; 
                setVector (50 , 0.98059540 , -0.21045316 , -0.05787959, XYZ) ; 
                setVector (51 , 0.44715250 , 0.87445989 , 0.21114563, XYZ) ; 
                setVector (52 , -0.88734207 , -0.47000845 , -0.03039787, XYZ) ; 
                setVector (53 , -0.68196989 , 0.01083880 , -0.73756865, XYZ) ; 
                setVector (54 , 0.99973042 , 0.05247791 , -0.08362471, XYZ) ; 
                setVector (55 , -0.66615239 , -0.65520718 , -0.36898844, XYZ) ; 
                setVector (56 , 0.06717364 , -0.93406746 , 0.36361190, XYZ) ; 
                setVector (57 , -0.71038993 , 0.54605033 , 0.45429407, XYZ) ; 
                setVector (58 , 0.27620830 , 0.54996228 , -0.79401410, XYZ) ; 
                setVector (59 , 0.21546691 , -0.79232003 , 0.57880126, XYZ) ; 
                setVector (60 , 0.49822966 , -0.83924234 , 0.23800723, XYZ) ; 
                setVector (61 , -0.85392226 , 0.07887933 , -0.52326165, XYZ) ; 
                setVector (62 , -0.25864433 , 0.29310891 , -0.92541785, XYZ) ; 
                
            case 2 
                setVector (0 , -0.08450892 , -0.90393574 , 0.43007715, XYZ) ; 
                setVector (1 , 0.15945125 , -0.33122778 , -0.93491785, XYZ) ; 
                setVector (2 , -0.86087242 , -0.32831427 , 0.40039524, XYZ) ; 
                setVector (3 , -0.72722530 , 0.68052271 , -0.13130169, XYZ) ; 
                setVector (4 , -0.55592114 , -0.40825268 , -0.73040356, XYZ) ; 
                setVector (5 , -0.37607233 , 0.16305484 , -0.91716446, XYZ) ; 
                setVector (6 , -0.22798026 , -0.68542687 , 0.69815681, XYZ) ; 
                setVector (7 , -0.68693882 , -0.31378955 , -0.66246442, XYZ) ; 
                setVector (8 , -0.99999601 , -0.06058949 , -0.07446357, XYZ) ; 
                setVector (9 , -0.48660623 , 0.87584131 , -0.07296775, XYZ) ; 
                setVector (10 , 0.70192792 , -0.28152284 , -0.66124883, XYZ) ; 
                setVector (11 , 0.29563627 , 0.91934077 , 0.27680258, XYZ) ; 
                setVector (12 , 0.48084693 , 0.67134140 , 0.57209691, XYZ) ; 
                setVector (13 , 0.98066407 , -0.00194375 , 0.21794986, XYZ) ; 
                setVector (14 , -0.98880966 , 0.08729485 , 0.15441177, XYZ) ; 
                setVector (15 , -0.06857654 , -0.77036534 , -0.64112591, XYZ) ; 
                setVector (16 , -0.16286219 , -0.03597366 , -0.99065116, XYZ) ; 
                setVector (17 , 0.19181588 , 0.30201422 , -0.93872360, XYZ) ; 
                setVector (18 , -0.69507062 , -0.15535628 , -0.70848372, XYZ) ; 
                setVector (19 , -0.13150398 , -0.93146263 , -0.35255072, XYZ) ; 
                setVector (20 , 0.94283019 , 0.33051460 , 0.10506792, XYZ) ; 
                setVector (21 , -0.66201565 , -0.46185844 , -0.59802172, XYZ) ; 
                setVector (22 , 0.87760820 , -0.12154070 , -0.47353948, XYZ) ; 
                setVector (23 , 0.80116242 , -0.25755698 , 0.54864478, XYZ) ; 
                setVector (24 , 0.74882669 , -0.64758362 , -0.17059303, XYZ) ; 
                setVector (25 , 0.62829342 , -0.66558697 , 0.41406436, XYZ) ; 
                setVector (26 , 0.78523749 , -0.13868472 , -0.61104547, XYZ) ; 
                setVector (27 , -0.21551638 , 0.93834976 , -0.28681063, XYZ) ; 
                setVector (28 , -0.08849023 , 0.25798742 , -0.96686085, XYZ) ; 
                setVector (29 , 0.80568635 , 0.58668959 , 0.12598714, XYZ) ; 
                setVector (30 , -0.06448169 , -0.96797240 , 0.26092043, XYZ) ; 
                setVector (31 , -0.36790768 , -0.77069687 , -0.52903516, XYZ) ; 
                setVector (32 , -0.00351686 , 0.27034348 , 0.96752776, XYZ) ; 
                setVector (33 , 0.75533486 , 0.15950321 , -0.64283428, XYZ) ; 
                setVector (34 , 0.20729234 , 0.92927175 , -0.32045565, XYZ) ; 
                setVector (35 , 0.46974483 , -0.67423100 , -0.57789296, XYZ) ; 
                setVector (36 , 0.34354689 , 0.23720419 , -0.91373827, XYZ) ; 
                setVector (37 , -0.71441818 , 0.42037226 , 0.56754009, XYZ) ; 
                setVector (38 , -0.59755336 , 0.45941235 , -0.66413719, XYZ) ; 
                setVector (39 , 0.49098039 , 0.82015263 , -0.30902403, XYZ) ; 
                setVector (40 , -0.70547682 , -0.70559005 , 0.11684633, XYZ) ; 
                setVector (41 , -0.20202143 , 0.14488072 , -0.97334724, XYZ) ; 
                setVector (42 , 0.20152856 , -0.75174756 , -0.63519270, XYZ) ; 
                setVector (43 , 0.40582090 , -0.33122826 , 0.85720777, XYZ) ; 
                setVector (44 , 0.43135057 , 0.35540370 , 0.83476513, XYZ) ; 
                setVector (45 , 0.74016630 , 0.00663501 , -0.67920377, XYZ) ; 
                setVector (46 , 0.34951621 , 0.86771375 , -0.36622287, XYZ) ; 
                setVector (47 , -0.50196199 , 0.77798405 , -0.38984985, XYZ) ; 
                setVector (48 , -0.26300476 , -0.86996838 , -0.42800869, XYZ) ; 
                setVector (49 , 0.93972586 , 0.13090271 , 0.33013286, XYZ) ; 
                setVector (50 , 0.48243336 , -0.58679695 , 0.65737000, XYZ) ; 
                setVector (51 , -0.41666174 , -0.83497576 , -0.37204356, XYZ) ; 
                setVector (52 , -0.42364894 , -0.77919129 , 0.47179490, XYZ) ; 
                setVector (53 , -0.88932982 , 0.20222631 , -0.42119463, XYZ) ; 
                setVector (54 , -0.68349236 , 0.73546043 , 0.03382429, XYZ) ; 
                setVector (55 , 0.32927296 , 0.94786026 , -0.04845802, XYZ) ; 
                setVector (56 , -0.41443707 , 0.90872611 , 0.10801245, XYZ) ; 
                setVector (57 , 0.13536776 , -0.16389398 , -0.98184635, XYZ) ; 
                setVector (58 , 0.12645050 , -0.98122603 , -0.17439517, XYZ) ; 
                setVector (59 , 0.06751367 , -0.76793888 , 0.64414261, XYZ) ; 
                setVector (60 , 0.21692059 , 0.80208254 , 0.56463881, XYZ) ; 
                setVector (61 , 0.94906120 , -0.24058584 , -0.22496491, XYZ) ; 
                setVector (62 , 0.79415551 , -0.61514693 , -0.01091827, XYZ) ; 
                
            case 3 
                setVector (0 , 0.20807479 , -0.87728276 , 0.44304376, XYZ) ; 
                setVector (1 , -0.79396610 , -0.20847763 , -0.57910521, XYZ) ; 
                setVector (2 , 0.55158286 , -0.83639958 , -0.07348489, XYZ) ; 
                setVector (3 , -0.42792147 , -0.48220954 , 0.77043176, XYZ) ; 
                setVector (4 , 0.07643444 , 0.72232571 , -0.69398218, XYZ) ; 
                setVector (5 , -0.87183064 , -0.47738820 , 0.14567004, XYZ) ; 
                setVector (6 , 0.68691273 , 0.30290079 , -0.66754022, XYZ) ; 
                setVector (7 , -0.71945717 , -0.58766704 , 0.38240915, XYZ) ; 
                setVector (8 , -0.06775283 , -0.19698873 , 0.98275781, XYZ) ; 
                setVector (9 , 0.02753000 , 0.99952665 , 0.09693557, XYZ) ; 
                setVector (10 , -0.90401998 , 0.40375111 , -0.17012011, XYZ) ; 
                setVector (11 , -0.09860575 , 0.99972115 , -0.00651548, XYZ) ; 
                setVector (12 , -0.81075286 , 0.41166642 , 0.42710478, XYZ) ; 
                setVector (13 , -0.80843569 , 0.27873535 , 0.52720610, XYZ) ; 
                setVector (14 , 0.16087211 , 0.97301863 , 0.19121411, XYZ) ; 
                setVector (15 , -0.36268972 , 0.79410252 , -0.49705662, XYZ) ; 
                setVector (16 , 0.23777205 , 0.84945633 , -0.48072479, XYZ) ; 
                setVector (17 , -0.52501719 , -0.54612642 , -0.65978088, XYZ) ; 
                setVector (18 , -0.72354629 , 0.42385354 , -0.55320599, XYZ) ; 
                setVector (19 , -0.92722171 , -0.19194785 , 0.33559476, XYZ) ; 
                setVector (20 , -0.60062177 , 0.55328289 , 0.58509782, XYZ) ; 
                setVector (21 , -0.96997410 , -0.19998627 , 0.16841518, XYZ) ; 
                setVector (22 , -0.73915004 , 0.61625155 , -0.28826926, XYZ) ; 
                setVector (23 , 0.82562646 , -0.55153885 , 0.15281941, XYZ) ; 
                setVector (24 , -0.77811542 , -0.36450960 , -0.52045853, XYZ) ; 
                setVector (25 , -0.97120854 , -0.05164325 , 0.25158473, XYZ) ; 
                setVector (26 , -0.57762924 , 0.28398308 , 0.77130150, XYZ) ; 
                setVector (27 , -0.29514347 , -0.28990817 , -0.91545153, XYZ) ; 
                setVector (28 , -0.68763165 , 0.65668513 , 0.32424570, XYZ) ; 
                setVector (29 , -0.34670057 , -0.91788714 , 0.21561503, XYZ) ; 
                setVector (30 , 0.84872423 , -0.51612978 , -0.14995050, XYZ) ; 
                setVector (31 , 0.98317597 , 0.20633209 , 0.00013473, XYZ) ; 
                setVector (32 , 0.30180659 , -0.93131189 , -0.22534167, XYZ) ; 
                setVector (33 , -0.58820757 , -0.80342478 , -0.13314806, XYZ) ; 
                setVector (34 , 0.02517170 , -0.42834994 , -0.90834500, XYZ) ; 
                setVector (35 , 0.58995128 , -0.66693559 , -0.46514764, XYZ) ; 
                setVector (36 , 0.99915910 , -0.09215661 , 0.04895082, XYZ) ; 
                setVector (37 , 0.07777428 , -0.98297844 , 0.19212621, XYZ) ; 
                setVector (38 , -0.20373191 , 0.68302845 , -0.70793599, XYZ) ; 
                setVector (39 , -0.91681262 , -0.03754967 , 0.40896525, XYZ) ; 
                setVector (40 , 0.35932651 , -0.86827476 , 0.35523418, XYZ) ; 
                setVector (41 , 0.06422143 , -0.70627799 , -0.71151594, XYZ) ; 
                setVector (42 , 0.79444233 , -0.53365131 , -0.30542690, XYZ) ; 
                setVector (43 , 0.33607442 , -0.86911326 , -0.37537188, XYZ) ; 
                setVector (44 , 0.46537458 , 0.32403962 , -0.82923625, XYZ) ; 
                setVector (45 , 0.11598874 , 0.81434960 , -0.57670553, XYZ) ; 
                setVector (46 , 0.62555291 , -0.76260360 , -0.19059702, XYZ) ; 
                setVector (47 , -0.61189368 , 0.78987589 , -0.10435580, XYZ) ; 
                setVector (48 , -0.36612299 , -0.06988998 , 0.93288653, XYZ) ; 
                setVector (49 , 0.45241209 , -0.78060659 , -0.44179699, XYZ) ; 
                setVector (50 , 0.13568493 , 0.36288034 , 0.92688477, XYZ) ; 
                setVector (51 , 0.04808653 , -0.57368433 , -0.82327511, XYZ) ; 
                setVector (52 , -0.83811878 , -0.02322378 , 0.55337645, XYZ) ; 
                setVector (53 , -0.90965959 , 0.30791104 , -0.29481884, XYZ) ; 
                setVector (54 , 0.61557678 , 0.20709795 , -0.76640955, XYZ) ; 
                setVector (55 , -0.06905549 , 0.82540970 , 0.56845237, XYZ) ; 
                setVector (56 , -0.62665069 , 0.73894625 , -0.26547179, XYZ) ; 
                setVector (57 , -0.17985065 , -0.98812946 , -0.02149062, XYZ) ; 
                setVector (58 , -0.27559496 , -0.43882872 , -0.86063041, XYZ) ; 
                setVector (59 , 0.11508701 , 0.59764423 , -0.79923981, XYZ) ; 
                setVector (60 , -0.94038074 , -0.34592893 , 0.07228528, XYZ) ; 
                setVector (61 , -0.49138102 , -0.17003632 , 0.85955818, XYZ) ; 
                setVector (62 , -0.15800405 , -0.45841427 , 0.87982896, XYZ) ; 
                
            case 4 
                setVector (0 , -0.28509126 , 0.22177012 , 0.93741608, XYZ) ; 
                setVector (1 , -0.44477856 , -0.19633056 , -0.87910993, XYZ) ; 
                setVector (2 , 0.77804504 , 0.59664014 , -0.21880217, XYZ) ; 
                setVector (3 , 0.48627071 , 0.67065353 , -0.56830675, XYZ) ; 
                setVector (4 , -0.15214934 , -0.20393948 , -0.97183702, XYZ) ; 
                setVector (5 , 0.85694968 , 0.17568381 , -0.49394370, XYZ) ; 
                setVector (6 , -0.91431570 , -0.33955518 , 0.24070109, XYZ) ; 
                setVector (7 , 0.19983088 , -0.62402150 , -0.76149374, XYZ) ; 
                setVector (8 , 0.31677782 , 0.94616582 , 0.11674751, XYZ) ; 
                setVector (9 , -0.56991675 , -0.71283620 , 0.41984208, XYZ) ; 
                setVector (10 , -0.08548018 , -0.86816194 , -0.49819264, XYZ) ; 
                setVector (11 , -0.24862411 , -0.57752999 , -0.78348778, XYZ) ; 
                setVector (12 , -0.73845525 , -0.51328372 , -0.44769588, XYZ) ; 
                setVector (13 , -0.94746278 , 0.33196681 , -0.03633515, XYZ) ; 
                setVector (14 , 0.49868719 , -0.69217086 , 0.53048894, XYZ) ; 
                setVector (15 , -0.75866533 , 0.13171734 , -0.64520182, XYZ) ; 
                setVector (16 , -0.61772996 , -0.60006339 , -0.51724419, XYZ) ; 
                setVector (17 , -0.57332814 , -0.25926066 , -0.78312622, XYZ) ; 
                setVector (18 , -0.91868956 , -0.29926530 , -0.27505946, XYZ) ; 
                setVector (19 , 0.62848292 , 0.58917108 , -0.51681194, XYZ) ; 
                setVector (20 , 0.03266019 , 0.36331733 , -0.93602444, XYZ) ; 
                setVector (21 , 0.71157332 , 0.70738405 , 0.04979086, XYZ) ; 
                setVector (22 , -0.88228762 , 0.48027530 , -0.01058871, XYZ) ; 
                setVector (23 , 0.50679991 , 0.01270033 , -0.86729492, XYZ) ; 
                setVector (24 , 0.59641438 , -0.42435750 , -0.68805416, XYZ) ; 
                setVector (25 , -0.04952299 , 0.64533529 , -0.76830838, XYZ) ; 
                setVector (26 , -0.66805373 , -0.45146401 , 0.59924318, XYZ) ; 
                setVector (27 , -0.35890493 , 0.91675006 , -0.19991126, XYZ) ; 
                setVector (28 , 0.50598695 , -0.20623497 , 0.84300194, XYZ) ; 
                setVector (29 , -0.55377325 , 0.06011570 , -0.83601987, XYZ) ; 
                setVector (30 , 0.64322548 , -0.17910130 , 0.75059419, XYZ) ; 
                setVector (31 , 0.82792251 , -0.37078755 , 0.43158874, XYZ) ; 
                setVector (32 , -0.46283370 , -0.89042288 , -0.04626004, XYZ) ; 
                setVector (33 , 0.19770237 , 0.97309063 , -0.15236902, XYZ) ; 
                setVector (34 , 0.46190352 , -0.85564551 , -0.25243579, XYZ) ; 
                setVector (35 , 0.58487093 , 0.79477447 , -0.18832811, XYZ) ; 
                setVector (36 , -0.00036841 , -0.96819872 , -0.26795340, XYZ) ; 
                setVector (37 , 0.33181175 , -0.58108080 , 0.74930234, XYZ) ; 
                setVector (38 , -0.43601322 , 0.26264790 , 0.86609265, XYZ) ; 
                setVector (39 , 0.46151267 , -0.41197412 , -0.79151204, XYZ) ; 
                setVector (40 , 0.65953055 , 0.70069243 , -0.28854375, XYZ) ; 
                setVector (41 , 0.08807534 , 0.64810657 , 0.76250150, XYZ) ; 
                setVector (42 , -0.29439876 , 0.44813864 , -0.84953461, XYZ) ; 
                setVector (43 , 0.95458721 , -0.14633253 , 0.27669114, XYZ) ; 
                setVector (44 , 0.13706668 , -0.40870427 , 0.90740370, XYZ) ; 
                setVector (45 , 0.76230864 , 0.45756328 , -0.46768505, XYZ) ; 
                setVector (46 , -0.33412669 , 0.52564601 , 0.78820274, XYZ) ; 
                setVector (47 , 0.17433064 , -0.55377691 , 0.81984626, XYZ) ; 
                setVector (48 , -0.69134922 , 0.30648167 , -0.66129658, XYZ) ; 
                setVector (49 , -0.94622796 , 0.10435511 , 0.32089031, XYZ) ; 
                setVector (50 , 0.40532807 , 0.50467346 , 0.76825893, XYZ) ; 
                setVector (51 , 0.88178415 , -0.39338294 , -0.27733465, XYZ) ; 
                setVector (52 , -0.80745276 , -0.59604331 , 0.04427584, XYZ) ; 
                setVector (53 , -0.66982944 , 0.13081330 , 0.73717320, XYZ) ; 
                setVector (54 , -0.56637012 , -0.77794904 , -0.28849282, XYZ) ; 
                setVector (55 , 0.83249143 , -0.47284916 , 0.30426898, XYZ) ; 
                setVector (56 , -0.22336504 , -0.70335758 , -0.68161873, XYZ) ; 
                setVector (57 , 0.31558352 , 0.39852808 , -0.86648158, XYZ) ; 
                setVector (58 , -0.56856576 , -0.38941955 , 0.73095370, XYZ) ; 
                setVector (59 , 0.38432224 , 0.62565857 , -0.68560608, XYZ) ; 
                setVector (60 , -0.52900422 , -0.73413789 , -0.43635309, XYZ) ; 
                setVector (61 , 0.04904060 , -0.08248976 , 0.99999919, XYZ) ; 
                setVector (62 , -0.04829269 , 0.85989217 , -0.51716650, XYZ) ; 
                
        end                
        
        
    case 64 
        
        switch (mdr_lBerhens_flag) 
            case 1 
                setVector (0 , 0.53230019 , -0.84537418 , 0.00000000, XYZ);
                setVector (1 , 0.08949872 , -0.90755578 , 0.40784004, XYZ);
                setVector (2 , 0.98472637 , 0.16827055 , 0.00000000, XYZ);
                setVector (3 , 0.49178125 , -0.78102390 , -0.38230075, XYZ);
                setVector (4 , 0.24116334 , -0.77921387 , 0.57677292, XYZ);
                setVector (5 , 0.03882314 , -0.45322020 , 0.88942971, XYZ);
                setVector (6 , 0.90788893 , -0.35940194 , 0.21111355, XYZ);
                setVector (7 , 0.68699651 , -0.72528394 , 0.00000000, XYZ);
                setVector (8 , 0.45322020 , 0.03882314 , 0.88942971, XYZ);
                setVector (9 , 0.37558875 , 0.48597224 , 0.78788645, XYZ);
                setVector (10 , 0.72528394 , 0.68699651 , 0.00000000, XYZ);
                setVector (11 , 0.10384662 , -0.16492432 , -0.97980449, XYZ);
                setVector (12 , 0.76767670 , -0.02000842 , -0.63896255, XYZ);
                setVector (13 , 0.65047105 , -0.64148449 , -0.40421042, XYZ);
                setVector (14 , 0.96803024 , 0.13186641 , -0.20863775, XYZ);
                setVector (15 , 0.35940194 , 0.90788893 , -0.21111355, XYZ);
                setVector (16 , 0.03512096 , -0.84667608 , 0.52906252, XYZ);
                setVector (17 , 0.25471318 , -0.60433774 , 0.75359014, XYZ);
                setVector (18 , 0.15092408 , -0.98194503 , 0.10491364, XYZ);
                setVector (19 , 0.05286024 , -0.99203372 , 0.10524205, XYZ);
                setVector (20 , 0.52973701 , -0.84130347 , 0.09791912, XYZ);
                setVector (21 , 0.02091525 , -0.22047740 , 0.97414232, XYZ);
                setVector (22 , 0.74834412 , -0.39758640 , 0.52906252, XYZ);
                setVector (23 , 0.52872604 , -0.25387224 , 0.80870183, XYZ);
                setVector (24 , 0.94914598 , -0.29289510 , 0.10646762, XYZ);
                setVector (25 , 0.60433774 , 0.25471318 , 0.75359014, XYZ);
                setVector (26 , 0.57438237 , -0.71530934 , 0.39549770, XYZ);
                setVector (27 , 0.55136780 , -0.67814234 , 0.48385691, XYZ);
                setVector (28 , 0.89576701 , 0.32199886 , 0.30318179, XYZ);
                setVector (29 , 0.79029522 , 0.37325185 , 0.48385691, XYZ);
                setVector (30 , 0.99647804 , 0.07094022 , 0.00000000, XYZ);
                setVector (31 , 0.84130347 , 0.52973701 , 0.09791912, XYZ);
                setVector (32 , 0.39057968 , 0.12195417 , 0.91135926, XYZ);
                setVector (33 , 0.39850663 , 0.25092457 , 0.88103934, XYZ);
                setVector (34 , 0.16879818 , 0.74279224 , 0.64634964, XYZ);
                setVector (35 , 0.43796269 , 0.65565191 , 0.61344132, XYZ);
                setVector (36 , 0.21626161 , 0.95165386 , 0.21351075, XYZ);
                setVector (37 , 0.56047656 , 0.82028053 , 0.10491364, XYZ);
                setVector (38 , 0.49462890 , 0.86795349 , 0.00000000, XYZ);
                setVector (39 , 0.65445407 , 0.75477870 , 0.00000000, XYZ);
                setVector (40 , 0.27877677 , 0.29951430 , 0.91135926, XYZ);
                setVector (41 , 0.03512096 , -0.84667608 , -0.52906252, XYZ);
                setVector (42 , 0.31144264 , -0.69443281 , -0.64712252, XYZ);
                setVector (43 , 0.02091525 , -0.22047740 , -0.97414232, XYZ);
                setVector (44 , 0.11131360 , -0.94313975 , -0.30999629, XYZ);
                setVector (45 , 0.34390222 , -0.74544747 , -0.56924540, XYZ);
                setVector (46 , 0.80019165 , -0.25995723 , -0.53862471, XYZ);
                setVector (47 , 0.67224833 , -0.07470960 , -0.73518818, XYZ);
                setVector (48 , 0.65565191 , -0.43796269 , -0.61344132, XYZ);
                setVector (49 , 0.78053707 , 0.11159136 , -0.61344132, XYZ);
                setVector (50 , 0.75187377 , -0.62509101 , -0.20481223, XYZ);
                setVector (51 , 0.62762561 , -0.60088612 , -0.49297359, XYZ);
                setVector (52 , 0.42668316 , -0.09696297 , -0.89807609, XYZ);
                setVector (53 , 0.24102090 , -0.19469257 , -0.94973929, XYZ);
                setVector (54 , 0.39316608 , -0.42512983 , -0.81405532, XYZ);
                setVector (55 , 0.89576701 , 0.32199886 , -0.30318179, XYZ);
                setVector (56 , 0.97595326 , 0.03196157 , -0.21093763, XYZ);
                setVector (57 , 0.92968916 , 0.35174319 , -0.09977873, XYZ);
                setVector (58 , 0.34022926 , 0.83943603 , -0.42141690, XYZ);
                setVector (59 , 0.14607057 , 0.98251159 , -0.10646762, XYZ);
                setVector (60 , 0.49087882 , 0.43490331 , -0.75359014, XYZ);
                setVector (61 , 0.56047656 , 0.82028053 , -0.10491364, XYZ);
                setVector (62 , 0.71898906 , 0.68636720 , -0.09977873, XYZ);
                setVector (63 , 0.27877677 , 0.29951430 , -0.91135926, XYZ);
                
            case 2 
                setVector (0 , 0.84537418 , 0.53230019 , 0.00000000, XYZ);
                setVector (1 , 0.16827055 , -0.98472637 , 0.00000000, XYZ);
                setVector (2 , 0.78102390 , 0.49178125 , 0.38230075, XYZ);
                setVector (3 , 0.56243320 , -0.12781192 , -0.81568008, XYZ);
                setVector (4 , 0.12867507 , -0.60056527 , 0.78788645, XYZ);
                setVector (5 , 0.29573014 , -0.46966474 , 0.83063814, XYZ);
                setVector (6 , 0.97404926 , -0.06826509 , 0.21111355, XYZ);
                setVector (7 , 0.96803024 , 0.13186641 , 0.20863775, XYZ);
                setVector (8 , 0.39193019 , -0.23088203 , 0.88942971, XYZ);
                setVector (9 , 0.55423733 , 0.59846051 , 0.57677292, XYZ);
                setVector (10 , 0.02707330 , 0.99863309 , 0.00000000, XYZ);
                setVector (11 , 0.29736400 , -0.86382265 , -0.40421042, XYZ);
                setVector (12 , 0.97404926 , -0.06826509 , -0.21111355, XYZ);
                setVector (13 , 0.67886335 , -0.70549291 , -0.19855809, XYZ);
                setVector (14 , 0.91708953 , 0.34281552 , -0.19855809, XYZ);
                setVector (15 , 0.55423733 , 0.59846051 , -0.57677292, XYZ);
                setVector (16 , 0.01070843 , -0.94892066 , 0.31214725, XYZ);
                setVector (17 , 0.21670153 , -0.92483040 , 0.30940292, XYZ);
                setVector (18 , 0.25796033 , -0.95990301 , 0.10021817, XYZ);
                setVector (19 , 0.07094022 , -0.99647804 , 0.00000000, XYZ);
                setVector (20 , 0.12195417 , -0.39057968 , 0.91135926, XYZ);
                setVector (21 , 0.15451859 , -0.24539918 , 0.95598340, XYZ);
                setVector (22 , 0.83943603 , -0.34022926 , 0.42141690, XYZ);
                setVector (23 , 0.43490331 , -0.49087882 , 0.75359014, XYZ);
                setVector (24 , 0.98251159 , -0.14607057 , 0.10646762, XYZ);
                setVector (25 , 0.74034766 , -0.59511022 , 0.30940292, XYZ);
                setVector (26 , 0.59085912 , -0.75028922 , 0.29317672, XYZ);
                setVector (27 , 0.52360889 , -0.63229297 , 0.56924540, XYZ);
                setVector (28 , 0.94818156 , 0.23875652 , 0.20481223, XYZ);
                setVector (29 , 0.82567519 , 0.27062615 , 0.49297359, XYZ);
                setVector (30 , 0.99117075 , -0.12482610 , 0.00000000, XYZ);
                setVector (31 , 0.42668316 , -0.09696297 , 0.89807609, XYZ);
                setVector (32 , 0.30149168 , 0.07140758 , 0.94973929, XYZ);
                setVector (33 , 0.53629970 , 0.33768766 , 0.77223744, XYZ);
                setVector (34 , 0.11139878 , 0.83395128 , 0.53862471, XYZ);
                setVector (35 , 0.47019216 , 0.54837328 , 0.69007759, XYZ);
                setVector (36 , 0.41975931 , 0.85109765 , 0.31214725, XYZ);
                setVector (37 , 0.50798552 , 0.80240514 , 0.30999629, XYZ);
                setVector (38 , 0.31579542 , 0.94777332 , 0.00000000, XYZ);
                setVector (39 , 0.78912893 , 0.61259818 , 0.00000000, XYZ);
                setVector (40 , 0.18036951 , 0.47990231 , 0.85741567, XYZ);
                setVector (41 , 0.01070843 , -0.94892066 , -0.31214725, XYZ);
                setVector (42 , 0.07140758 , -0.30149168 , -0.94973929, XYZ);
                setVector (43 , 0.15451859 , -0.24539918 , -0.95598340, XYZ);
                setVector (44 , 0.19453419 , -0.89073889 , -0.40834015, XYZ);
                setVector (45 , 0.46944685 , -0.74555347 , -0.47092534, XYZ);
                setVector (46 , 0.83395128 , -0.11139878 , -0.53862471, XYZ);
                setVector (47 , 0.93611167 , -0.13548493 , -0.32146504, XYZ);
                setVector (48 , 0.59889207 , -0.34781937 , -0.71996595, XYZ);
                setVector (49 , 0.69772332 , 0.18701897 , -0.69007759, XYZ);
                setVector (50 , 0.59085912 , -0.75028922 , -0.29317672, XYZ);
                setVector (51 , 0.52360889 , -0.63229297 , -0.56924540, XYZ);
                setVector (52 , 0.36568934 , -0.01177102 , -0.92958794, XYZ);
                setVector (53 , 0.46552603 , 0.16875696 , -0.86764486, XYZ);
                setVector (54 , 0.18978485 , -0.11414660 , -0.97414232, XYZ);
                setVector (55 , 0.95990301 , 0.25796033 , -0.10021817, XYZ);
                setVector (56 , 0.99218297 , -0.04795153 , -0.10618195, XYZ);
                setVector (57 , 0.88955872 , 0.44385360 , -0.09838828, XYZ);
                setVector (58 , 0.39758640 , 0.74834412 , -0.52906252, XYZ);
                setVector (59 , 0.29289510 , 0.94914598 , -0.10646762, XYZ);
                setVector (60 , 0.59511022 , 0.74034766 , -0.30940292, XYZ);
                setVector (61 , 0.64747453 , 0.75414458 , -0.10021817, XYZ);
                setVector (62 , 0.77117370 , 0.60349000 , -0.19771682, XYZ);
                setVector (63 , 0.19469257 , 0.24102090 , -0.94973929, XYZ);
                
            case 3 
                setVector (0 , 0.37639307 , -0.59776982 , 0.70639967, XYZ);
                setVector (1 , 0.49178125 , -0.78102390 , 0.38230075, XYZ);
                setVector (2 , 0.32351070 , 0.20370246 , 0.92295565, XYZ);
                setVector (3 , 0.77974385 , -0.47293449 , -0.40784004, XYZ);
                setVector (4 , 0.13186641 , -0.96803024 , 0.20863775, XYZ);
                setVector (5 , 0.10384662 , -0.16492432 , 0.97980449, XYZ);
                setVector (6 , 0.77921387 , 0.24116334 , 0.57677292, XYZ);
                setVector (7 , 0.86382265 , 0.29736400 , 0.40421042, XYZ);
                setVector (8 , 0.27382222 , -0.06222560 , 0.95872330, XYZ);
                setVector (9 , 0.35940194 , 0.90788893 , 0.21111355, XYZ);
                setVector (10 , 0.23088203 , 0.39193019 , 0.88942971, XYZ);
                setVector (11 , 0.13186641 , -0.96803024 , -0.20863775, XYZ);
                setVector (12 , 0.90788893 , -0.35940194 , -0.21111355, XYZ);
                setVector (13 , 0.39193019 , -0.23088203 , -0.88942971, XYZ);
                setVector (14 , 0.70290294 , 0.44259143 , -0.55501466, XYZ);
                setVector (15 , 0.37558875 , 0.48597224 , -0.78788645, XYZ);
                setVector (16 , 0.11159136 , -0.78053707 , 0.61344132, XYZ);
                setVector (17 , 0.32199886 , -0.89576701 , 0.30318179, XYZ);
                setVector (18 , 0.37325185 , -0.79029522 , 0.48385691, XYZ);
                setVector (19 , 0.35174319 , -0.92968916 , 0.09977873, XYZ);
                setVector (20 , 0.07140758 , -0.30149168 , 0.94973929, XYZ);
                setVector (21 , 0.05217454 , -0.08286116 , 0.99418955, XYZ);
                setVector (22 , 0.90403321 , -0.05597096 , 0.42141690, XYZ);
                setVector (23 , 0.49163278 , -0.58097388 , 0.64712252, XYZ);
                setVector (24 , 0.78053707 , 0.11159136 , 0.61344132, XYZ);
                setVector (25 , 0.75187377 , -0.62509101 , 0.20481223, XYZ);
                setVector (26 , 0.87175942 , -0.47640384 , 0.10524205, XYZ);
                setVector (27 , 0.60349000 , -0.77117370 , 0.19771682, XYZ);
                setVector (28 , 0.85708832 , 0.42124579 , 0.29317672, XYZ);
                setVector (29 , 0.74544747 , 0.34390222 , 0.56924540, XYZ);
                setVector (30 , 0.92968916 , 0.35174319 , 0.09977873, XYZ);
                setVector (31 , 0.33486156 , -0.14742790 , 0.92958794, XYZ);
                setVector (32 , 0.34687291 , -0.35337315 , 0.86764486, XYZ);
                setVector (33 , 0.22047740 , 0.02091525 , 0.97414232, XYZ);
                setVector (34 , 0.07470960 , 0.67224833 , 0.73518818, XYZ);
                setVector (35 , 0.51801762 , 0.69504842 , 0.49655458, XYZ);
                setVector (36 , 0.29289510 , 0.94914598 , 0.10646762, XYZ);
                setVector (37 , 0.56038637 , 0.71918457 , 0.40834015, XYZ);
                setVector (38 , 0.60088612 , 0.62762561 , 0.49297359, XYZ);
                setVector (39 , 0.04795153 , 0.99218297 , 0.10618195, XYZ);
                setVector (40 , 0.04474807 , 0.51072203 , 0.85741567, XYZ);
                setVector (41 , 0.05483151 , -0.69039376 , -0.71996595, XYZ);
                setVector (42 , 0.12195417 , -0.39057968 , -0.91135926, XYZ);
                setVector (43 , 0.05217454 , -0.08286116 , -0.99418955, XYZ);
                setVector (44 , 0.39682224 , -0.82711225 , -0.39549770, XYZ);
                setVector (45 , 0.41147360 , -0.65348308 , -0.63375889, XYZ);
                setVector (46 , 0.74279224 , -0.16879818 , -0.64634964, XYZ);
                setVector (47 , 0.90272820 , -0.28238809 , -0.32146504, XYZ);
                setVector (48 , 0.54837328 , -0.47019216 , -0.69007759, XYZ);
                setVector (49 , 0.85066236 , 0.16675767 , -0.49655458, XYZ);
                setVector (50 , 0.57438237 , -0.71530934 , -0.39549770, XYZ);
                setVector (51 , 0.86629019 , -0.45060805 , -0.21093763, XYZ);
                setVector (52 , 0.33486156 , -0.14742790 , -0.92958794, XYZ);
                setVector (53 , 0.52955884 , 0.08417288 , -0.84290174, XYZ);
                setVector (54 , 0.22047740 , 0.02091525 , -0.97414232, XYZ);
                setVector (55 , 0.98194503 , 0.15092408 , -0.10491364, XYZ);
                setVector (56 , 0.82567519 , 0.27062615 , -0.49297359, XYZ);
                setVector (57 , 0.80897266 , 0.50937952 , -0.28999440, XYZ);
                setVector (58 , 0.22315866 , 0.63851356 , -0.73518818, XYZ);
                setVector (59 , 0.43796269 , 0.65565191 , -0.61344132, XYZ);
                setVector (60 , 0.67750557 , 0.66863145 , -0.30318179, XYZ);
                setVector (61 , 0.67814234 , 0.55136780 , -0.48385691, XYZ);
                setVector (62 , 0.78470967 , 0.61036995 , -0.09838828, XYZ);
                setVector (63 , 0.30476978 , 0.44117234 , -0.84290174, XYZ);
                
            case 4 
                setVector (0 , 0.59776982 , 0.37639307 , 0.70639967, XYZ);
                setVector (1 , 0.20370246 , -0.32351070 , 0.92295565, XYZ);
                setVector (2 , 0.47293449 , 0.77974385 , 0.40784004, XYZ);
                setVector (3 , 0.90755578 , 0.08949872 , -0.40784004, XYZ);
                setVector (4 , 0.29736400 , -0.86382265 , 0.40421042, XYZ);
                setVector (5 , 0.88116348 , -0.20024280 , 0.42597503, XYZ);
                setVector (6 , 0.60056527 , 0.12867507 , 0.78788645, XYZ);
                setVector (7 , 0.91708953 , 0.34281552 , 0.19855809, XYZ);
                setVector (8 , 0.46966474 , 0.29573014 , 0.83063814, XYZ);
                setVector (9 , 0.06826509 , 0.97404926 , 0.21111355, XYZ);
                setVector (10 , 0.06222560 , 0.27382222 , 0.95872330, XYZ);
                setVector (11 , 0.34281552 , -0.91708953 , -0.19855809, XYZ);
                setVector (12 , 0.59846051 , -0.55423733 , -0.57677292, XYZ);
                setVector (13 , 0.45322020 , 0.03882314 , -0.88942971, XYZ);
                setVector (14 , 0.82913056 , 0.52207219 , -0.19489523, XYZ);
                setVector (15 , 0.53727691 , 0.81596863 , -0.20863775, XYZ);
                setVector (16 , 0.05483151 , -0.69039376 , 0.71996595, XYZ);
                setVector (17 , 0.23875652 , -0.94818156 , 0.20481223, XYZ);
                setVector (18 , 0.27062615 , -0.82567519 , 0.49297359, XYZ);
                setVector (19 , 0.43462039 , -0.87750452 , 0.19771682, XYZ);
                setVector (20 , 0.08417288 , -0.52955884 , 0.84290174, XYZ);
                setVector (21 , 0.80019165 , -0.25995723 , 0.53862471, XYZ);
                setVector (22 , 0.84667608 , 0.03512096 , 0.52906252, XYZ);
                setVector (23 , 0.90272820 , -0.28238809 , 0.32146504, XYZ);
                setVector (24 , 0.69039376 , 0.05483151 , 0.71996595, XYZ);
                setVector (25 , 0.66863145 , -0.67750557 , 0.30318179, XYZ);
                setVector (26 , 0.86629019 , -0.45060805 , 0.21093763, XYZ);
                setVector (27 , 0.68636720 , -0.71898906 , 0.09977873, XYZ);
                setVector (28 , 0.82711225 , 0.39682224 , 0.39549770, XYZ);
                setVector (29 , 0.74555347 , 0.46944685 , 0.47092534, XYZ);
                setVector (30 , 0.87750452 , 0.43462039 , 0.19771682, XYZ);
                setVector (31 , 0.36568934 , -0.01177102 , 0.92958794, XYZ);
                setVector (32 , 0.44117234 , -0.30476978 , 0.84290174, XYZ);
                setVector (33 , 0.18978485 , -0.11414660 , 0.97414232, XYZ);
                setVector (34 , 0.22315866 , 0.63851356 , 0.73518818, XYZ);
                setVector (35 , 0.25387224 , 0.52872604 , 0.80870183, XYZ);
                setVector (36 , 0.14607057 , 0.98251159 , 0.10646762, XYZ);
                setVector (37 , 0.71530934 , 0.57438237 , 0.39549770, XYZ);
                setVector (38 , 0.67814234 , 0.55136780 , 0.48385691, XYZ);
                setVector (39 , 0.12482610 , 0.99117075 , 0.00000000, XYZ);
                setVector (40 , 0.35337315 , 0.34687291 , 0.86764486, XYZ);
                setVector (41 , 0.11159136 , -0.78053707 , -0.61344132, XYZ);
                setVector (42 , 0.16875696 , -0.46552603 , -0.86764486, XYZ);
                setVector (43 , 0.21670153 , -0.92483040 , -0.30940292, XYZ);
                setVector (44 , 0.42124579 , -0.85708832 , -0.29317672, XYZ);
                setVector (45 , 0.43462039 , -0.87750452 , -0.19771682, XYZ);
                setVector (46 , 0.84667608 , 0.03512096 , -0.52906252, XYZ);
                setVector (47 , 0.95165386 , -0.21626161 , -0.21351075, XYZ);
                setVector (48 , 0.52872604 , -0.25387224 , -0.80870183, XYZ);
                setVector (49 , 0.58651698 , 0.00043530 , -0.80870183, XYZ);
                setVector (50 , 0.71918457 , -0.56038637 , -0.40834015, XYZ);
                setVector (51 , 0.87175942 , -0.47640384 , -0.10524205, XYZ);
                setVector (52 , 0.30149168 , 0.07140758 , -0.94973929, XYZ);
                setVector (53 , 0.53828433 , 0.21345942 , -0.81405532, XYZ);
                setVector (54 , 0.13568270 , -0.03083365 , -0.98926260, XYZ);
                setVector (55 , 0.94313975 , 0.11131360 , -0.30999629, XYZ);
                setVector (56 , 0.79029522 , 0.37325185 , -0.48385691, XYZ);
                setVector (57 , 0.84130347 , 0.52973701 , -0.09791912, XYZ);
                setVector (58 , 0.07470960 , 0.67224833 , -0.73518818, XYZ);
                setVector (59 , 0.34781937 , 0.59889207 , -0.71996595, XYZ);
                setVector (60 , 0.62509101 , 0.75187377 , -0.20481223, XYZ);
                setVector (61 , 0.60088612 , 0.62762561 , -0.49297359, XYZ);
                setVector (62 , 0.09696297 , 0.42668316 , -0.89807609, XYZ);
                setVector (63 , 0.35337315 , 0.34687291 , -0.86764486, XYZ);
                
            case 5 
                setVector (0 , 0.97416289 , -0.22137674 , 0.00000000, XYZ);
                setVector (1 , 0.90755578 , 0.08949872 , 0.40784004, XYZ);
                setVector (2 , 0.12781192 , 0.56243320 , 0.81568008, XYZ);
                setVector (3 , 0.32351070 , 0.20370246 , -0.92295565, XYZ);
                setVector (4 , 0.34281552 , -0.91708953 , 0.19855809, XYZ);
                setVector (5 , 0.70092911 , -0.31372959 , 0.63896255, XYZ);
                setVector (6 , 0.65047105 , -0.64148449 , 0.40421042, XYZ);
                setVector (7 , 0.70290294 , 0.44259143 , 0.55501466, XYZ);
                setVector (8 , 0.16492432 , 0.10384662 , 0.97980449, XYZ);
                setVector (9 , 0.64148449 , 0.65047105 , 0.40421042, XYZ);
                setVector (10 , 0.12867507 , -0.60056527 , -0.78788645, XYZ);
                setVector (11 , 0.44259143 , -0.70290294 , -0.55501466, XYZ);
                setVector (12 , 0.48597224 , -0.37558875 , -0.78788645, XYZ);
                setVector (13 , 0.27382222 , -0.06222560 , -0.95872330, XYZ);
                setVector (14 , 0.31372959 , 0.70092911 , -0.63896255, XYZ);
                setVector (15 , 0.64148449 , 0.65047105 , -0.40421042, XYZ);
                setVector (16 , 0.18701897 , -0.69772332 , 0.69007759, XYZ);
                setVector (17 , 0.42124579 , -0.85708832 , 0.29317672, XYZ);
                setVector (18 , 0.34390222 , -0.74544747 , 0.56924540, XYZ);
                setVector (19 , 0.44385360 , -0.88955872 , 0.09838828, XYZ);
                setVector (20 , 0.16875696 , -0.46552603 , 0.86764486, XYZ);
                setVector (21 , 0.74279224 , -0.16879818 , 0.64634964, XYZ);
                setVector (22 , 0.59889207 , -0.34781937 , 0.71996595, XYZ);
                setVector (23 , 0.93611167 , -0.13548493 , 0.32146504, XYZ);
                setVector (24 , 0.69772332 , 0.18701897 , 0.69007759, XYZ);
                setVector (25 , 0.75414458 , -0.64747453 , 0.10021817, XYZ);
                setVector (26 , 0.91546263 , -0.38555747 , 0.10618195, XYZ);
                setVector (27 , 0.61036995 , -0.78470967 , 0.09838828, XYZ);
                setVector (28 , 0.89073889 , 0.19453419 , 0.40834015, XYZ);
                setVector (29 , 0.65348308 , 0.41147360 , 0.63375889, XYZ);
                setVector (30 , 0.88955872 , 0.44385360 , 0.09838828, XYZ);
                setVector (31 , 0.24102090 , -0.19469257 , 0.94973929, XYZ);
                setVector (32 , 0.39316608 , -0.42512983 , 0.81405532, XYZ);
                setVector (33 , 0.13568270 , -0.03083365 , 0.98926260, XYZ);
                setVector (34 , 0.39758640 , 0.74834412 , 0.52906252, XYZ);
                setVector (35 , 0.49087882 , 0.43490331 , 0.75359014, XYZ);
                setVector (36 , 0.59511022 , 0.74034766 , 0.30940292, XYZ);
                setVector (37 , 0.75028922 , 0.59085912 , 0.29317672, XYZ);
                setVector (38 , 0.63229297 , 0.52360889 , 0.56924540, XYZ);
                setVector (39 , 0.09696297 , 0.42668316 , 0.89807609, XYZ);
                setVector (40 , 0.30476978 , 0.44117234 , 0.84290174, XYZ);
                setVector (41 , 0.18701897 , -0.69772332 , -0.69007759, XYZ);
                setVector (42 , 0.08417288 , -0.52955884 , -0.84290174, XYZ);
                setVector (43 , 0.23875652 , -0.94818156 , -0.20481223, XYZ);
                setVector (44 , 0.05286024 , -0.99203372 , -0.10524205, XYZ);
                setVector (45 , 0.35174319 , -0.92968916 , -0.09977873, XYZ);
                setVector (46 , 0.90403321 , -0.05597096 , -0.42141690, XYZ);
                setVector (47 , 0.85109765 , -0.41975931 , -0.31214725, XYZ);
                setVector (48 , 0.69504842 , -0.51801762 , -0.49655458, XYZ);
                setVector (49 , 0.60433774 , 0.25471318 , -0.75359014, XYZ);
                setVector (50 , 0.80240514 , -0.50798552 , -0.30999629, XYZ);
                setVector (51 , 0.91546263 , -0.38555747 , -0.10618195, XYZ);
                setVector (52 , 0.39057968 , 0.12195417 , -0.91135926, XYZ);
                setVector (53 , 0.39850663 , 0.25092457 , -0.88103934, XYZ);
                setVector (54 , 0.24539918 , 0.15451859 , -0.95598340, XYZ);
                setVector (55 , 0.89073889 , 0.19453419 , -0.40834015, XYZ);
                setVector (56 , 0.74544747 , 0.34390222 , -0.56924540, XYZ);
                setVector (57 , 0.25995723 , 0.80019165 , -0.53862471, XYZ);
                setVector (58 , 0.13548493 , 0.93611167 , -0.32146504, XYZ);
                setVector (59 , 0.47019216 , 0.54837328 , -0.69007759, XYZ);
                setVector (60 , 0.75028922 , 0.59085912 , -0.29317672, XYZ);
                setVector (61 , 0.63229297 , 0.52360889 , -0.56924540, XYZ);
                setVector (62 , 0.01177102 , 0.36568934 , -0.92958794, XYZ);
                setVector (63 , 0.42512983 , 0.39316608 , -0.81405532, XYZ);
                
            case 6 
                setVector (0 , 0.22137674 , 0.97416289 , 0.00000000, XYZ);
                setVector (1 , 0.77974385 , -0.47293449 , 0.40784004, XYZ);
                setVector (2 , 0.57732144 , 0.81529194 , 0.00000000, XYZ);
                setVector (3 , 0.78102390 , 0.49178125 , -0.38230075, XYZ);
                setVector (4 , 0.44259143 , -0.70290294 , 0.55501466, XYZ);
                setVector (5 , 0.76767670 , -0.02000842 , 0.63896255, XYZ);
                setVector (6 , 0.81596863 , -0.53727691 , 0.20863775, XYZ);
                setVector (7 , 0.99863309 , -0.02707330 , 0.00000000, XYZ);
                setVector (8 , 0.20024280 , 0.88116348 , 0.42597503, XYZ);
                setVector (9 , 0.53727691 , 0.81596863 , 0.20863775, XYZ);
                setVector (10 , 0.24116334 , -0.77921387 , -0.57677292, XYZ);
                setVector (11 , 0.52207219 , -0.82913056 , -0.19489523, XYZ);
                setVector (12 , 0.60056527 , 0.12867507 , -0.78788645, XYZ);
                setVector (13 , 0.46966474 , 0.29573014 , -0.83063814, XYZ);
                setVector (14 , 0.20024280 , 0.88116348 , -0.42597503, XYZ);
                setVector (15 , 0.70549291 , 0.67886335 , -0.19855809, XYZ);
                setVector (16 , 0.00043530 , -0.58651698 , 0.80870183, XYZ);
                setVector (17 , 0.39682224 , -0.82711225 , 0.39549770, XYZ);
                setVector (18 , 0.46944685 , -0.74555347 , 0.47092534, XYZ);
                setVector (19 , 0.50937952 , -0.80897266 , 0.28999440, XYZ);
                setVector (20 , 0.21345942 , -0.53828433 , 0.81405532, XYZ);
                setVector (21 , 0.83395128 , -0.11139878 , 0.53862471, XYZ);
                setVector (22 , 0.65565191 , -0.43796269 , 0.61344132, XYZ);
                setVector (23 , 0.95165386 , -0.21626161 , 0.21351075, XYZ);
                setVector (24 , 0.58651698 , 0.00043530 , 0.80870183, XYZ);
                setVector (25 , 0.82028053 , -0.56047656 , 0.10491364, XYZ);
                setVector (26 , 0.86795349 , -0.49462890 , 0.00000000, XYZ);
                setVector (27 , 0.75477870 , -0.65445407 , 0.00000000, XYZ);
                setVector (28 , 0.94313975 , 0.11131360 , 0.30999629, XYZ);
                setVector (29 , 0.97595326 , 0.03196157 , 0.21093763, XYZ);
                setVector (30 , 0.80897266 , 0.50937952 , 0.28999440, XYZ);
                setVector (31 , 0.29951430 , -0.27877677 , 0.91135926, XYZ);
                setVector (32 , 0.52955884 , 0.08417288 , 0.84290174, XYZ);
                setVector (33 , 0.24539918 , 0.15451859 , 0.95598340, XYZ);
                setVector (34 , 0.34022926 , 0.83943603 , 0.42141690, XYZ);
                setVector (35 , 0.58097388 , 0.49163278 , 0.64712252, XYZ);
                setVector (36 , 0.62509101 , 0.75187377 , 0.20481223, XYZ);
                setVector (37 , 0.47640384 , 0.87175942 , 0.10524205, XYZ);
                setVector (38 , 0.77117370 , 0.60349000 , 0.19771682, XYZ);
                setVector (39 , 0.14742790 , 0.33486156 , 0.92958794, XYZ);
                setVector (40 , 0.42512983 , 0.39316608 , 0.81405532, XYZ);
                setVector (41 , 0.16675767 , -0.85066236 , -0.49655458, XYZ);
                setVector (42 , 0.21345942 , -0.53828433 , -0.81405532, XYZ);
                setVector (43 , 0.32199886 , -0.89576701 , -0.30318179, XYZ);
                setVector (44 , 0.03196157 , -0.97595326 , -0.21093763, XYZ);
                setVector (45 , 0.44385360 , -0.88955872 , -0.09838828, XYZ);
                setVector (46 , 0.83943603 , -0.34022926 , -0.42141690, XYZ);
                setVector (47 , 0.94892066 , 0.01070843 , -0.31214725, XYZ);
                setVector (48 , 0.49163278 , -0.58097388 , -0.64712252, XYZ);
                setVector (49 , 0.69443281 , 0.31144264 , -0.64712252, XYZ);
                setVector (50 , 0.82028053 , -0.56047656 , -0.10491364, XYZ);
                setVector (51 , 0.68636720 , -0.71898906 , -0.09977873, XYZ);
                setVector (52 , 0.51072203 , -0.04474807 , -0.85741567, XYZ);
                setVector (53 , 0.53629970 , 0.33768766 , -0.77223744, XYZ);
                setVector (54 , 0.08286116 , 0.05217454 , -0.99418955, XYZ);
                setVector (55 , 0.82711225 , 0.39682224 , -0.39549770, XYZ);
                setVector (56 , 0.74555347 , 0.46944685 , -0.47092534, XYZ);
                setVector (57 , 0.11139878 , 0.83395128 , -0.53862471, XYZ);
                setVector (58 , 0.28238809 , 0.90272820 , -0.32146504, XYZ);
                setVector (59 , 0.25387224 , 0.52872604 , -0.80870183, XYZ);
                setVector (60 , 0.71530934 , 0.57438237 , -0.39549770, XYZ);
                setVector (61 , 0.45060805 , 0.86629019 , -0.21093763, XYZ);
                setVector (62 , 0.14742790 , 0.33486156 , -0.92958794, XYZ);
                setVector (63 , 0.11414660 , 0.18978485 , -0.97414232, XYZ);
                
            case 7 
                setVector (0 , 0.37639307 , -0.59776982 , -0.70639967, XYZ);
                setVector (1 , 0.56243320 , -0.12781192 , 0.81568008, XYZ);
                setVector (2 , 0.08949872 , -0.90755578 , -0.40784004, XYZ);
                setVector (3 , 0.12781192 , 0.56243320 , -0.81568008, XYZ);
                setVector (4 , 0.35714786 , -0.93297717 , 0.00000000, XYZ);
                setVector (5 , 0.48597224 , -0.37558875 , 0.78788645, XYZ);
                setVector (6 , 0.67886335 , -0.70549291 , 0.19855809, XYZ);
                setVector (7 , 0.93297717 , 0.35714786 , 0.00000000, XYZ);
                setVector (8 , 0.31372959 , 0.70092911 , 0.63896255, XYZ);
                setVector (9 , 0.70549291 , 0.67886335 , 0.19855809, XYZ);
                setVector (10 , 0.03882314 , -0.45322020 , -0.88942971, XYZ);
                setVector (11 , 0.70092911 , -0.31372959 , -0.63896255, XYZ);
                setVector (12 , 0.77921387 , 0.24116334 , -0.57677292, XYZ);
                setVector (13 , 0.16492432 , 0.10384662 , -0.97980449, XYZ);
                setVector (14 , 0.02000842 , 0.76767670 , -0.63896255, XYZ);
                setVector (15 , 0.23088203 , 0.39193019 , -0.88942971, XYZ);
                setVector (16 , 0.16675767 , -0.85066236 , 0.49655458, XYZ);
                setVector (17 , 0.19453419 , -0.89073889 , 0.40834015, XYZ);
                setVector (18 , 0.41147360 , -0.65348308 , 0.63375889, XYZ);
                setVector (19 , 0.26398035 , -0.96349124 , 0.00000000, XYZ);
                setVector (20 , 0.25092457 , -0.39850663 , 0.88103934, XYZ);
                setVector (21 , 0.67224833 , -0.07470960 , 0.73518818, XYZ);
                setVector (22 , 0.54837328 , -0.47019216 , 0.69007759, XYZ);
                setVector (23 , 0.94892066 , 0.01070843 , 0.31214725, XYZ);
                setVector (24 , 0.85066236 , 0.16675767 , 0.49655458, XYZ);
                setVector (25 , 0.80240514 , -0.50798552 , 0.30999629, XYZ);
                setVector (26 , 0.94777332 , -0.31579542 , 0.00000000, XYZ);
                setVector (27 , 0.61259818 , -0.78912893 , 0.00000000, XYZ);
                setVector (28 , 0.98194503 , 0.15092408 , 0.10491364, XYZ);
                setVector (29 , 0.99203372 , 0.05286024 , 0.10524205, XYZ);
                setVector (30 , 0.96349124 , 0.26398035 , 0.00000000, XYZ);
                setVector (31 , 0.47990231 , -0.18036951 , 0.85741567, XYZ);
                setVector (32 , 0.46552603 , 0.16875696 , 0.86764486, XYZ);
                setVector (33 , 0.08286116 , 0.05217454 , 0.99418955, XYZ);
                setVector (34 , 0.05597096 , 0.90403321 , 0.42141690, XYZ);
                setVector (35 , 0.28238809 , 0.90272820 , 0.32146504, XYZ);
                setVector (36 , 0.67750557 , 0.66863145 , 0.30318179, XYZ);
                setVector (37 , 0.45060805 , 0.86629019 , 0.21093763, XYZ);
                setVector (38 , 0.71898906 , 0.68636720 , 0.09977873, XYZ);
                setVector (39 , 0.01177102 , 0.36568934 , 0.92958794, XYZ);
                setVector (40 , 0.11414660 , 0.18978485 , 0.97414232, XYZ);
                setVector (41 , 0.00043530 , -0.58651698 , -0.80870183, XYZ);
                setVector (42 , 0.25092457 , -0.39850663 , -0.88103934, XYZ);
                setVector (43 , 0.25796033 , -0.95990301 , -0.10021817, XYZ);
                setVector (44 , 0.27062615 , -0.82567519 , -0.49297359, XYZ);
                setVector (45 , 0.50937952 , -0.80897266 , -0.28999440, XYZ);
                setVector (46 , 0.74834412 , -0.39758640 , -0.52906252, XYZ);
                setVector (47 , 0.98251159 , -0.14607057 , -0.10646762, XYZ);
                setVector (48 , 0.43490331 , -0.49087882 , -0.75359014, XYZ);
                setVector (49 , 0.74034766 , -0.59511022 , -0.30940292, XYZ);
                setVector (50 , 0.75414458 , -0.64747453 , -0.10021817, XYZ);
                setVector (51 , 0.60349000 , -0.77117370 , -0.19771682, XYZ);
                setVector (52 , 0.47990231 , -0.18036951 , -0.85741567, XYZ);
                setVector (53 , 0.44117234 , -0.30476978 , -0.84290174, XYZ);
                setVector (54 , 0.92483040 , 0.21670153 , -0.30940292, XYZ);
                setVector (55 , 0.85708832 , 0.42124579 , -0.29317672, XYZ);
                setVector (56 , 0.65348308 , 0.41147360 , -0.63375889, XYZ);
                setVector (57 , 0.16879818 , 0.74279224 , -0.64634964, XYZ);
                setVector (58 , 0.21626161 , 0.95165386 , -0.21351075, XYZ);
                setVector (59 , 0.51801762 , 0.69504842 , -0.49655458, XYZ);
                setVector (60 , 0.56038637 , 0.71918457 , -0.40834015, XYZ);
                setVector (61 , 0.47640384 , 0.87175942 , -0.10524205, XYZ);
                setVector (62 , 0.04474807 , 0.51072203 , -0.85741567, XYZ);
                setVector (63 , 0.03083365 , 0.13568270 , -0.98926260, XYZ);
                
            case 8 
                setVector (0 , 0.59776982 , 0.37639307 , -0.70639967, XYZ);
                setVector (1 , 0.81529194 , -0.57732144 , 0.00000000, XYZ);
                setVector (2 , 0.20370246 , -0.32351070 , -0.92295565, XYZ);
                setVector (3 , 0.47293449 , 0.77974385 , -0.40784004, XYZ);
                setVector (4 , 0.52207219 , -0.82913056 , 0.19489523, XYZ);
                setVector (5 , 0.59846051 , -0.55423733 , 0.57677292, XYZ);
                setVector (6 , 0.91225616 , -0.40717281 , 0.00000000, XYZ);
                setVector (7 , 0.82913056 , 0.52207219 , 0.19489523, XYZ);
                setVector (8 , 0.02000842 , 0.76767670 , 0.63896255, XYZ);
                setVector (9 , 0.40717281 , 0.91225616 , 0.00000000, XYZ);
                setVector (10 , 0.29573014 , -0.46966474 , -0.83063814, XYZ);
                setVector (11 , 0.88116348 , -0.20024280 , -0.42597503, XYZ);
                setVector (12 , 0.81596863 , -0.53727691 , -0.20863775, XYZ);
                setVector (13 , 0.86382265 , 0.29736400 , -0.40421042, XYZ);
                setVector (14 , 0.06826509 , 0.97404926 , -0.21111355, XYZ);
                setVector (15 , 0.06222560 , 0.27382222 , -0.95872330, XYZ);
                setVector (16 , 0.31144264 , -0.69443281 , 0.64712252, XYZ);
                setVector (17 , 0.11131360 , -0.94313975 , 0.30999629, XYZ);
                setVector (18 , 0.03196157 , -0.97595326 , 0.21093763, XYZ);
                setVector (19 , 0.44687586 , -0.89347802 , 0.00000000, XYZ);
                setVector (20 , 0.33768766 , -0.53629970 , 0.77223744, XYZ);
                setVector (21 , 0.63851356 , -0.22315866 , 0.73518818, XYZ);
                setVector (22 , 0.69504842 , -0.51801762 , 0.49655458, XYZ);
                setVector (23 , 0.85109765 , -0.41975931 , 0.31214725, XYZ);
                setVector (24 , 0.69443281 , 0.31144264 , 0.64712252, XYZ);
                setVector (25 , 0.71918457 , -0.56038637 , 0.40834015, XYZ);
                setVector (26 , 0.62762561 , -0.60088612 , 0.49297359, XYZ);
                setVector (27 , 0.92483040 , 0.21670153 , 0.30940292, XYZ);
                setVector (28 , 0.95990301 , 0.25796033 , 0.10021817, XYZ);
                setVector (29 , 0.99218297 , -0.04795153 , 0.10618195, XYZ);
                setVector (30 , 0.89347802 , 0.44687586 , 0.00000000, XYZ);
                setVector (31 , 0.51072203 , -0.04474807 , 0.85741567, XYZ);
                setVector (32 , 0.53828433 , 0.21345942 , 0.81405532, XYZ);
                setVector (33 , 0.25995723 , 0.80019165 , 0.53862471, XYZ);
                setVector (34 , 0.34781937 , 0.59889207 , 0.71996595, XYZ);
                setVector (35 , 0.13548493 , 0.93611167 , 0.32146504, XYZ);
                setVector (36 , 0.64747453 , 0.75414458 , 0.10021817, XYZ);
                setVector (37 , 0.38555747 , 0.91546263 , 0.10618195, XYZ);
                setVector (38 , 0.78470967 , 0.61036995 , 0.09838828, XYZ);
                setVector (39 , 0.19469257 , 0.24102090 , 0.94973929, XYZ);
                setVector (40 , 0.03083365 , 0.13568270 , 0.98926260, XYZ);
                setVector (41 , 0.25471318 , -0.60433774 , -0.75359014, XYZ);
                setVector (42 , 0.33768766 , -0.53629970 , -0.77223744, XYZ);
                setVector (43 , 0.15092408 , -0.98194503 , -0.10491364, XYZ);
                setVector (44 , 0.37325185 , -0.79029522 , -0.48385691, XYZ);
                setVector (45 , 0.52973701 , -0.84130347 , -0.09791912, XYZ);
                setVector (46 , 0.63851356 , -0.22315866 , -0.73518818, XYZ);
                setVector (47 , 0.94914598 , -0.29289510 , -0.10646762, XYZ);
                setVector (48 , 0.69039376 , 0.05483151 , -0.71996595, XYZ);
                setVector (49 , 0.66863145 , -0.67750557 , -0.30318179, XYZ);
                setVector (50 , 0.55136780 , -0.67814234 , -0.48385691, XYZ);
                setVector (51 , 0.61036995 , -0.78470967 , -0.09838828, XYZ);
                setVector (52 , 0.29951430 , -0.27877677 , -0.91135926, XYZ);
                setVector (53 , 0.34687291 , -0.35337315 , -0.86764486, XYZ);
                setVector (54 , 0.94818156 , 0.23875652 , -0.20481223, XYZ);
                setVector (55 , 0.99203372 , 0.05286024 , -0.10524205, XYZ);
                setVector (56 , 0.87750452 , 0.43462039 , -0.19771682, XYZ);
                setVector (57 , 0.05597096 , 0.90403321 , -0.42141690, XYZ);
                setVector (58 , 0.41975931 , 0.85109765 , -0.31214725, XYZ);
                setVector (59 , 0.58097388 , 0.49163278 , -0.64712252, XYZ);
                setVector (60 , 0.50798552 , 0.80240514 , -0.30999629, XYZ);
                setVector (61 , 0.38555747 , 0.91546263 , -0.10618195, XYZ);
                setVector (62 , 0.18036951 , 0.47990231 , -0.85741567, XYZ);
                setVector (63 , 0.04795153 , 0.99218297 , -0.10618195, XYZ);
                
                
        end
        
        
    otherwise
        setVector (0 , 0.59776982 , 0.37639307 , -0.70639967, XYZ);
        setVector (1 , 0.81529194 , -0.57732144 , 0.00000000, XYZ);
        setVector (2 , 0.20370246 , -0.32351070 , -0.92295565, XYZ);
        setVector (3 , 0.47293449 , 0.77974385 , -0.40784004, XYZ);
        setVector (4 , 0.52207219 , -0.82913056 , 0.19489523, XYZ);
        setVector (5 , 0.59846051 , -0.55423733 , 0.57677292, XYZ);
        setVector (6 , 0.91225616 , -0.40717281 , 0.00000000, XYZ);
        setVector (7 , 0.82913056 , 0.52207219 , 0.19489523, XYZ);
        setVector (8 , 0.02000842 , 0.76767670 , 0.63896255, XYZ);
        setVector (9 , 0.40717281 , 0.91225616 , 0.00000000, XYZ);
        setVector (10 , 0.29573014 , -0.46966474 , -0.83063814, XYZ);
        setVector (11 , 0.88116348 , -0.20024280 , -0.42597503, XYZ);
        setVector (12 , 0.81596863 , -0.53727691 , -0.20863775, XYZ);
        setVector (13 , 0.86382265 , 0.29736400 , -0.40421042, XYZ);
        setVector (14 , 0.06826509 , 0.97404926 , -0.21111355, XYZ);
        setVector (15 , 0.06222560 , 0.27382222 , -0.95872330, XYZ);
        setVector (16 , 0.31144264 , -0.69443281 , 0.64712252, XYZ);
        setVector (17 , 0.11131360 , -0.94313975 , 0.30999629, XYZ);
        setVector (18 , 0.03196157 , -0.97595326 , 0.21093763, XYZ);
        setVector (19 , 0.44687586 , -0.89347802 , 0.00000000, XYZ);
        setVector (20 , 0.33768766 , -0.53629970 , 0.77223744, XYZ);
        setVector (21 , 0.63851356 , -0.22315866 , 0.73518818, XYZ);
        setVector (22 , 0.69504842 , -0.51801762 , 0.49655458, XYZ);
        setVector (23 , 0.85109765 , -0.41975931 , 0.31214725, XYZ);
        setVector (24 , 0.69443281 , 0.31144264 , 0.64712252, XYZ);
        setVector (25 , 0.71918457 , -0.56038637 , 0.40834015, XYZ);
        setVector (26 , 0.62762561 , -0.60088612 , 0.49297359, XYZ);
        setVector (27 , 0.92483040 , 0.21670153 , 0.30940292, XYZ);
        setVector (28 , 0.95990301 , 0.25796033 , 0.10021817, XYZ);
        setVector (29 , 0.99218297 , -0.04795153 , 0.10618195, XYZ);
        setVector (30 , 0.89347802 , 0.44687586 , 0.00000000, XYZ);
        setVector (31 , 0.51072203 , -0.04474807 , 0.85741567, XYZ);
        setVector (32 , 0.53828433 , 0.21345942 , 0.81405532, XYZ);
        setVector (33 , 0.25995723 , 0.80019165 , 0.53862471, XYZ);
        setVector (34 , 0.34781937 , 0.59889207 , 0.71996595, XYZ);
        setVector (35 , 0.13548493 , 0.93611167 , 0.32146504, XYZ);
        setVector (36 , 0.64747453 , 0.75414458 , 0.10021817, XYZ);
        setVector (37 , 0.38555747 , 0.91546263 , 0.10618195, XYZ);
        setVector (38 , 0.78470967 , 0.61036995 , 0.09838828, XYZ);
        setVector (39 , 0.19469257 , 0.24102090 , 0.94973929, XYZ);
        setVector (40 , 0.03083365 , 0.13568270 , 0.98926260, XYZ);
        setVector (41 , 0.25471318 , -0.60433774 , -0.75359014, XYZ);
        setVector (42 , 0.33768766 , -0.53629970 , -0.77223744, XYZ);
        setVector (43 , 0.15092408 , -0.98194503 , -0.10491364, XYZ);
        setVector (44 , 0.37325185 , -0.79029522 , -0.48385691, XYZ);
        setVector (45 , 0.52973701 , -0.84130347 , -0.09791912, XYZ);
        setVector (46 , 0.63851356 , -0.22315866 , -0.73518818, XYZ);
        setVector (47 , 0.94914598 , -0.29289510 , -0.10646762, XYZ);
        setVector (48 , 0.69039376 , 0.05483151 , -0.71996595, XYZ);
        setVector (49 , 0.66863145 , -0.67750557 , -0.30318179, XYZ);
        setVector (50 , 0.55136780 , -0.67814234 , -0.48385691, XYZ);
        setVector (51 , 0.61036995 , -0.78470967 , -0.09838828, XYZ);
        setVector (52 , 0.29951430 , -0.27877677 , -0.91135926, XYZ);
        setVector (53 , 0.34687291 , -0.35337315 , -0.86764486, XYZ);
        setVector (54 , 0.94818156 , 0.23875652 , -0.20481223, XYZ);
        setVector (55 , 0.99203372 , 0.05286024 , -0.10524205, XYZ);
        setVector (56 , 0.87750452 , 0.43462039 , -0.19771682, XYZ);
        setVector (57 , 0.05597096 , 0.90403321 , -0.42141690, XYZ);
        setVector (58 , 0.41975931 , 0.85109765 , -0.31214725, XYZ);
        setVector (59 , 0.58097388 , 0.49163278 , -0.64712252, XYZ);
        setVector (60 , 0.50798552 , 0.80240514 , -0.30999629, XYZ);
        setVector (61 , 0.38555747 , 0.91546263 , -0.10618195, XYZ);
        setVector (62 , 0.18036951 , 0.47990231 , -0.85741567, XYZ);
        setVector (63 , 0.04795153 , 0.99218297 , -0.10618195, XYZ);
        
        
    end
    
    
    % void mat44_to_quatern( mat44 R ,
    %                    float *qb, float *qc, float *qd,
    %                    float *qx, float *qy, float *qz,
    %                    float *dx, float *dy, float *dz, float *qfac )
    function [qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac] =  mat44_to_quatern( R )
    
    %double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
    %double xd,yd,zd , a,b,c,d ;
    %mat33 P,Q ;
    nothing = 0;
    %NOTE: in matlab conversion, 1 has been added to all the indices
    
    %/* offset outputs are read write out of input matrix  */
    
    %ASSIF(qx,R.m[0,3]) ; ASSIF(qy,R.m[1,3]) ; ASSIF(qz,R.m[2,3]) ;
    qx = R(1,4);
    qy = R(2,4);
    qz = R(3,4);
    
    %/* load 3x3 matrix into local variables */
    
    r11 = R(1,1) ; r12 = R(1,2) ; r13 = R(1,3) ;
    r21 = R(2,1) ; r22 = R(2,2) ; r23 = R(2,3) ;
    r31 = R(3,1) ; r32 = R(3,2) ; r33 = R(3,3) ;
    
    %/* compute lengths of each column; these determine grid spacings  */
    
    xd = sqrt( r11*r11 + r21*r21 + r31*r31 ) ;
    yd = sqrt( r12*r12 + r22*r22 + r32*r32 ) ;
    zd = sqrt( r13*r13 + r23*r23 + r33*r33 ) ;
    
    %/* if a column length is zero, patch the trouble */
    
    if ( xd == 0.0 )
        r11 = 1.0 ; 
        r21 = 0.0 ;
        r31 = 0.0;
        xd = 1.0 ; 
    end
    if ( yd == 0.0 )
        r22 = 1.0 ; 
        r12 = 0.0;
        r32 = 0.0 ; 
        yd = 1.0 ; 
    end
    if ( zd == 0.0 )
        r33 = 1.0 ; 
        r13 = 0.0;
        r23 = 0.0 ; zd = 1.0 ; 
    end
    
    %/* assign the output lengths */
    
    dx = xd; %ASSIF(dx,xd) ; 
    dy = yd; %ASSIF(dy,yd) ; 
    dz = zd; %ASSIF(dz,zd) ;
    
    %/* normalize the columns */
    
    r11 = r11/ xd ; r21 = r21/ xd ; r31 = r31/ xd ;
    r12 = r12/ yd ; r22 = r22/ yd ; r32 = r32/ yd ;
    r13 = r13/ zd ; r23 = r23/ zd ; r33 = r33/ zd ;
    
    %/* At this point, the matrix has normal columns, but we have to allow
    %   for the fact that the hideous user may not have given us a matrix
    %   with orthogonal columns.
    %
    %   So, now find the orthogonal matrix closest to the current matrix.
    %
    %   One reason for using the polar decomposition to get this
    %   orthogonal matrix, rather than just directly orthogonalizing
    %   the columns, is so that inputting the inverse matrix to R
    %   will result in the inverse orthogonal matrix at this point.
    %   If we just orthogonalized the columns, this wouldn't necessarily hold. */
    
    Q(1,1) = r11 ; Q(1,2) = r12 ; Q(1,3) = r13 ; %/* load Q */
    Q(2,1) = r21 ; Q(2,2) = r22 ; Q(2,3) = r23 ;
    Q(3,1) = r31 ; Q(3,2) = r32 ; Q(3,3) = r33 ;
    
    P = mat33_polar(Q) ;  %/* P is orthog matrix closest to Q */
    
    r11 = P(1,1) ; r12 = P(1,2) ; r13 = P(1,3) ; %/* unload */
    r21 = P(2,1) ; r22 = P(2,2) ; r23 = P(2,3) ;
    r31 = P(3,1) ; r32 = P(3,2) ; r33 = P(3,3) ;
    
    %/*                            ( r11 r12 r13 )               */
    %/* at this point, the matrix  ( r21 r22 r23 ) is orthogonal */
    %/*                            ( r31 r32 r33 )               */
    
    %/* compute the determinant to determine if it is proper */
    
    zd = r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  %/* should be -1 or 1 */
    
    if( zd > 0 )             %/* proper */
        qfac =  1.0;
    else                   %/* improper ==> flip 3rd column */
        qfac = -1.0;
        r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
    end
    
    %/* now, compute quaternion parameters */
    
    a = r11 + r22 + r33 + 1.0 ;
    
    if( a > 0.5 )                %/* simplest case */
        a = 0.5 * sqrt(a) ;
        b = 0.25 * (r32-r23) / a ;
        c = 0.25 * (r13-r31) / a ;
        d = 0.25 * (r21-r12) / a ;
    else                        %/* trickier case */
        xd = 1.0 + r11 - (r22+r33) ;  %/* 4*b*b */
        yd = 1.0 + r22 - (r11+r33) ;  %/* 4*c*c */
        zd = 1.0 + r33 - (r11+r22) ;  %/* 4*d*d */
        if( xd > 1.0 )
            b = 0.5 * sqrt(xd) ;
            c = 0.25* (r12+r21) / b ;
            d = 0.25* (r13+r31) / b ;
            a = 0.25* (r32-r23) / b ;
        elseif ( yd > 1.0 )
            c = 0.5 * sqrt(yd) ;
            b = 0.25* (r12+r21) / c ;
            d = 0.25* (r23+r32) / c ;
            a = 0.25* (r13-r31) / c ;
        else 
            d = 0.5 * sqrt(zd) ;
            b = 0.25* (r13+r31) / d ;
            c = 0.25* (r23+r32) / d ;
            a = 0.25* (r21-r12) / d ;
        end
        if( a < 0.0 )
            b=-b ; c=-c ; d=-d; a=-a; 
        end
    end
    
    qb = b; %ASSIF(qb,b) ; 
    qc = c; %ASSIF(qc,c) ; 
    qd = d; %ASSIF(qd,d) ;
    return 
    
    
    
    function Z =  mat33_polar( A )
    
    %mat33 X , Y , Z ;
    %float alp,bet,gam,gmi , dif=1.0 ;
    dif = 1.0;
    %int 
    k=0 ;
    
    X = A ;
    
    %/* force matrix to be nonsingular */
    
    %gam = mat33_determ(X) ;
    gam = det(X);
    while( gam == 0.0 )        %/* perturb matrix */
        gam = 0.00001 * ( 0.001 + mat33_rownorm(X) ) ;
        X(1,1) = X(1,1) + gam ; 
        X(2,2) = X(2,2) + gam ; 
        X(3,3) = X(3,3) + gam ;
        gam = mat33_determ(X) ;
    end
    
    while(1)
        %Y = mat33_inverse(X) ;
        Y = inv(X);
        if( dif > 0.3 )     %/* far from convergence */
            alp = sqrt( mat33_rownorm(X) * mat33_colnorm(X) ) ;
            bet = sqrt( mat33_rownorm(Y) * mat33_colnorm(Y) ) ;
            gam = sqrt( bet / alp ) ;
            gmi = 1.0 / gam ;
        else 
            gam = 1.0; 
            gmi = 1.0;  %/* close to convergence */
        end
        Z(1,1) = 0.5 * ( gam*X(1,1) + gmi*Y(1,1) ) ;
        Z(1,2) = 0.5 * ( gam*X(1,2) + gmi*Y(2,1) ) ;
        Z(1,3) = 0.5 * ( gam*X(1,3) + gmi*Y(3,1) ) ;
        Z(2,1) = 0.5 * ( gam*X(2,1) + gmi*Y(1,2) ) ;
        Z(2,2) = 0.5 * ( gam*X(2,2) + gmi*Y(2,2) ) ;
        Z(2,3) = 0.5 * ( gam*X(2,3) + gmi*Y(3,2) ) ;
        Z(3,1) = 0.5 * ( gam*X(3,1) + gmi*Y(1,3) ) ;
        Z(3,2) = 0.5 * ( gam*X(3,2) + gmi*Y(2,3) ) ;
        Z(3,3) = 0.5 * ( gam*X(3,3) + gmi*Y(3,3) ) ;
        
        dif = abs(Z(1,1)-X(1,1))+abs(Z(1,2)-X(1,2)) ...
            +abs(Z(1,3)-X(1,3))+abs(Z(2,1)-X(2,1)) ...
            +abs(Z(2,2)-X(2,2))+abs(Z(2,3)-X(2,3)) ...
            +abs(Z(3,1)-X(3,1))+abs(Z(3,2)-X(3,2)) ...
            +abs(Z(3,3)-X(3,3))                          ;
        
        k = k+1 ;
        if( k > 100 | dif < 3.e-6 ) 
            break ;  %/* convergence or exhaustion */
            X = Z ;
        end
        
        return %Z 
    end
    
    
    function r1 = mat33_colnorm(  A )  %/* max column norm of 3x3 matrix */
    
    %   float r1,r2,r3 ;
    
    r1 = abs(A(1,1))+abs(A(2,1))+abs(A(3,1)) ;
    r2 = abs(A(1,2))+abs(A(2,2))+abs(A(3,2)) ;
    r3 = abs(A(1,3))+abs(A(2,3))+abs(A(3,3)) ;
    if( r1 < r2 ) 
        r1 = r2 ;
    end
    if( r1 < r3 ) 
        r1 = r3 ;
    end
    return %r1
    
    
    function r1 = mat33_rownorm(  A )  %/* max row norm of 3x3 matrix */
    
    %   float r1,r2,r3 ;
    
    r1 = abs(A(1,1))+abs(A(1,2))+abs(A(1,3)) ;
    r2 = abs(A(2,1))+abs(A(2,2))+abs(A(2,3)) ;
    r3 = abs(A(3,1))+abs(A(3,2))+abs(A(3,3)) ;
    if( r1 < r2 ) 
        r1 = r2 ;
    end
    if( r1 < r3 ) 
        r1 = r3 ;
    end
    return %r1
    
    function R = quatern_to_mat44( qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac )
    
    %mat44 R ;
    a = qb;
    b = qb;
    c = qc;
    d = qd;
    %, xd,yd,zd ;
    
    %/* last row is always [ 0 0 0 1 ] */
    
    R(3,0) = 0.0;
    R(3,1) = 0.0;
    R(3,2) = 0.0 ; 
    R(3,3) = 1.0 ;
    
    %/* compute a parameter from b,c,d */
    
    a = 1.0 - (b*b + c*c + d*d) ;
    if( a < 1.e-7 )                   %/* special case */
        a = 1.0 / sqrt(b*b+c*c+d*d) ;
        b = b * a ; 
        c = c * a ; 
        d = d * a ;        %/* normalize (b,c,d) vector */
        a = 0.0 ;          %/* a = 0 ==> 180 degree rotation */
    else
        a = sqrt(a) ;      %               /* angle = 2*arccos(a) */
    end
    
    %/* load rotation matrix, including scaling factors for voxel sizes */
    
    if (dx > 0.0)
        xd = dx;
    else
        xd = 1.0;
    end
    %xd = (dx > 0.0) ? dx : 1.0 ;       /* make sure are positive */
    if (dy > 0.0)
        yd = dy;
    else
        yd = 1.0;
    end
    %yd = (dy > 0.0) ? dy : 1.0 ;
    if (dz > 0.0)
        zd = dz;
    else
        zd = 1.0;
    end
    %zd = (dz > 0.0) ? dz : 1.0 ;
    
    if( qfac < 0.0 ) 
        zd = -zd ;         %/* left handedness? */
    end
    
    R(1,1) =       (a*a+b*b-c*c-d*d) * xd ;
    R(1,2) = 2.0 * (b*c-a*d        ) * yd ;
    R(1,3) = 2.0 * (b*d+a*c        ) * zd ;
    R(2,1) = 2.0 * (b*c+a*d        ) * xd ;
    R(2,2) =       (a*a+c*c-b*b-d*d) * yd ;
    R(2,3) = 2.0 * (c*d-a*b        ) * zd ;
    R(3,1) = 2.0 * (b*d-a*c        ) * xd ;
    R(3,2) = 2.0 * (c*d+a*b        ) * yd ;
    R(3,3) =       (a*a+d*d-c*c-b*b) * zd ;
    
    %/* load offsets */
    
    R(1,4) = qx ; 
    R(2,4) = qy ; 
    R(3,4) = qz ;
    
    return 
    
    
