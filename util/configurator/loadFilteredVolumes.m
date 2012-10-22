% check out a folder and load all volumes which fulfill the filemask and
% the header mask at the same time.
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


function [data,sampleheader,alldescrip] = loadFilteredVolumes(dir_images,mask_file,mask_header);



if nargin == 0 %call the function without argument to manually specify them
    dir_images = uigetdir;
    mask_file=input('File Mask: ','s');
end
if nargin ==1
    mask_file = '*';
    mask_header = [];
end
if nargin == 2
    mask_header = [];
end


imglist=dir(fullfile(dir_images,[mask_file,'*.img']));
niilist=dir(fullfile(dir_images,[mask_file,'*.nii']));
filelist=[imglist;niilist];
filecount = size (filelist,1);
if filecount == 0
    error(['No *.img / *.nii files found. Dir: ',dir_images,' Mask File: ',mask_file,' Mask header: ',mask_header])
end


sampleheader = load_nii(fullfile(dir_images,filelist(1).name));
dimensions = size(sampleheader.img);
if filecount > 1
    
    sampleheader=rmfield(sampleheader,'img');
    index = 1;
    for i = 1:filecount
        tempheader = load_nii_hdr(fullfile(dir_images,filelist(i).name));
        if numel(mask_header) > 0
            if strfind(tempheader.hist.descrip,mask_header);
      
                files2readin(index) = i;
                alldescrip(index,1:size(tempheader.hist.descrip,2)) = tempheader.hist.descrip;
                index = index + 1;
            end
        else
            files2readin(index) = i;
            alldescrip(index,1:size(tempheader.hist.descrip,2)) = tempheader.hist.descrip;
            index = index + 1;
        end
    end
    
    
    volumes = size(files2readin,2);
    data = zeros([dimensions volumes]);
    
    for i = 1:volumes
        niidata = load_nii(fullfile(dir_images,filelist(files2readin(i)).name));
        data(:,:,:,i) = niidata.img;
    end
else
    data = sampleheader.img;
    sampleheader=rmfield(sampleheader,'img');
    alldescrip(1,1:size(sampleheader.hdr.hist.descrip,2)) = sampleheader.hdr.hist.descrip;
end
       
        
%clearvars -except data sampleheader alldescrip

