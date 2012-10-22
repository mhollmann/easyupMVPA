function [outgoing] = manipulateBox(incoming,searchlightdiameter,manipulation)

searchlightradius = floor(searchlightdiameter/2);

switch manipulation
    case 'increase'
        switch size(size(incoming),2)
            case 3
                hsize = size(incoming,1);
                ksize = size(incoming,2);
                lsize = size(incoming,3);
                outgoing = zeros(hsize+2*searchlightradius,ksize+2*searchlightradius,lsize+2*searchlightradius);
                outgoing(searchlightradius+1:hsize+searchlightradius,searchlightradius+1:ksize+searchlightradius,searchlightradius+1:lsize+searchlightradius) = incoming;
            case 4
                hsize = size(incoming,1);
                ksize = size(incoming,2);
                lsize = size(incoming,3);
                volumes = size(incoming,4);
                outgoing = zeros(hsize+2*searchlightradius,ksize+2*searchlightradius,lsize+2*searchlightradius,volumes);
                outgoing(searchlightradius+1:hsize+searchlightradius,searchlightradius+1:ksize+searchlightradius,searchlightradius+1:lsize+searchlightradius,:) = incoming;
        end
        
    case 'decrease'
        switch size(size(incoming),2)
            case 3
                hsize = size(incoming,1);
                ksize = size(incoming,2);
                lsize = size(incoming,3);
                outgoing = zeros(hsize-2*searchlightradius,ksize-2*searchlightradius,lsize-2*searchlightradius);
                outgoing = incoming(searchlightradius+1:hsize-searchlightradius,searchlightradius+1:ksize-searchlightradius,searchlightradius+1:lsize-searchlightradius);
            case 4
                hsize = size(incoming,1);
                ksize = size(incoming,2);
                lsize = size(incoming,3);
                volumes = size(incoming,4);
                outgoing = zeros(hsize-2*searchlightradius,ksize-2*searchlightradius,lsize-2*searchlightradius,volumes);
                outgoing = incoming(searchlightradius+1:hsize-searchlightradius,searchlightradius+1:ksize-searchlightradius,searchlightradius+1:lsize-searchlightradius,:);
        end
end

clearvars -except outgoing