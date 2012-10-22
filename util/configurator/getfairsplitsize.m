% gets fair distrib parts



function [split_inds] = getfairsplitsize(samplecount,splits)

splitpart_sizes = ones(splits,1);
baseval = floor(samplecount/splits);
baseval = baseval * splitpart_sizes;

if mod(samplecount,splits)
    wind = 1;
    while sum(splitpart_sizes) ~= samplecount
        splitpart_sizes = baseval;
        splitpart_sizes = round(rand(splits,1)) + splitpart_sizes;
        wind = wind + 1;
        if wind > 100000
            error('problem in achieving a fair split. get in contact with easyupmvpa admins')
        end
    end
else
    splitpart_sizes = baseval;
end



split_inds = zeros(splits,2);
split_inds(1,1) = 1;
split_inds(1,2) = splitpart_sizes(1);

for s=2:splits
    split_inds(s,1) = sum(splitpart_sizes(1:s-1))+1;
    split_inds(s,2) = sum(splitpart_sizes(1:s));
end


    
