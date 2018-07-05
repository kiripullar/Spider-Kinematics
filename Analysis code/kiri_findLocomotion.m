function [first,last]=kiri_findLocomotion(cut,sum_vector)
% kiri_findLocomotion takes vector containing speed in all frames and a
% minimum value for speeds desired for analysis. Starts at frame in the
% middle of the sequence and works in forward and backward direction to 
% determine the frames where speed drops under the threshold value.
%
% Kiri Pullar, masters thesis 2009

first=2;
last=length(sum_vector)-2;
index=round((first+last)/2);

for n=index:-1:first;
    if(cut>sum_vector(n))
        first=n;
        break;
    end
end

for n=index:last;
    if(cut>sum_vector(n))
        last=n;
        break;
    end
end

if first==last %if first equals last start at 1/3 sequence rather than 1/2
first=2;
last=length(sum_vector)-2;
index=round((first+last)/3);
for n=index:-1:first;
    if(cut>sum_vector(n))
        first=n;
        break;
    end
end

for n=index:last;
    if(cut>sum_vector(n))
        last=n;
        break;
    end
end
end