function [ out] = imagefilt(ima,filter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% apply 

if nargin<2
    filter = 'hann';
end
ima = abs(ima);
ks = fftshift(fftn(ima));
fks = ks.*sr_filt3d(size(ima,1),filter);
out = ifftn((fks));

end

