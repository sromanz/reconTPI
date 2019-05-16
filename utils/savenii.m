function [ out ] = savenii( ima,filename,res)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%
if nargin <2
    filename = 'ma_untitled.nii';
    res = 1.0;
end
N = size(ima);
m = max(abs(ima(:)));
normFactor= 8192/m;
normFactor=1.0;
save_nii(make_nii(flipdim(permute(normFactor*abs(squeeze(ima)),[2 1 3]),3),[res res res],N/2 ),filename);

out = 1;
end

