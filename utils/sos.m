function [ out] = sos(ima);
%UNTITLED2 Summary of this function goes here
%   filter types = 'hann', 'hamming', 'hamming2'

%%
out = sqrt(sum(abs(ima).^2,4));
    
end

