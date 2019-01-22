function [image, image_size] = fft3D(k_space, filter, pad_size, os, invert)
%function [image, image_size] = fft3D(k_space, filter, pad_size, os, invert)
% 
% Perform 3D FFT of the give k_space data. Does correctly use the
% fftshift/ifftshift functions and allows error checked zero padding and
% easy application of k_space filters.
%
% input:
%  -    k_space: 3-D Array of (complex) k_space data
%  -   [filter]: string - {hann_sep},'hann_rad', empty argument indicates
%                no filtering
%  - [pad_size]: bool|3-element vector: If scalar, the argument will be
%                interpreted as boolean value (irrespective of the value!!!).
%                If true the data will be zero padded to the next power of
%                2 in each dimension. 
%                If a 3 element vector is provided each element correspond to
%                the new size in this dimension. Zero padding is performed
%                after filtering. Empty argument means no padding.
%   -      [os]: scalar - Will cut os fold oversampling (along first dimension) 
%                from the data. If os is true, 2.0 is the standard os value
%   -  [invert]: Perform inverse fourier transform instead of forward FFT.
%
% output:
%   -        image: 3-D complex image
%   - [image_size]: 2-D vector representing the size of the image
%
% requirements
%   * If filtering is performed needs hannFilt.m
%
% see also
%   hannFilt fft2D reco_cart rawread
%


TAG = '[ FFT3D ]';

if nargin < 1
    error([TAG 'Sorry! Please provide necessary input arguments!']);
end

if ndims(k_space)~=3
    error([TAG 'Sorry! This function only supports 3D data!']);
end

% Check for filter and filter the k-space correspondingly
if nargin > 1 && ~isempty(filter)
   if ischar(filter)
       switch lower(filter)
           case {'hann','hann_sep'}
               k_space = hannFilt(k_space,'sep');
           case {'hann_rad'}
               k_space = hannFilt(k_space,'rad');      
           otherwise
               error([TAG 'Unsupported Filter argument!'])
       end
   else
       error([TAG 'Unsupported Filter argument!'])
   end
end

% Check for zero padding
sz = size(k_space);
if nargin > 2 && ~isempty(pad_size)
    if islogical(pad_size) && pad_size
       pad_size(1) = 2^nextpow2(sz(1));
       pad_size(2) = 2^nextpow2(sz(2));
       pad_size(3) = 2^nextpow2(sz(3));
    elseif numel(pad_size) == 3 
       if any(pad_size < sz)
           error([TAG 'Invalid argument for pad_size!']);
       end
    else    
        error([TAG 'Invalid argument for pad_size!']);
    end
else    
    pad_size = sz;
end

% Perform the FT - 'direction' depends on the invert argument
if nargin < 5 || isempty(invert) || invert == false;
    image = ifftshift(fftn(fftshift(k_space),pad_size));
else
    image = ifftshift(ifftn(fftshift(k_space),pad_size));
end

image_size = size(image);

% Finally remove the oversampling
if (nargin > 3) && ~isempty(os) && ( (islogical(os) && os) || (os > 1.0))
   if islogical(os) && os
       os = 2.0; % Standard is twofold oversampling
   end
   bnd = (os-1)/os*1/2;
   image([1:round(bnd*image_size(1)) round((1-bnd)*image_size(1))+1:end],:,:) = [];
end
   

end