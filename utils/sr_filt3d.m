function [res] = sr_filt3d(n,window)
res=ones(n,n,n);

if ~rem(n,2)
    % Even length window
    half = n/2;
    w = calc_window(half,n,window);
    w = [w; w(end:-1:1)];
else
    % Odd length window
    half = (n+1)/2;
    w = calc_window(half,n,window);
    w = [w; w(end-1:-1:1)];
end

a = w*w';

for k=1:n
res(:,:,k) = a.*w(k);
end

%---------------------------------------------------------------------
function w = calc_window(m,n,window)
%CALC_WINDOW   Calculate the generalized cosine window samples.
%   CALC_WINDOW Calculates and returns the first M points of an N point
%   generalized cosine window determined by the 'window' string.

x = (0:m-1)'/(n-1);

switch window
    
    case 'hann'
    % Hann window
    % w = 0.5 * (1 - cos(2*pi*(0:m-1)'/(n-1)));     
    w = 0.5 - 0.5*cos(2*pi*x);   
case 'hamming'
    % Hamming window
    % w = (54 - 46*cos(2*pi*(0:m-1)'/(n-1)))/100;
    w = 0.54 - 0.46*cos(2*pi*x);
case 'hamming2'
    w = 0.80 - 0.20*cos(2*pi*x);%was 65-35
case 'blackman'
    % Blackman window
    % Force end points to zero to avoid close-to-zero negative values caused
    % by roundoff errors.
    % w = (42 - 50*cos(2*pi*(0:m-1)/(n-1)) + 8*cos(4*pi*(0:m-1)/(n-1)))'/100;
    w = 0.42 - 0.5*cos(2*pi*x) + 0.08*cos(4*pi*x);
    w(1) = 0;    
case 'flattopwin'
    % Flattop window
    % Coefficients as defined in the reference [1] (see flattopwin.m)
    a0 = 0.21557895;
    a1 = 0.41663158;
    a2 = 0.277263158;
    a3 = 0.083578947;
    a4 = 0.006947368;
    w = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x) - a3*cos(6*pi*x) + ...
      a4*cos(8*pi*x);    
end
