function [ima] = reconTPI(nVoxel,filter,filename,outfolder)
%reconTPI reconstruction of TPI data from scanner based on VB17D 
%   Author Dr.Sandro Romanzetti - University of Aachen
%   usage:
%   reconTPI(nvoxels, filter , filename)
%   
%   nvoxel = final image resolution
%   filter = filter to apply in k-space. It can take one of the following:
%            ['none'| 'hann' | 'hamming' | 'hamming2' | 'blackmann' | 'flattopwin']
%            default = 'hamming2'
%   filename = fullpath to meas.dat datafile
%  
%   the function saves data in nifti format and returns complex images in a matrix

%include paths:
currentFolder = pwd;
path(path,strcat(currentFolder,'/readsiemens'));
path(path,strcat(currentFolder,'/gridding'));
path(path,strcat(currentFolder,'/utils'));
path(path,strcat(currentFolder,'/NIfTI_20140122'));


%%
if nargin < 1
    nVoxel = 64;
end

if nargin < 2
    filter = 'hamming2';
end

if nargin < 3
    info = 'Please select binary file to read';
    [file,pathname]=uigetfile('*.dat',info);
    filename = fullfile(pathname, file);
end

if nargin < 4
info = 'Please select folder for reconstructed data';
    outfolder=uigetdir(pathname,info);
end


disp('running rsr_pref');

%% Read meas.dat protocol information

mrprot = readVB17Header(filename);

%% read protocol

ro       = mrprot.Meas.BaseResolution;
nChannels = mrprot.Meas.iMaxNoOfRxChannels;
res      = mrprot.MeasYaps.sWiPMemBlock.adFree{1}; % [mm]
p        = mrprot.MeasYaps.sWiPMemBlock.adFree{2}; % [%] 
fov      = mrprot.MeasYaps.sSliceArray.asSlice{1,1}.dReadoutFOV; % [mm]
kmax =   1/(2*res); %1/mm



nEchos    = mrprot.MeasYaps.lContrasts;
res       = mrprot.MeasYaps.sWiPMemBlock.adFree{1};               % [mm]
fov       = mrprot.MeasYaps.sSliceArray.asSlice{1,1}.dReadoutFOV; % [mm]

%% load trajectory. Attention this works only with sr_ sequences

crd = readTrajectory(filename,mrprot);

data = readTPI(filename);

%% Calculate density compensation function (dcf) based on trajectory

dcf=dcftpi(crd,mrprot);

%% 

rokern = generateRokern(nVoxel); 


%% averaging

datavg = reshape(data, size(data,1), size(data,2)/mrprot.MeasYaps.lAverages,mrprot.MeasYaps.lAverages,[]);
%%
data = squeeze(sum(datavg,3));


%% create array to store image
ima=single(zeros(nVoxel,nVoxel,nVoxel,nEchos));


%% shft 
if isfield(mrprot.MeasYaps.sSliceArray.asSlice{1,1},'sPosition')
    dx = (mrprot.MeasYaps.sSliceArray.asSlice{1,1}.sPosition.dCor)/res;
else 
    dx=0.0;
end

if isfield(mrprot.MeasYaps.sSliceArray.asSlice{1,1},'sPosition')
    dy = (mrprot.MeasYaps.sSliceArray.asSlice{1,1}.sPosition.dTra/res)-4.0;
else 
    dy =0.0;
end

if isfield(mrprot.MeasYaps.sSliceArray.asSlice{1,1},'sPosition')
    dz =(mrprot.MeasYaps.sSliceArray.asSlice{1,1}.sPosition.dSag)/res;
else 
    dz =0.0;
end


%%
kx = reshape(crd(1,:,1,:),1,[]);
ky = reshape(crd(2,:,1,:),1,[]);
kz = reshape(crd(3,:,1,:),1,[]);
%%

for echo = 1:nEchos    
    data_echo(:,:)=data(:,echo:nEchos:end);
    Signal = reshape(data_echo(:,:),1,[]);
    
    % Apply a phase shift based on the values read from the protocol
    Shifted = Signal.*exp(2*pi.*1i.*(kx.*dx + ky.*dy + kz.*dz));
   
    shifted_data(1,:) = double(real(Shifted(:)));
    shifted_data(2,:) = double(imag(Shifted(:)));

    dat=reshape(shifted_data,2,ro,1,[]);
    
    %% regridding
    N=nVoxel*1.5;
    gdat=grid3_MAT(dat,crd,dcf,nVoxel*1.5,1);
    % rearrange gridded data
    
    gdat = squeeze(gdat(1,:,:,:) + 1j*gdat(2,:,:,:));

    % Apply filter
    
    if ~strcmp(filter,'none')
    
        gdat = gdat.*sr_filt3d(size(gdat,1),filter);
    
    end
    
%% Fourier Transform

gdat = fftn(gdat);
gdat = fftshift(gdat,1);
gdat = fftshift(gdat,2);
gdat = fftshift(gdat,3);

%% extract data
   
gdat(rokern > 0) = gdat(rokern > 0) ./ rokern(rokern > 0);
%%
xs = floor(N/2 - N/1.5/2)+1;
xe = floor(N/2 + N/1.5/2);
gdat = gdat(xs:xe,xs:xe,xs:xe);

   
    ima(:,:,:,echo)= gdat; 

%%
outfilename=strcat(outfolder,'/',strrep(file,'.dat',''),'_filter_',filter,'_echo_',num2str(echo),'.nii');
savenii(gdat,strcat(outfolder,'/',oufilename,'.nii',4));
%%

outfilename_imag = strcat(outfolder,'/','im_',strrep(file,'.dat',''),'_filter_',filter,'_echo_',num2str(echo),'.nii');
savenii_imag(gdat,strcat(outfolder,'/',outfilename_imag,'.nii',4));

outfilename_real = strcat(outfolder,'/',strrep(file,'.dat',''),'_filter_',filter,'_echo_',num2str(echo),'.nii');
savenii_real(gdat,strcat(outfolder,'/','re_',outfilename_real,'.nii',4);

end


end
