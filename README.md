# reconTPI package

roman.sandro@gmail.com
Tue Jan 22 10:02:34 2019

## Installation

copy the entire folder in your favourite location

## MEX compilation 
compilation steps all inside matlab (you have to have a gcc compiler installed)

1) go into readsiemens folder and type the following commands in order 

mex -c CFLAGS="-g -Wall -std=c99" siemens_meas.c
mex mexLoadSiemensTraces.c siemens_meas.o
mex mexLoadSiemensScandata.c siemens_meas.o

2) go into gridding folder and type 

mex grid3_MAT.c

## Usage


reconTPI(nvoxels, filter , filename, outfolder)

where: 

nvoxel = final image resolution

filter = filter to apply in k-space. It can take one of the following:
         ['none'| 'hann' | 'hamming' | 'hamming2' | 'blackmann' | 'flattopwin']
         default = 'hamming2'

filename = fullpath to meas.dat datafile

outfolder = location to store data

the function saves data in nifti format and returns complex images in a matrix
