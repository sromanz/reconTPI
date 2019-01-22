function data=readTPI(filename)
%% read all data in from data file data = readoutlength,[ntrajectories,channels] complex
data = mexLoadSiemensTraces(filename);
%% somehow the data reader detects two channels. take the second one.
data=data(:,:,2);
