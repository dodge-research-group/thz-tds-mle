function Data = dataexp
%% DATAEXP loads TDMLE THz data

% Assign path to data
curDir = fileparts(mfilename('fullpath'));
fs = strfind(curDir, filesep);
dataDir = fullfile(curDir(1:fs(end)),...
    'dat', '2017-03-20', 'Air Scans', 'Normal');

% Load data
RawData = thzload(dataDir);
Data.AirScans.RawData = RawData;

% Assign time base and photocurrent
t = [RawData.Time];
Data.AirScans.t = t(:,1);

gain = 5e7;  %gain is 5e7 V/A
Data.AirScans.X = [RawData.Amplitude].*1e12/gain; % in pA

end
