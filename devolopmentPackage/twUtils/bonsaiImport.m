function [scaleFactor,vid_ts,vid_x,vid_y,vid_headdir,d_epoch] = bonsaiImport(basePath,fileName)
% Takes in tracking saved as csv files from Bonsai
%
% Structure of file:
% Centroid.X Centroid.Y Orientation MajorAxisLength MinorAxisLength Area
% ###        ###        ###         ###             ###             ###
% ....
% where ### is a number
%
% Note: Each line (excluding first line with variable names) contains a
% stray space at the end. This causes matlab to read 7 columns instead of 6
% because the file is space delimited. The last column, however, contains
% no data.
%
% basePath: full file path where the csv file is contained
% fileName: full name of csv file WITH .csv extension


% default overwrite to 0
if ~exist('overwrite','var'),
  overwrite = 0;
end

%% Read in file excluding the first line (with variable names) and assign variables
posdata = dlmread(fullfile(basePath,fileName),' ',1,0); 
posdata(:,end) = []; 
vid_x = posdata(:,1);
vid_y = posdata(:,2);
vid_headdir = posdata(:,3);

metadata = dlmread(bonsaifile,' ',1,0); 
metadata(:,end) = []; 
frames = metadata(:,2);
ts = metadata(:,3);

% get bad frames
badinds = find(frames == 0);
vid_x(badinds) = [];
vid_y(badinds) = [];
vid_headdir(badinds) = [];
ts(badinds) = [];
frames(badinds) = [];

scaleFactor = 1; % ???

% extract time data from flea3 timestamps
ts = dec2bin(ts,32); % 32 bit time stamps
sec = bin2dec(:,1:7); % first 7 bits are seconds
cycles = bin2dec(:,8:20); % next 13 bits are cycles (8 kHz)
offset = bin2dec(:,21:28); % last 12 bits are cycle offsets, but last 4 bits are inaccurate (1/125 us)

wrapinds = find(diff(sec)<0); % find where time stamps wrap (wraps every 128 seconds)
wrapinds = [0; wrapinds; length(sec)];
vid_ts = zeros(length(sec),1);
for iW = 1:length(wrapinds)-1,
  for iN = wrapinds(iW)+1:wrapinds(iW+1),
    vid_ts(iN) = 128*(iW-1)+sec(iN)+cycles(iN)/8000;
  end
end
vid_ts = vid_ts-vid_ts(1);

d_epoch = [0 vid_ts(end)];


