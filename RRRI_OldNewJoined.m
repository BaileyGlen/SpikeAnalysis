%function RRRI_OldNewJoined()
%RRRI_OLDNEWJOINED joining "ChrisStyle" and "OldStruct"

%Check the OS and get the relative paths - Reuse this code!
if ispc 
    userdir= getenv('USERPROFILE');
    boxdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\_mat\'];
    outputdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\z3s_press\'];
else
    userdir= getenv('HOME');
    boxdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/_mat/'];
    outputdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/z3s_press/'];
end
    datadir = [boxdir 'OldStruct'];




% Get the current folder name
originaldir = cd;
%set to the data dir
cd (datadir);
%fileNameStruct = dir(['*D' timepointList{timepointIDX} '*.mat']);
% Get the names of the files in the directory
fileStruct = dir('*.mat');

if isempty (fileStruct) 
   error ('The fileStruct is empty');
end

%Enter the main loop, once per file
for fileIDX = 1:length(fileStruct)
    %Get the Filename
    fileName = fileStruct(fileIDX).name;
    %update display
    display(['Joining ' fileName]);
    
    % opem up the existing data struct for the current animal/session
    load([boxdir 'OldStruct\' fileName]);
    data.Lapish = load([boxdir 'ChrisStyle\' fileName]);
    spikeDataset = importSpikesTXT([boxdir 'Spikes\' fileName(1:end-4) '.txt']);
    
    %Add data.SessionType
    
    data.SessionType = 'RRRI';
    
    %Add data.SessionLength
    
    data.SessionLength = 30;
    
    % Add data.AnimalID
    % format - MXX
    data.AnimalID = ['M' fileName(2:3)];
    % Add data.DayVar
    % format - DXX
    data.timepoint = fileName(5:6);
  

    % Add spikeDataset
    spikeDataset = sortrows(spikeDataset);
    
    % put the raw spike dataset into the data struct
    data.spikes=spikeDataset;
    
    % resave the dataStruct
    save([datadir fileName(1:end-4) '_final.mat'],'data');
    
    display(['Completed joining ' fileName]);
    
end

cd (originaldir);