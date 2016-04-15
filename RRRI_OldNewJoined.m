function RRRI_OldNewJoined()
%RRRI_OLDNEWJOINED joining "ChrisStyle" and "OldStruct"

%Check the OS and get the relative paths - Reuse this code!
if ispc 
    userdir= getenv('USERPROFILE');
    boxdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\'];
else
    userdir= getenv('HOME');
    boxdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/'];
end
    datadir = [boxdir '_mat'];
    outputdir = [boxdir 'z3s_press'];



% Get the current folder name
originaldir = cd;
%set to the data dir
cd (datadir);
%fileNameStruct = dir(['*D' timepointList{timepointIDX} '*.mat']);
% Get the names of the files in the directory
fileStruct = dir('*.mat');

%Enter the main loop, once per file
for fileIDX = 1:length(fileStruct)
    %Get the Filename
    fileName = fileStruct(fileIDX).name;
    %update display
    display(['Joining ' fileName]);
    
    % opem up the existing data struct for the current animal/session
    load([datadir 'OldStruct\' fileName]);
    data.Lapish = load([datadir 'ChrisStyle\' fileName]);
    spikeDataset = importSpikesTXT([datadir 'Spikes\' fileName]);
    
    %Add data.SessionType
    
    data.SessionType = 'RRRI';
    
    %Add data.SessionLength
    
    data.SessionLength = 30;
    
    % Add data.AnimalID
    % format - MXX
    data.AnimalID = fileName(startIDX:endIDX);
    % Add data.DayVar
    % format - DXX

    % Add spikeDataset
    spikeDataset = sortrows(spikeDataset);
    
    % put the raw spike dataset into the data struct
    data.spikes=spikeDataset;
    
    % resave the dataStruct
    save([datadir fileName '_final.mat'],'data');
    
    display(['Completed joining ' filenName]);
    
end

cd (originaldir);