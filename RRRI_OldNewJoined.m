function fileJoiner(filenameTXT,filenameMAT)
%RRRI_OLDNEWJOINED joining "ChrisStyle" and "OldStruct"

%Check the OS and get the relative paths
if ispc; userdir= getenv('USERPROFILE');
else userdir= getenv('HOME');
end
boxdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\'];
datadir = [userdir boxdir '_mat'];
outputdir = [userdir boxdir 'z3s_press'];
% Get the current folder name
originalFolder = cd;
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
    
    %Add data.SessionType
    
    data.SessionType = 'RRRI';
    
    %Add data.SessionLength
    
    data.SessionLength = 30;
    
    %Add data.AnimalID
    % Add data.DayVar
    
    
    % http://regexr.com/
    % Parse the Filename
    [startIDX, endIDX] = regexp(filenameTXT,'_[mM]\d{2,3}_');
    startIDX = startIDX+1;
    endIDX = endIDX-1;
    % Read in spike data from offline spike sorter as a dataset
    spikeDataset = importSpikesTXT(filenameTXT);
    spikeDataset = sortrows(spikeDataset);
    
    % opem up the existing data struct for the current animal/session
    load(filenameMAT);
    
    %Save AnimalID
    data.AnimalID = filenameMAT(startIDX:endIDX);
    
    % put the raw spike dataset into the data struct
    data.spikes=spikeDataset;
    
    % put the lapisch type spike mtx into data struct
    data.Lapish.STMtx = getSTMtx(spikeDataset);
    
    
    % resave the dataStruct
    save([filenameMAT(startIDX:end-4) '_final.mat'],'data');
    
    display(['Completed joining ' filenameMAT]);
    
end
% Helper Function to get chris sytle spike MTX
    function STMtx = getSTMtx(spikeDataset)
        uniqueIDs = unique(spikeDataset.UID);
        numUniqueCells = length(uniqueIDs);
        maxSpikesPerCell = max(arrayfun (@(x) length(spikeDataset(spikeDataset.UID==x,:)),uniqueIDs));
        STMtx = nan(maxSpikesPerCell,numUniqueCells);
        temp=arrayfun (@(x) spikeDataset.Timestamp(spikeDataset.UID==x,:),uniqueIDs,'UniformOutput',false);
        for x=1:numUniqueCells
            y=length(temp{x});
            STMtx(1:y,x)=temp{x};
        end
    end

end
