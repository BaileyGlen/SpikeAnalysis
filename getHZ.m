function datasetHZ = getHZ(DirList,datasetHZ, timepoint)
    eventStruct.AnimalIDList={};
    eventStruct.AnimalIDList_Dropped = {};
    numFiles = length(find(~cellfun(@isempty,DirList)));
    for XX=1:numFiles
        %% Load data
        tempdata = struct();
        load(DirList{XX});
        uniqueIDs = unique(data.spikes.UID);
        %% create temp data struct for appending
        tempdata.SessionType = cell(length(uniqueIDs),1);
        tempdata.timepoint = cell(length(uniqueIDs),1);
        tempdata.AnimalID = cell(length(uniqueIDs),1);
        tempdata.CellID = [1:length(uniqueIDs)]';

        tempdata.SessionType(cellfun('isempty',tempdata.SessionType)) = {data.SessionType}; 
        tempdata.timepoint(cellfun('isempty',tempdata.timepoint)) = {timepoint};
        tempdata.AnimalID(cellfun('isempty',tempdata.AnimalID)) = {data.AnimalID};

        goodTime = sum(data.goodDetectionMask)/400;
        tempdata.meanHZ = arrayfun (@(x) length(data.spikes(data.spikes.UID==x,:)),uniqueIDs)/goodTime;

        datasetHZ = [datasetHZ;struct2dataset(tempdata)];
    end
end