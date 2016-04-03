function fileJoiner(filenameTXT,filenameMAT)                                 
%fileJoiner Joins spikesorting TXT and MCD
%   This version currently requires a header. For the header Cell and
%   delimiter, you can leave these empty []

% if isempty(headerCell)
%headerCell = {'Channel','Unit','Timestamp'};
% end
% if isempty(delimiter)
% end

% fileID = fopen(filenameTXT,'r');
% assert(fileID>=3, 'Did not Open Correctly');
% temptext = strsplit(fgetl(fileID),delimiter);
% headerList = parseHeader(temptext, headerCell);
%
%     function list = parseHeader(fileHeaders, templateHeaders)
%         list = cellfun(@(x) find(strcmp(x,fileHeaders)),templateHeaders, ...
%             'UniformOutput', false);
%         assert(~isempty(list), ...
%             ['Not all templateHeaders were found in the fileHeader.' ...
%             'Is there a header?']);


display(['Joining ' filenameMAT]);

% http://regexr.com/
% Parse the Filename
[startIDX, endIDX] = regexp(filenameTXT,'_M\d{2,3}_');
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
