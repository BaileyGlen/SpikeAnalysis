sortVer = 2;
shortFileName = 'bintest03_M16D09';
filename = 'm22D09_03';
delimiter = ',';

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen([filename '.txt'],'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
tempdata = [dataArray{1:end-1}];
% Clear temporary variables
clearvars delimiter formatSpec fileID dataArray ans;
%%
% celltest = cell(1);
% totcount = 1;
% for numchan = unique(tempdata(:,1)).'
%     for numcell=  unique(tempdata(tempdata(:,1)==numchan,2)).'
%         shorttemp = tempdata(tempdata(:,1)==numchan & tempdata(:,2)==numcell,:);
%         celltest{totcount,1} = shorttemp(:,3).';
%         totcount = totcount +1;
%     end
% end
%
% clearvars numtest totcount numchan numcell;

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
mask = data.goodDetectionMask;
%set the first and last .2s to bad to remove filtfilt ripples
mask(1:400*.2) = 0;
mask(end-(400*.1):end) = 0;
%find the transitions to 0, and remove the
zerovals = find(diff(mask) == -1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x)-(.05*400):zerovals(x)) = 0;
end
onevals = find(diff(mask) == +1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x):zerovals(x)+(.05*400)) = 0;
end
BeamBreak_IDX = round(data.BeamBreakTS(data.BeamBreakTS<=1800) .* 400);
goodBeamBreakTS = ones(length(data.BeamBreakTS),1);
% goodRL_TS = ones(length(data.RL_TS),1);
% Previously Missing Initialization %
celltest = cell (1);
%
%         if (animalNumber==3 && dayNumber==1)
%             goodLL_TS(1:end-2) = arrayfun (@(x) min(mask(x-(400*.5):x+(400*.5)))  , data.LL_IDX(1:end-2), 'UniformOutput', true);
%             goodLL_TS(end-1:end) = 0;
%         else
%             goodLL_TS = arrayfun (@(x) min(mask(x-(400*.5):x+(400*.5)))  , data.LL_IDX, 'UniformOutput', true);
%         end
%goodRL_TS = arrayfun (@(x) min(mask(x-(400*.5):x+(400*.5))) , data.RL_IDX, 'UniformOutput', true);
goodBeamBreakTS = arrayfun (@(x) min(mask(x:x+2)) , BeamBreak_IDX, 'UniformOutput', true);
totcount = 1;
for numchan = unique(tempdata(:,1)).'
    for numcell=  unique(tempdata(tempdata(:,1)==numchan,2)).'
        shorttemp = tempdata(tempdata(:,1)==numchan & tempdata(:,2)==numcell,:);
        celltest{totcount,1} = shorttemp(:,3).';
        totcount = totcount +1;
    end
end

nz=max(cellfun(@numel,celltest));
STMtx=cell2mat(cellfun(@(x) [x,nan(1,nz-numel(x))],celltest,'uni',false)).';

%% Get the edited mask
mask = data.goodDetectionMask;
%set the first and last .2s to bad to remove filtfilt ripples
mask(1:400*.2) = 0;
mask(end-(400*.1):end) = 0;
%find the transitions to 0, and remove the
zerovals = find(diff(mask) == -1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x)-(.05*400):zerovals(x)) = 0;
end
onevals = find(diff(mask) == +1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x):zerovals(x)+(.05*400)) = 0;
end
maskup = upsample(mask,100);
for x = 1:720000
    maskup(x*100-99:x*100) = mask (x);
end

mask = data.goodDetectionMask;
%set the first and last .2s to bad to remove filtfilt ripples
mask(1:400*.2) = 0;
mask(end-(400*.1):end) = 0;
%find the transitions to 0, and remove the points before. then do
%to 1
zerovals = find(diff(mask) == -1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x)-(.05*400):zerovals(x)) = 0;
end
onevals = find(diff(mask) == +1)+1;
for x = 1:length(zerovals)
    mask(zerovals(x):zerovals(x)+(.05*400)) = 0;
end
maskup = upsample(mask,100);
for x = 1:720000
    maskup(x*100-99:x*100) = mask (x);
end
mask = maskup;
chandata (~maskup,:) = 0;

%% Save the raw events

%% Save the Event Bursts

%% Set the 5 event types
%        behavEvents_Chris = cell(0);
%        tempCell = cell(0);
m=1;
n=1;
%b=1;
lEvt = cell(length(data.LL_TS),1);
lBEvt = cell(0);
lBTm = [];
for i = 1:length(data.LL_TS)
    if data.LL_Rf_Mask(i) == 0
        lEvt{m} = 'LL_U';
        m = m+1;
    end
    if data.LL_Rf_Mask(i) == 1
        lEvt{m} = 'LL_R';
        m = m+1;
 %       lfind(data.BeamBreakTS>data.LL_TS(i),1);
    end
    if data.burstInfo.LL_First(i) == 1
        lBEvt{n,1} = 'lB_F';
        lBTm(n,1)=data.LL_TS(i);
        n=n+1;
    end
    if data.burstInfo.LL_Last(i) == 1
        lBEvt{n,1} = 'lB_L';
        lBTm(n,1)=data.LL_TS(i);
        n=n+1;
    end
end


m=1;
n=1;
%b=1;
rEvt = cell(length(data.RL_TS),1);
rBEvt = cell(0);
rBTm = [];
for i = 1:length(data.RL_TS)
    if data.RL_Rf_Mask(i) == 0
        rEvt{m} = 'RL_U';
        m = m+1;
    end
    if data.RL_Rf_Mask(i) == 1
        rEvt{m} = 'RL_R';
        m = m+1;
    end
    if data.burstInfo.RL_First(i) == 1
        rBEvt{n,1} = 'rB_F';
        rBTm(n,1)=data.RL_TS(i);
        n=n+1;
    end
    if data.burstInfo.RL_Last(i) == 1
        rBEvt{n,1} = 'rB_L';
        rBTm(n,1)=data.RL_TS(i);
        n=n+1;  
    end
end

m=1;
bbEvt = cell(length(data.BeamBreakTS),1);
for i = 1:length(data.BeamBreakTS)
    bbEvt{m} = 'BB';
    m = m+1;
end




Evt = [lEvt;lBEvt;rEvt;rBEvt;bbEvt;];
EvtTm = [data.LL_TS;lBTm;data.RL_TS;rBTm;data.BeamBreakTS];
[EvtTm,k] = sort(EvtTm);
Evt = Evt(k);


behaveEvt_Raw =Evt;
behaveEvtTm_Raw = EvtTm;

save([shortFileName '_03.mat'], 'STMtx', 'mask', 'behaveEvt_Raw', 'behaveEvtTm_Raw', 'sortVer');