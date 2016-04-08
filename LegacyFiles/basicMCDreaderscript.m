filenameMCD = '2014-11-10_M14_RRRI.mcd';
goodChannels = [1 2 3 4]%[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

[pathname, fileName, ext]=fileparts(which('nsMCDLibrary64.dll'));
ns_SetLibrary([pathname filesep fileName ext]);
[nsresult, hfile] = ns_OpenFile(filenameMCD);
[nsresult, FileInfo] = ns_GetFileInfo(hfile);
[nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);
freq = 1/FileInfo.TimeStampResolution;

%[nsresult,data.EventData,~,~] = ns_GetEventData(hfile,2,[1:2:EntityInfo(2).ItemCount]);
%startIDX = data.EventData(1) * 40000;

% data.EventData = data.EventData (2:length(data.EventData),:) - data.EventData(1);
 channelList = goodChannels;
% %reference acquisition
% x = 1;
% tempArray = zeros (72000000,1);
% sumArray = zeros(72000000,1);
% 
% data.chandata = zeros (72000000, length(channelList));
x=1;
data.chandata = nan(EntityInfo(18).ItemCount,length(goodChannels));
%         if (strcmp(method, 'rawOnly'))
for chanNumber = channelList
    [nsresult,~,data.chandata(:,x)]= ns_GetAnalogData(hfile,chanNumber+17,1,EntityInfo(chanNumber+17).ItemCount);
    %data.Nex = nexAddContinuous(data.Nex, 0, freq, tempdata, 'chanNameHere');
    
    x = x+1;
end
ns_CloseFile(hfile)