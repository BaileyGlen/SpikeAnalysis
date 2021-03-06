clear all;
cd 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat';
scalefactor = (10^8);


% Read txt into cell A
fid = fopen('singleofi.ofi','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);






animalArray = [13 14 16 18 21 22 24 25];
for animalIDX = 1:length(animalArray)
    
    cd (['D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\' num2str(animalArray(animalIDX))]);
    fileList = dir('02_min*.mat');
    for sessionIDX = 1:length(fileList)
        display(fileList(sessionIDX).name);
        load(fileList(sessionIDX).name);
        display(['raw_' fileList(sessionIDX).name(12:end)]);
        load(['raw_' fileList(sessionIDX).name(12:end)]);
        
        channelList = data.ChannelList;
        %filename_MedPC = 'D:\System Files\Users\Mouse Maze.MouseMaze-PC\Documents\Dropbox\Matlab\test1_7_16,mcd';
        % filenameMCD = 'D:\Users\Bailey\Dropbox\Matlab\test1_7_16.mcd';
        % %data= readRawFiles(channelList , filename_MCD, filename_MedPC);
        % [pathname, fileName, ext]=fileparts(which('nsMCDLibrary64.dll'));
        % ns_SetLibrary([pathname filesep fileName ext]);
        % [nsresult, hfile] = ns_OpenFile(filenameMCD);
        % [nsresult, FileInfo] = ns_GetFileInfo(hfile);
        % [nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);
        % freq = 1/FileInfo.TimeStampResolution;
        
        % [nsresult,data.EventData,~,~] = ns_GetEventData(hfile,2,[1:2:EntityInfo(2).ItemCount]);
        % startIDX = data.EventData(1) * 40000;
        
        %data.EventData = data.EventData (2:length(data.EventData),:) - data.EventData(1);
        %channelList = [1 2 5 6 7 10 11 12 16];
        %reference acquisition
        
        % chanData = zeros (EntityInfo(18).ItemCount, length(channelList));
        % x=1;
        % %         if (strcmp(method, 'rawOnly'))
        % for chanNumber = channelList
        %     [nsresult,~,chanData(:,x)]= ns_GetAnalogData(hfile,chanNumber+17,1,EntityInfo(18).ItemCount);
        %     %data.Nex = nexAddContinuous(data.Nex, 0, freq, tempdata, 'chanNameHere');
        %
        %     x = x+1;
        % end
        chandata = bsxfun(@minus,chandata, mean(chandata,2));
        %         freq = 40000;
        %
        load ('somebandpass.mat');
        for x = 1:size(chandata,2)
            chandata(:,x) = filtfilt(SOS,G,chandata(:,x));
        end
        %chandata = filter(Hd,chandata);
        %try 11 and 12, both should work
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
        
        chandata (~maskup,:) = 0;
        chandata = int16(chandata*scalefactor).';
        fid = fopen(['bintest03_' fileList(sessionIDX).name(12:end-4) '.ofi'], 'w');
        A(1,2) = cellstr(sprintf('MaxMV=%.6f',(2^15-1)*scalefactor*1000));
        A(1,3) = cellstr([ 'NChan=' num2str(length(data.ChannelList))]);
        %
        % for chanNumber = channelList
        %     data.Nex = myNexAddContinuous(data.Nex, 0, 40000, data.chanData(:,x), sprintf('Chan %d',chanNumber));
        %     x = x+1;
        % end
        % myWriteNexFile(data.Nex, '2014-11-04_M13_RRRI.myNex');
        %         fwrite(fid, chandata,'int16');
        % for x = 1:EntityInfo(18).ItemCount
        %     fwrite(fid, int16(chanData(x,:)/(1*10^-9)), 'int16');
        % end
        for i = 1:numel(A)
            if A{i+1} == -1
                fprintf(fid,'%s', A{i});
                break
            else
                fprintf(fid,'%s\n', A{i});
            end
        end
        fclose (fid);
        fid = fopen(['bintest03_' fileList(sessionIDX).name(12:end-4) '.bin'], 'w');
        fwrite(fid, chandata,'int16');
        fclose (fid);
        %clear chandata;
        %clear data;
    end
end