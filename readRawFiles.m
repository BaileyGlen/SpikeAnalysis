function [ data  ] = readRawFiles( badChannels, varargin)
%READRAWFILES Reads inn .mcd and MedPC and returns data struct
%   Currently set up to read 30min RIRR sessions in any order



%% Parse Inputs -- Raw Only isn't working, this is kindof a mess but functional
%% !!Needs updating to select run type or id or run type!!
if nargin >= 3
    filenameMCD = varargin{1};
    filenameMedPC = varargin{2};
    if nargin == 4 && strcmp (varargin{3}, 'raw')
        method='rawOnly';
    end
elseif nargin == 1
    [filename,pathname] = uigetfile('*.mcd');
    filenameMCD = [pathname filename];
    testfileName=[filenameMCD(1:end-4) '.txt'];
    if exist(testfileName,'file')
        filenameMedPC = [testfileName];
    else
        [filename,pathname] = uigetfile('*.txt');
        filenameMedPC = [ pathname filename];
    end
else
    display ('You did not provide a correct number of inputs');
end

data.MCD = filenameMCD;
data.MedPCFile = filenameMedPC;


%% FileName to settings
% This loop is criticl for handling multiple types of files. Based upon the
% filename, it will define a variety of settings
if strfind (filenameMedPC,'Deval')
    data.SessionType='Deval';
    numLevers=1;
    data.SessionLength=10;
elseif strfind (filenameMedPC,'ContDeg')
    data.SessionType='ContDeg';
    numLevers=1;
    data.SessionLength=15;
elseif (strfind (filenameMedPC,'RRRI') || strfind (filenameMedPC,'RIRR'))
    data.SessionType='RRRI';
    numLevers=2;
    data.SessionLength=30;
else error('Unknown File ID, check the name');
end

%% Default Variables
totChannelList = 1:1:16; %should probably be calculated by readMCD
data.ChannelList =  totChannelList(~ismember(totChannelList, badChannels));
OpenedFile = filenameMCD(1:end-4);
outputVer = '_01';
newMCDStructName=[OpenedFile outputVer '.mat'];
newMCDRawName=[OpenedFile '_raw' outputVer '.mat'];
newBinName=[OpenedFile outputVer '.bin'];
newOfiName=[OpenedFile outputVer '.ofi'];

%% ---------Generating data struct from Files------------------
% Because these functions are nested functions, the data struct can be
% edited directly

%Function calls
readMedPC();
readMCD();
IDClippingArtifact(); %remove this for future files
MCDMedPCSync();

%% ----------Saving Files-------------------------------------------------
%% Create bin file and ofi file
WriteBinary()
%% save data Struct
save(newMCDStructName, 'data');
%% ----------END OF MAIN FUNCTION-----------------------------------------





%% Nested Functions
%These are nested functions so they can get variables from the parent
    function readMedPC ()
        %READMEDPC Helper Function for readRawFiles. Nested.
        %   This would be much smarter if it looked for the name of the
        %   variable, instead of having it predefined.
        % current version only accepts a single file name
        %  val = isspace(test); gets logical array for spaces. if first char is a
        %  space, then the line is a continuation of previous variable
        % need to make a struct for each file that contains...
        % Start Time (3)
        % End Time (3)
        % File run (MSN) (1) i think
        % Single variables (A B D F and sometimes J(licks))
        % skip to arrays (C
        
        %% Open MedPC file
        fid = fopen(filenameMedPC);
        %% start parsing
        %% Text Variables
        %Start Date
        %Cont Deg
        %15 min 1 lever
        %reward is delivered, unpaired from presses
        %U sets delivery rate, E is delivery timestamp
        
        %deval
        %10 min 1 lever
        %no reward delivered
        
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'Start Date'));
                    data.StartDate = strtrim(curLineSplit{2});
                    success = 1;
                end
            end
        end
        %End Date
        success = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'End Date'));
                    data.MedPC.EndDate = strtrim(curLineSplit{2});
                    success = 1;
                end
            end
        end
        
        %Start Time
        success = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'Start Time'));
                    data.MedPC.StartTime{1} = str2double(curLineSplit{2});
                    data.MedPC.StartTime{2} = str2double(curLineSplit{3});
                    data.MedPC.StartTime{3} = str2double(curLineSplit{4});
                    success = 1;
                end
            end
        end
        
        %End Time
        success = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'End Time'));
                    data.MedPC.EndTime{1} = str2double(curLineSplit{2});
                    data.MedPC.EndTime{2} = str2double(curLineSplit{3});
                    data.MedPC.EndTime{3} = str2double(curLineSplit{4});
                    success = 1;
                end
            end
        end
        
        
        %MSN (also extract from this which lever is done first)
        success = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'MSN'));
                    data.MedPC.MSN = strtrim(curLineSplit{2});
                    success = 1;
                    %%% Extract here the actual protocol for future deviations%%%
                end
            end
        end
        curLine = fgetl(fid);
        while ~stop
            parseNextVar()
        end
        
        
        function parseNextVar()
            success=0;
            while ~stop && ~success
                if curLine == -1;
                    stop = 1;
                else
                    curLineSplit = strsplit(curLine, ':');
                    if ( length(curLineSplit{1})==1 && isletter(curLineSplit{1}))
                        curFieldName=curLineSplit{1};
                        %                                  above line needs to be opened. recognize
                        %                                  var, but not predetermined
                        %                                  Else needs to be flexibile enough to
                        %                                  handle array type or value type
                        tempDouble = str2double(curLineSplit{2});
                        if isnumeric(tempDouble) && ~isnan(tempDouble)
                            data.MedPC.(curFieldName)=tempDouble;
                            success = 1;
                            curLine = fgetl(fid);
                        else
                            tempdataArray = zeros(20000,1);  %why 20000????
                            arrayX = 1;
                            while (~stop && ~success)
                                curLine = fgetl(fid);
                                if curLine == -1
                                    data.MedPC.(curFieldName) = nonzeros (tempdataArray);
                                    stop = 1;
                                    success=1;%!!!!!!!!This is busted above, needs
                                    %                             checking!!!!!!!!!!
                                else
                                    curLineSplit = strsplit(curLine);  %white space this time
                                    if isempty(curLineSplit{1})
                                        
                                        for x = 1:length(curLineSplit) - 2
                                            tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                            arrayX = arrayX +1;
                                        end
                                    else
                                        data.MedPC.(curFieldName) = nonzeros (tempdataArray);
                                        success = 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        %% close file
        fclose(fid);
    end
    function readMCD
        %% ----------Open .mcd----------------
        %[pathname, fileName, ext]=fileparts(which('nsMCDLibrary64.dll'));
        [pathname, fileName, ext]=fileparts(which('nsMCDLibrary.so'));
        if ns_SetLibrary([pathname filesep fileName ext]) < 0 
            [pathname, fileName, ext]=fileparts(which('nsMCDLibrary64.dll'));
            if ns_SetLibrary([pathname filesep fileName ext]) < 0 
                error('MCDImport not working, missing dll');
            end
        end
        [nsresult, hfile] = ns_OpenFile(filenameMCD);
        [nsresult, FileInfo] = ns_GetFileInfo(hfile);
        [nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);
        %% ----------Read .mcd----------------
        % Get Event TS
        
        
        %there are 33 entities
        %the 1st entity is the sync pulse signal
        %the 18th-33rd are the individual channels
        data.Freq = 1/FileInfo.TimeStampResolution;
        [nsresult,data.EventData,~,~] = ns_GetEventData(hfile,3,[1:2:EntityInfo(3).ItemCount]);
        startIDX = data.EventData(1) * data.Freq;
        % Read in Sync Single
        data.EventData = data.EventData (2:length(data.EventData),:) - data.EventData(1);
        timeIdx=data.Freq*data.SessionLength*60;
        % Read in data for good channels
        chandata = nan(timeIdx,length(data.ChannelList));
        for x = 1:length(data.ChannelList)
            [~,~,chandata(:,x)]= ns_GetAnalogData(hfile,data.ChannelList(x)+17,startIDX,timeIdx);
        end
        %common average referencing
        data.refChan=mean(chandata,2);
        data.downDataRaw=downsample(chandata,100);
        %% memory fix
        % because of memory contraints, it is best to save the raw data to
        % make room for filtering and spectral analysis. The data can later
        % be reloaded to subtract off the artifact and create the .bin
        save(newMCDRawName, 'chandata', '-v7.3');
        clear ('chandata');
        %% hard Coded array declarations, needs to change
        %based on parssing of filename
        
        
        %         tempArray = zeros (timeIdx,1);
        %         sumArray = zeros(timeIdx,1);
        %
        %         for chanNumber = data.ChannelList
        %             [nsresult,~,tempArray]= ns_GetAnalogData(hfile,chanNumber+17,startIDX,timeIdx);
        %             sumArray = sumArray + tempArray;
        %             %data.Nex = nexAddContinuous(data.Nex, 0, data.Freq, tempdata, 'chanNameHere');
        %             x = x+1;
        %         end
        %         data.refChan = sumArray/length(data.ChannelList);
        %         clear ('sumArray');
        %         data.chandata = zeros (timeIdx, length(data.ChannelList));
        
        %tempSum=[];
        %         if (strcmp(method, 'rawOnly'))
        
        %Hd = load('strong bandpass coeff.mat');
        %         else
        %             for chanNumber = channelList
        %                 [nsresult,~,tempArray]= ns_GetAnalogdata(hfile,chanNumber+17,startIDX,72000000);
        %                 %data.Nex = nexAddContinuous(data.Nex, 0, data.Freq, tempdata, 'chanNameHere');
        %                 tempArray = tempArray - data.refChan;
        %                 data.chandata(:,x) = filtfilt(Hd.SOS,Hd.G, tempArray);
        %                 x = x+1;
        %             end
        %         end
        %         %reference subtraction
        %         meandata = mean(data.chandata,2, 'single');
        %         for x = 1:length(channelList)
        %             data.chandata (:,x) = data.chandata (:,x) - meandata;
    end
    function MCDMedPCSync()
        %For combinding event data between file types
        %Final goal is to...
        %  1)refer to a med pc time point (reward etc)
        %  2)convert that into a time point in the recording file.
        %This is based on the lever press timestamp data
        
        %% ID first lever
        if numLevers==2
            if data.MedPC.C(2) < data.MedPC.G(2); firstArrayTemp = data.MedPC.C; data.FirstArray = 'LL';
            else firstArrayTemp = data.MedPC.G; data.FirstArray = 'RR';
            end;
        else
            if ~isempty(data.MedPC.C)
                data.FirstArray = 'LL';
                firstArrayTemp = data.MedPC.C;
            elseif ~isempty(data.MedPC.G)
                data.FirstArray = 'RR';
                firstArrayTemp = data.MedPC.G;
            else error('There are no lever presses');
            end
        end
        %% calculate linear offset
        lastIDXOfFirst = length(firstArrayTemp);
        if isempty(data.EventData)
            temp = load('defaultOffset.mat');
            data.MedPC.offsetP = temp;
            clear temp;
        else
            data.MedPC.offsetP = polyfit (firstArrayTemp, data.EventData(1:lastIDXOfFirst,1),1);            
        end
        data.Lapish.behaveEvt_Raw ={};
        data.Lapish.behaveEvtTm_Raw = [];
        curLever=data.FirstArray;
        curIDX = 1;
        lastIDX= lastIDXOfFirst;
        for x = 1:numLevers
            if x==2
                if strcmp(curLever,'LL')
                    curLever='RR';
                else
                    curLever='LL';
                end
                curIDX = 1+(lastIDXOfFirst);
                lastIDX=length(data.EventData);
            end
            BaileyDataType();

            LapishDataType();
        end
        function BaileyDataType()
            % set up the key variables
            %Presses extracted from events
            BBVar='K';
            if strcmp(curLever,'LL')
                LeverVar='C';
                RfVar='I';
            else
                LeverVar='G';
                RfVar='E';
            end
            %Need to remove second lever press
            
            %            data.([curLever '_TS'])=data.EventData(lastIDXOfFirst+1:length(data.EventData),1);
            %                 data.LL_TS = data.EventData(1:lastIDXOfFirst);
            %                 data.RL_TS = data.EventData(lastIDXOfFirst+1:length(data.EventData),1);
            %
            %                 data.LL_TS = data.EventData(lastIDXOfFirst+1:length(data.EventData),1);
            %                 data.RL_TS = data.EventData(1:lastIDXOfFirst);
            %             end
            % get Press TS from event array
            if isempty(data.EventData)
                 data.([curLever '_TS']) =  data.MedPC.offsetP(1)*data.MedPC.(LeverVar) +  data.MedPC.offsetP(2);
            else
                data.([curLever '_TS'])=data.EventData(curIDX:lastIDX,1);
            end
            data.([curLever '_Rf_TS']) =  data.MedPC.offsetP(1)*data.MedPC.(RfVar) +  data.MedPC.offsetP(2);
            
            LevIDX = round(data.([curLever '_TS']) .* 400);
            RfIDX = round(data.([curLever '_Rf_TS']) .* 400);
            if ~strcmp(data.SessionType,'ContDeg')
                data.([curLever '_Rf_Mask'])=ismember(round(LevIDX) ,round(RfIDX));
                %data.([curLever '_IDX'])=round(data.([curLever '_TS']) .* 400);
                %data.([curLever '_Rf_IDX'])=round(RefTS .* 400);
                %if these are done at 40000, they don't come out correctly with
                %RF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                %data.LL_IDX = round(data.LL_TS .* 400);
                %data.RL_IDX = round(data.RL_TS .* 400);
                %         data.LeftReinfTS = offsetP(1)*data.I + offsetP(2);
                %         data.RightReinfTS = offsetP(1)*data.E + offsetP(2);
                %         data.BeamBreakTS = offsetP(1)*data.K + offsetP(2);
                %         idx = ismember(round(data.LL_TS*100)/100,round(data.LeftReinfTS*100)/100) ;
                %         data.LeftUnRefTS = data.LL_TS(~idx);
                %         idx = ismember(round(data.RL_TS*100)/100,round(data.RightReinfTS*100)/100) ;
                %         data.RightUnRefTS = data.RL_TS(~idx);
                
                %data.LL_Rf_IDX = round(LeftRefTS .* 400);
                %data.RL_Rf_IDX = round(RightRefTS .* 400);
            end
            if x==1
                data.BeamBreakTS =  data.MedPC.offsetP(1)*data.MedPC.(BBVar) +  data.MedPC.offsetP(2);
                data.Version = '04';
            end
        end
        %Lapish function not working at all
        function LapishDataType()
            % not sure if any of this works with onesided data
            % Set the 5 event types
            %        behavEvents_Chris = cell(0);
            %        tempCell = cell(0);
            %m=1;
            %n=1;
            %b=1;
            
            %lBEvt = cell(0);
            %lBTm = [];
            if ~strcmp(data.SessionType,'ContDeg')
                Evt{1}(data.([curLever '_Rf_Mask'])==1,1)={[curLever '_R']};
                Evt{1}(data.([curLever '_Rf_Mask'])==0,1)={[curLever '_U']};
                EvtTm{1}=data.([curLever '_TS']);
                if x == 1
                    Evt{2}(1:length(data.BeamBreakTS),1)={'BB'};
                    EvtTm{2}=[data.BeamBreakTS];
                end
            else
                Evt{1}(1:length(data.([curLever '_TS'])),1)={[curLever '_U']};
                EvtTm{1}=[data.([curLever '_TS'])];
                Evt{2}(1:length(data.([curLever '_Rf_TS'])),1)={[curLever '_U']};
                EvtTm{2}=[data.([curLever '_Rf_TS'])];
                if x == 1
                    Evt{3}(1:length(data.BeamBreakTS),1)={'BB'};
                    EvtTm{3}=[data.BeamBreakTS];
                end
            end
            Evt=cat(1,Evt{:});
            EvtTm=cat(1,EvtTm{:});
            data.Lapish.behaveEvt_Raw=vertcat(data.Lapish.behaveEvt_Raw, Evt);
            data.Lapish.behaveEvtTm_Raw = vertcat(data.Lapish.behaveEvtTm_Raw, EvtTm);
            [data.Lapish.behaveEvtTm_Raw,k] = sort(data.Lapish.behaveEvtTm_Raw);
            data.Lapish.behaveEvt_Raw = data.Lapish.behaveEvt_Raw(k);
            %             if ~strcmp(data.SessionType,'ContDeg')
            %                 Evt=cell(length(data.([curLever '_TS'])),1);
            %                 EvtTm=data.([curLever '_TS']);
            %                 for i = 1:length(data.([curLever '_TS']))
            %                     if data.([curLever '_Rf_Mask'])(i) == 0
            %                         Evt{m} = [curLever '_U'];
            %                         m = m+1;
            %                     elseif data.([curLever '_Rf_Mask'])(i) == 1
            %                         Evt{m} = [curLever '_R'];
            %                         m = m+1;
            %                         %       lfind(data.BeamBreakTS>data.LL_TS(i),1);
            %                     end
            %                 end
            %             else
            %                 Evt=cell(length(data.([curLever '_TS']))+length(data.([curLever 'Rf_TS'])),1  ,1);
            %                 EvtTm=[data.([curLever '_TS']) data.([curLever 'Rf_TS'])];
            %                 for i = 1:length(data.([curLever '_TS']))
            %                     Evt{m} = [curLever '_U'];
            %                     m = m+1;
            %                 end
            %                 for i = 1:length(data.([curLever '_Rf_TS']))
            %                     Evt{m} = [curLever '_R'];
            %                     m = m+1;
            %                 end
            %             end
            %             [EvtTm,k] = sort(EvtTm);
            %             Evt = Evt(k);
            %             data.Lapish.behaveEvt_Raw=[data.Lapish.behaveEvt_Raw Evt];
            %             data.Lapish.behaveEvtTm_Raw = [data.Lapish.behaveEvtTm_Raw EvtTm];
        end
    end
    function IDClippingArtifact()
        %data.goodDetectionMask = ones(,1);
        %% create filter for Clipping ID
        Fs = 400;  % Sampling data.Frequency
        Fstop1 = 2;       % First Stopband data.Frequency
        Fpass1 = 5;       % First Passband data.Frequency
        Fpass2 = 100;     % Second Passband data.Frequency
        Fstop2 = 110;     % Second Stopband data.Frequency
        Astop1 = 60;      % First Stopband Attenuation (dB)
        Apass  = 1;       % Passband Ripple (dB)
        Astop2 = 80;      % Second Stopband Attenuation (dB)
        match  = 'both';  % Band to match exactly
        % Construct an FDESIGN object and call its ELLIP method.
        h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
            Astop2, Fs);
        Hd = design(h, 'ellip', 'MatchExactly', match);
        G = Hd.ScaleValues;
        SOS = Hd.sosmatrix;
        
        %% Implement mtspecgramc from chronux
        params.Fs = 400;
        params.tapers = [2 1];
        params.fpass = [1 100];
        %data2run=mean(bsxfun(@minus,data.downDataRaw,downsample(data.refChan,100)),2);
        data2run=data.downDataRaw(:,6)-downsample(data.refChan,100);
        [S, t, f] = mtspecgramc(data2run,[.05 .01],params);
        plot_matrix(S,t,f);
        % still unclear if this was CAR or not, and stupid to do on just
        % one chan, wtf was i thinking.
        %[S, ~, ~] = mtspecgramc(data.Field(:,6),[.05 .01],params);
        %% First pass Identification
        logS = 10*log10(S);
        goodDetectionMask = ~((logS(:,4)) <  (mean(logS(:,4))- std(logS(:,4))*2.5) |  (logS(:,4)) >  (mean(logS(:,4))+ std(logS(:,4))*2));
        badIDXs = find(~goodDetectionMask);
        %conversion = round (length(data.Field(:,6))/(length(goodDetectionMask)));
        %% Use nature of artifact to improve ID
        method = 0;
        lastIDX = length(goodDetectionMask);
        for x = 1:length(badIDXs)
            switch method
                case  0 %first detection
                    if x +2 < length(badIDXs) %not the last IDX in badIDXs
                        if badIDXs(x+3) <= 10+badIDXs(x) %the next 2 bad idx are nearby
                            method = 1;
                            goodDetectionMask(badIDXs(x):badIDXs(x+3)) = 0; %not sure if this is the right idx
                            if badIDXs(x)-10 >= 1
                                goodDetectionMask(badIDXs(x)-10:badIDXs(x)) = 0;
                            else
                                goodDetectionMask(1:badIDXs(x)) = 0;
                                %the next idx is not nearby
                                goodDetectionMask(badIDXs(x)) = 1; %not sure if this is the right idx
                            end
                            x = x+1;
                        else
                            %is the last value in bad IDX
                            if badIDXs(x) + 30 < lastIDX
                                goodDetectionMask(badIDXs(x)) = 1; %not sure if this is the right idx
                            end
                        end
                    end
                case 1
                    if x < length(badIDXs) %not the last IDX in badIDXs
                        if badIDXs(x+1) <= badIDXs(x) + 189 %the next bad idx is in the next .3 seconds
                            goodDetectionMask(badIDXs(x):badIDXs(x+1)) = 0; %not sure if this is the right idx
                        else
                            method = 0;
                            %the next idx is not nearby, set the next .3 to 0, the
                            %rest to zero
                            if badIDXs(x) + 10 <= lastIDX
                                goodDetectionMask(badIDXs(x):badIDXs(x) +10) = 0; %not sure if this is the right idx
                                
                            else
                                goodDetectionMask(badIDXs(x):lastIDX) = 0;
                            end
                        end
                    else
                        %is the last value in bad IDX
                        if badIDXs(x) + 10 <= lastIDX
                            goodDetectionMask(badIDXs(x):badIDXs(x) +10) = 0; %not sure if this is the right idx
                            
                        else
                            goodDetectionMask(badIDXs(x):lastIDX) = 0;
                        end
                    end
            end
        end
        %% Creat Final Variable
        first = goodDetectionMask(1);
        last = goodDetectionMask(end);
        goodDetectionMask = double([first; first; goodDetectionMask; last; last]);
        data.goodDetectionMask = logical(round(interp (goodDetectionMask,4,1,0.5)));
    end
%% Write Binary for Spike Sorter
    function WriteBinary()
        % Load the filter
        SOS=[]; %wtf?
        G=[]; %wtf?
        load ('SpikesFilter01.mat');
        scalefactor = (10^8);
        %% Open and Read template ofi
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
        %% need to reload the raw complete data
        load(newMCDRawName);
        %out of place but needed
        chandata=bsxfun(@minus,chandata,data.refChan);
        %% filter data for spikes
        %load ('SpikesFilter01.mat');
        for x = 1:size(chandata,2)
            chandata(:,x) = filtfilt(SOS,G,chandata(:,x));
        end
        %% upsample goodDetectionMask
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
        %% these onevals aren't being used below and that makes no sense
        onevals = find(diff(mask) == +1)+1;
        for x = 1:length(onevals)
            mask(onevals(x):onevals(x)+(.05*400)) = 0;
        end
        maskup = upsample(mask,100);
        for x = 1:data.SessionLength*60*data.Freq/100
            maskup(x*100-99:x*100) = mask (x);
        end
        %% set bad data to 0
        chandata (~maskup,:) = 0;
        chandata = int16(chandata*scalefactor).';
        %% write ofi
        fid = fopen(newOfiName, 'w');
        A(1,2) = cellstr(sprintf('MaxMV=%.6f',(2^15-1)*scalefactor*1000));
        A(1,3) = cellstr([ 'NChan=' num2str(length(data.ChannelList))]);
        for i = 1:numel(A)
            if A{i+1} == -1
                fprintf(fid,'%s', A{i});
                break
            else
                fprintf(fid,'%s\n', A{i});
            end
        end
        fclose (fid);
        %% write bin
        fid = fopen(newBinName, 'w');
        fwrite(fid, chandata,'int16');
        fclose (fid);
    end


% -----------UnUsed for Spike Sorting-------------------------------------
%% Preprocessing
%write bin and the other file for reading into spike sorter
%filtering etc.
%% Spectrum code, currently unused
%PreProcessFields();
%LeverPressesSpectrum();
%OverallSpectrum();
% ------------------------------------------------------------------------

%% Lever Press Basic Spectrum (and some of the detailed spectrum)
    function PreProcessFields()
        %data.Field = downsample(data.chandata,downRate);
        %data=rmfield (data, 'chandata');
        
        data.nacidx = (data.ChannelList == 2 | data.ChannelList == 3);
        data.pfcidx = ~(data.nacidx);
        Fs = 40000;  % Sampling data.Frequency
        
        Fstop1 = 0.01;    % First Stopband data.Frequency
        Fpass1 = 0.1;     % First Passband data.Frequency
        Fpass2 = 200;     % Second Passband data.Frequency
        Fstop2 = 250;     % Second Stopband data.Frequency
        Astop1 = 60;      % First Stopband Attenuation (dB)
        Apass  = 1;       % Passband Ripple (dB)
        Astop2 = 80;      % Second Stopband Attenuation (dB)
        match  = 'both';  % Band to match exactly
        
        % Construct an FDESIGN object and call its ELLIP method.
        [m, n] = size(data.chandata);
        scaleFactor= 100;
        h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
            Astop2, Fs);
        Hd = design(h, 'ellip', 'MatchExactly', match);
        tempRefed = zeros([m/scaleFactor n]);
        tempUnRefed = zeros([m/scaleFactor n]);
        
        %% needs vectorization - improve speed???
        for x = 1:n
            temp = data.chandata(:,x) - data.refChan;
            temp2 = filter(Hd, temp);
            tempRefed(:,x) = downsample (temp2,scaleFactor);
            
            temp = data.chandata(:,x);
            temp2 = filter(Hd, temp);
            tempUnRefed(:,x) = downsample (temp2,scaleFactor);
        end
        %%
        data.Field = tempRefed;
        data.FieldUnRefed = tempUnRefed;
        
        data=rmfield (data, 'chandata');
        data = rmfield (data, 'refChan');
        %        save(sprintf ('%sfields_%s_%s.mat', location,versionNumber,filename ),'data');
        %       clearvars ('-except', 'dateArray');
        %         data.LL_IDX = data.LL_IDX./100;
        %         data.RL_IDX = data.RL_IDX./100;
        %         data.LL_Rf_IDX = data.LL_Rf_IDX./100;
        %         data.RL_Rf_IDX = data.RL_Rf_IDX./100;
        
        %% Find bad Events
        data.goodDetectionMask = ones(720000,1);
        %if isfield(data, 'refChanDown')
        %rawdata = (data.Field(:,6))+ data.refChanDown;
        %else
        Fs = 400;  % Sampling data.Frequency
        
        Fstop1 = 2;       % First Stopband data.Frequency
        Fpass1 = 5;       % First Passband data.Frequency
        Fpass2 = 100;     % Second Passband data.Frequency
        Fstop2 = 110;     % Second Stopband data.Frequency
        Astop1 = 60;      % First Stopband Attenuation (dB)
        Apass  = 1;       % Passband Ripple (dB)
        Astop2 = 80;      % Second Stopband Attenuation (dB)
        match  = 'both';  % Band to match exactly
        
        % Construct an FDESIGN object and call its ELLIP method.
        h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
            Astop2, Fs);
        Hd = design(h, 'ellip', 'MatchExactly', match);
        G = Hd.ScaleValues;
        SOS = Hd.sosmatrix;
        params.Fs = 400;
        params.tapers = [2 1];
        params.fpass = [1 100];
        [S t f] = mtspecgramc(data.Field(:,6),[.05 .01],params);
        logS = 10*log10(S);
        goodDetectionMask = ~((logS(:,4)) <  (mean(logS(:,4))- std(logS(:,4))*2.5) |  (logS(:,4)) >  (mean(logS(:,4))+ std(logS(:,4))*2));
        badIDXs = find(~goodDetectionMask);
        conversion = round (length(data.Field(:,6))/(length(goodDetectionMask)));
        method = 0;
        lastIDX = length(goodDetectionMask);
        %% Needs vectorization -- improve speed???
        for x = 1:length(badIDXs)
            switch method
                case  0 %first detection
                    if x +2 < length(badIDXs) %not the last IDX in badIDXs
                        if badIDXs(x+3) <= 10+badIDXs(x) %the next 2 bad idx are nearby
                            method = 1;
                            goodDetectionMask(badIDXs(x):badIDXs(x+3)) = 0; %not sure if this is the right idx
                            if badIDXs(x)-10 >= 1
                                goodDetectionMask(badIDXs(x)-10:badIDXs(x)) = 0;
                            else
                                goodDetectionMask(1:badIDXs(x)) = 0;
                                %the next idx is not nearby
                                goodDetectionMask(badIDXs(x)) = 1; %not sure if this is the right idx
                            end
                            x = x+1;
                        else
                            %is the last value in bad IDX
                            if badIDXs(x) + 30 < lastIDX
                                goodDetectionMask(badIDXs(x)) = 1; %not sure if this is the right idx
                            end
                        end
                    end
                case 1
                    if x < length(badIDXs) %not the last IDX in badIDXs
                        if badIDXs(x+1) <= badIDXs(x) + 189 %the next bad idx is in the next .3 seconds
                            goodDetectionMask(badIDXs(x):badIDXs(x+1)) = 0; %not sure if this is the right idx
                        else
                            method = 0;
                            %the next idx is not nearby, set the next .3 to 0, the
                            %rest to zero
                            if badIDXs(x) + 10 <= lastIDX
                                goodDetectionMask(badIDXs(x):badIDXs(x) +10) = 0; %not sure if this is the right idx
                                
                            else
                                goodDetectionMask(badIDXs(x):lastIDX) = 0;
                            end
                        end
                    else
                        %is the last value in bad IDX
                        if badIDXs(x) + 10 <= lastIDX
                            goodDetectionMask(badIDXs(x):badIDXs(x) +10) = 0; %not sure if this is the right idx
                            
                        else
                            goodDetectionMask(badIDXs(x):lastIDX) = 0;
                        end
                    end
            end
        end
        %%
        first = goodDetectionMask(1);
        last = goodDetectionMask(end);
        goodDetectionMask = double([first; first; goodDetectionMask; last; last]);
        data.goodDetectionMask = logical(round(interp (goodDetectionMask,4,1,0.5)));
    end
    function OverallSpectrum()
        
        params.tapers = [5 9];
        params.Fs = 400;
        if data.FirstArray == 'LL'
            LLrange = logical([data.goodDetectionMask(1:360000); zeros(360000,1)]); RLrange = logical( [zeros(360000,1) ; data.goodDetectionMask(360001:720000) ]);
        else RLrange = logical([zeros(360000,1) ; data.goodDetectionMask(360001:720000) ]); LLrange = logical([data.goodDetectionMask(1:360000); zeros(360000,1)]);
        end
        %% Needs vectorization - improve speed, remove errors
        for y = 1:size(data.ChannelList,2) %x = channel
            %for x=4
            %for y = 1:size(tempDataRefed{fnCount},3) %y = event
            %        %for z = 1:6
            %         %   params.fpass = data.Freqs{z};
            if data.pfcidx (y)
                %normval = bandpower (data.Field(find(data.goodDetectionMask),y), Fs, data.Freqs{z});
                %                overallPFCRI {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RIrange)),y), Fs, data.Freqs{z});%/normval;
                %               overallPFCRR {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RRrange)),y), Fs, data.Freqs{z});%/normval;
                %                %[S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RIrange)),y), params);%/normval;
                %                %overallPFCRI {x} (y,z) = mean(S);
                %                %[S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RRrange)),y), params);%/normval;
                %                %overallPFCRR {x} (y,z) = mean(S);
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask,y), params);%/normval;
                NormPFC (y,:) = S./data.SimpleNormArray(1,1,y);
                NormPFCdata.Freq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(LLrange),y), params);%/normval;
                data.overallPFCLL  (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallPFCLLdata.Freq  (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(RLrange),y), params);%/normval;
                data.overallPFCRL (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallPFCRLdata.Freq (y,:) = f;
            else
                %normval = bandpower (data.Field(find(data.goodDetectionMask),y), Fs, data.Freqs{z});
                %               overallNACRI {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RIrange)),y), Fs, data.Freqs{z});%/normval;
                %               overallNACRR {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RRrange)),y), Fs, data.Freqs{z});%//normval;
                %                [S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RIrange)),y), params);%/normval;
                %                overallNACRI {x} (y,z) = mean(S);
                %                [S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RRrange)),y), params);%/normval;
                %                overallNACRR {x} (y,z) = mean(S);
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(),y), params);%/normval;
                NormNAC (y,:) = S./data.SimpleNormArray(1,1,y);
                NormNACdata.Freq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(LLrange),y), params);%/normval;
                data.overallNACLL (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallNACLLdata.Freq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(RLrange),y), params);%/normval;
                data.overallNACRL  (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallNACRLdata.Freq (y,:) = f;
            end
            %tempChanDataBasicNAC {fnCount} (x,y,z) = bandpower (tempDataRefed{fnCount}(:,nacidx,y), Fs, data.Freqs{z});
            %[S, t, f] = mtspecgramc (tempDataRefed{fnCount}(:,x,y), [1 .1], params);
            %tempChanDataBasic{fnCount} (:,:,y,x) = S;
        end
    end
    function LeverPressesSpectrum()
        % Several changes need to bee implemented to the press spectrum
        % stuff.
        % 1) I need to take longer chunks, so I can drop off parts (handle
        % the tails of the spectrum
        % 2) I need to vectorize
        % 3) I need to allow for different kinds of normalization to do to
        % this is probably
        % 4) I need to allow for different was to select which events to
        % use
        % The smart way to do this is probably to pull out ALL of the
        % events data, in a wide range, and then save it into either a) its
        % own file b) the main struct, c) or the raw data struct
        % There is the small issue of deciding how to handle different
        % versions (ie. filters, referencing, denoising)
        
        % Perhaos most important, I need to make this more modular
        
        
        timeRange = 2;
        data.LL_TSGood = logical(zeros(size(data.LL_TS,1),1));
        data.RL_TSGood = logical(zeros(size(data.RL_TS,1),1));
        %%
        % Need to vectorize, pretty simple, need a fx to get the range,
        % and arrayfun tp run it. not sure if it will be faster or slower
        
        
        % the final goal is to achieve pull out the submatrix from the
        % overall matrix.
        eventCount = 1;
        for eventIDX = 1:size(data.LL_TS,1)
            if data.LL_TS(eventIDX) >= timeRange && data.LL_TS(eventIDX) <=1797
                dataRangeIDX =[round((data.LL_TS(eventIDX) -2) *400): round(((data.LL_TS(eventIDX) +2) *400-1))];
                leverPreRange = dataRangeIDX (1,1:800);
                leverPostRange = dataRangeIDX (1,801:1600);
                if min(data.goodDetectionMask(dataRangeIDX,1))
                    data.LL_TSGood(eventIDX) = 1;
                    %temp = data.Field([round((data.LL_TS(eventIDX) -2) *400): round(((data.LL_TS(eventIDX) +.5) *400-1))],:);
                    data.LLPre_DataChunks(eventCount,:,:) =  data.Field(leverPreRange,:);
                    data.LLPost_DataChunks(eventCount,:,:) =  data.Field(leverPostRange,:);
                    eventCount = eventCount + 1;
                end
            end
        end
        params.tapers = [5 9];
        params.Fs = 400;
        params.trialave = 0;
        
        [S,f]= mtspectrumc (squeeze(data.LLPre_DataChunks(1,:,:)), params);
        data.LLPre_SimpleSpectrum = zeros(size(data.LLPre_DataChunks,1),size(S,1),size(S,2));
        data.LLPre_SimpleSpectrum(1,:,:) = S;
        data.LLPre_SimpleF = f;
        [S,f]= mtspectrumc (squeeze(data.LLPost_DataChunks(1,:,:)), params);
        data.LLPost_SimpleSpectrum = zeros(size(data.LLPost_DataChunks,1),size(S,1),size(S,2));
        data.LLPost_SimpleSpectrum(1,:,:) = S;
        data.LLPost_SimpleF = f;
        [S, t, f] = mtspecgramc (  squeeze([data.LLPre_DataChunks(1,:,:) data.LLPost_DataChunks(1,:,:)]) , [1  .1], params);
        data.LLAll_ComplexSpectrum = zeros(size(data.LLPre_DataChunks,1),size(S,1),size(S,2), size(S,3));
        data.LLAll_ComplexSpectrum (1,:,:,:) = S;
        data.LLAll_ComplexF = f;
        data.LLAll_ComplexT = t;
        for goodEventIDX = 2:size(data.LLPre_DataChunks,1)
            [S,f]= mtspectrumc (squeeze( data.LLPre_DataChunks(goodEventIDX,:,:)), params);
            data.LLPre_SimpleSpectrum(goodEventIDX,:,:) = S;
            [S,f]= mtspectrumc (squeeze( data.LLPost_DataChunks(goodEventIDX,:,:)), params);
            data.LLPost_SimpleSpectrum(goodEventIDX,:,:) = S;
            [S, t, f] = mtspecgramc (  squeeze([data.LLPre_DataChunks(goodEventIDX,:,:) data.LLPost_DataChunks(goodEventIDX,:,:)]) , [1  .1], params);
            data.LLAll_ComplexSpectrum (goodEventIDX,:,:,:) = S;
        end
        
        
        eventCount = 1;
        for eventIDX = 1:size(data.RL_TS,1)
            if data.RL_TS(eventIDX) >= timeRange && data.RL_TS(eventIDX) <=1797
                dataRangeIDX =[round((data.RL_TS(eventIDX) -2) *400): round(((data.RL_TS(eventIDX) +2) *400-1))];
                leverPreRange = dataRangeIDX (1,1:800);
                leverPostRange = dataRangeIDX (1,801:1600);
                if min(data.goodDetectionMask(dataRangeIDX,1))
                    data.RL_TSGood(eventIDX) = 1;
                    %temp = data.Field([round((data.RL_TS(eventIDX) -2) *400): round(((data.RL_TS(eventIDX) +.5) *400-1))],:);
                    data.RLPre_DataChunks(eventCount,:,:) =  data.Field(leverPreRange,:);
                    data.RLPost_DataChunks(eventCount,:,:) =  data.Field(leverPostRange,:);
                    eventCount = eventCount + 1;
                end
            end
        end
        
        [S,f]= mtspectrumc (squeeze(data.RLPre_DataChunks(1,:,:)), params);
        data.RLPre_SimpleSpectrum = zeros(size(data.RLPre_DataChunks,1),size(S,1),size(S,2));
        data.RLPre_SimpleSpectrum(1,:,:) = S;
        data.RLPre_SimpleF = f;
        [S,f]= mtspectrumc (squeeze(data.RLPost_DataChunks(1,:,:)), params);
        data.RLPost_SimpleSpectrum = zeros(size(data.RLPost_DataChunks,1),size(S,1),size(S,2));
        data.RLPost_SimpleSpectrum(1,:,:) = S;
        data.RLPost_SimpleF = f;
        [S, t, f] = mtspecgramc (  squeeze([data.RLPre_DataChunks(1,:,:) data.RLPost_DataChunks(1,:,:)]) , [1  .1], params);
        data.RLAll_ComplexSpectrum = zeros(size(data.RLPre_DataChunks,1),size(S,1),size(S,2), size(S,3));
        data.RLAll_ComplexSpectrum (1,:,:,:) = S;
        data.RLAll_ComplexF = f;
        data.RLAll_ComplexT = t;
        for goodEventIDX = 2:size(data.RLPre_DataChunks,1)
            [S,f]= mtspectrumc (squeeze( data.RLPre_DataChunks(goodEventIDX,:,:)), params);
            data.RLPre_SimpleSpectrum(goodEventIDX,:,:) = S;
            [S,f]= mtspectrumc (squeeze( data.RLPost_DataChunks(goodEventIDX,:,:)), params);
            data.RLPost_SimpleSpectrum(goodEventIDX,:,:) = S;
            [S, t, f] = mtspecgramc (  squeeze([data.RLPre_DataChunks(goodEventIDX,:,:) data.RLPost_DataChunks(goodEventIDX,:,:)]) , [1  .1], params);
            data.RLAll_ComplexSpectrum (goodEventIDX,:,:,:) = S;
        end
        
        data.NacIDX = (data.ChannelList == 2 | data.ChannelList == 3);
        data.PfcIDX = ~(data.NacIDX);
        
        [S f] = mtspectrumc(data.Field(logical(data.goodDetectionMask),:),params);
        data.SimpleNormArray (1,1,:) = mean(S).';
        data.ComplexNormArray (1,1,1,:) = mean(S).';
        data.LLPre_SNorm = bsxfun(@rdivide ,data.LLPre_SimpleSpectrum, data.SimpleNormArray);
        data.LLPost_SNorm = bsxfun(@rdivide ,data.LLPost_SimpleSpectrum, data.SimpleNormArray);
        data.LLAll_CNorm = bsxfun(@rdivide ,data.LLAll_ComplexSpectrum, data.ComplexNormArray);
        %old old         data.RIRefPFCSimpleRaw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.PfcIDX);
        %new old        data.LLPre_Rf_PFC_Simple_Raw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.PfcIDX);
        data.LLPre_Rf_PFC_Simple_Raw = data.LLPre_SNorm(data.LL_Rf_Mask(data.LL_TSGood),:,data.PfcIDX);
        data.LLPost_Rf_PFC_Simple_Raw = data.LLPost_SNorm(data.LL_Rf_Mask(data.LL_TSGood),:,data.PfcIDX);
        data.LLAll_Rf_PFC_Complex_Raw = data.LLAll_CNorm(data.LL_Rf_Mask(data.LL_TSGood),:,:,data.PfcIDX);
        %         data.RIUnRefPFCSimpleRaw = data.LLSNorm(find(data.RightUnRefIDX(data.RL_TSGood)),:,data.PfcIDX);
        data.LLPre_URf_PFC_Simple_Raw = data.LLPre_SNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,data.PfcIDX);
        data.LLPost_URf_PFC_Simple_Raw = data.LLPost_SNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,data.PfcIDX);
        data.LLAll_URf_PFC_Complex_Raw = data.LLAll_CNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,:,data.PfcIDX);
        %         data.RIRefNACSimpleRaw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.NacIDX);
        data.LLPre_Rf_NAC_Simple_Raw = data.LLPre_SNorm(data.LL_Rf_Mask(data.LL_TSGood),:,data.NacIDX);
        data.LLPost_Rf_NAC_Simple_Raw = data.LLPost_SNorm(data.LL_Rf_Mask(data.LL_TSGood),:,data.NacIDX);
        data.LLAll_Rf_NAC_Complex_Raw = data.LLAll_CNorm(data.LL_Rf_Mask(data.LL_TSGood),:,:,data.NacIDX);
        %         data.RIUnRefNACSimpleRaw = data.LLSNorm(find(data.RightUnRefIDX(data.RL_TSGood)),:,data.NacIDX);
        data.LLPre_URf_NAC_Simple_Raw = data.LLPre_SNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,data.NacIDX);
        data.LLPost_URf_NAC_Simple_Raw = data.LLPost_SNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,data.NacIDX);
        data.LLAll_URf_NAC_Complex_Raw = data.LLAll_CNorm(~data.LL_Rf_Mask(data.LL_TSGood),:,:,data.NacIDX);
        
        
        
        
        data.RLPre_SNorm = bsxfun(@rdivide ,data.RLPre_SimpleSpectrum, data.SimpleNormArray);
        data.RLPost_SNorm = bsxfun(@rdivide ,data.RLPost_SimpleSpectrum, data.SimpleNormArray);
        data.RLAll_CNorm = bsxfun(@rdivide ,data.RLAll_ComplexSpectrum, data.ComplexNormArray);
        %old old         data.RIRefPFCSimpleRaw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.PfcIDX);
        %new old        data.LLPre_Rf_PFC_Simple_Raw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.PfcIDX);
        data.RLPre_Rf_PFC_Simple_Raw = data.RLPre_SNorm(data.RL_Rf_Mask(data.RL_TSGood),:,data.PfcIDX);
        data.RLPost_Rf_PFC_Simple_Raw = data.RLPost_SNorm(data.RL_Rf_Mask(data.RL_TSGood),:,data.PfcIDX);
        data.RLAll_Rf_PFC_Complex_Raw = data.RLAll_CNorm(data.RL_Rf_Mask(data.RL_TSGood),:,:,data.PfcIDX);
        %         data.RIUnRefPFCSimpleRaw = data.LLSNorm(find(data.RightUnRefIDX(data.RL_TSGood)),:,data.PfcIDX);
        data.RLPre_URf_PFC_Simple_Raw = data.RLPre_SNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,data.PfcIDX);
        data.RLPost_URf_PFC_Simple_Raw = data.RLPost_SNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,data.PfcIDX);
        data.RLAll_URf_PFC_Complex_Raw = data.RLAll_CNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,:,data.PfcIDX);
        %         data.RIRefNACSimpleRaw = data.LLSNorm(find(data.RL_Rf_Mask(data.RL_TSGood)),:,data.NacIDX);
        data.RLPre_Rf_NAC_Simple_Raw = data.RLPre_SNorm(data.RL_Rf_Mask(data.RL_TSGood),:,data.NacIDX);
        data.RLPost_Rf_NAC_Simple_Raw = data.RLPost_SNorm(data.RL_Rf_Mask(data.RL_TSGood),:,data.NacIDX);
        data.RLAll_Rf_NAC_Complex_Raw = data.RLAll_CNorm(data.RL_Rf_Mask(data.RL_TSGood),:,:,data.NacIDX);
        %         data.RIUnRefNACSimpleRaw = data.LLSNorm(find(data.RightUnRefIDX(data.RL_TSGood)),:,data.NacIDX);
        data.RLPre_URf_NAC_Simple_Raw = data.RLPre_SNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,data.NacIDX);
        data.RLPost_URf_NAC_Simple_Raw = data.RLPost_SNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,data.NacIDX);
        data.RLAll_URf_NAC_Complex_Raw = data.RLAll_CNorm(~data.RL_Rf_Mask(data.RL_TSGood),:,:,data.NacIDX);
        %data.NacIDX = (data.ChannelList == 2 | data.ChannelList == 3);
        %data.PfcIDX = ~(data.NacIDX);
        %params.tapers = [5 9];
        %params.Fs = 400;
        % [S f] = mtspectrumc(data.Field(logical(data.goodDetectionMask),:),params);
        
        %
        %save ('compiled.mat', 'NormPFC', 'fileList', 'NormPFCdata.Freq', 'overallPFCLL', 'overallPFCLLdata.Freq', 'overallPFCRL', 'overallPFCRLdata.Freq', 'NormNAC', 'NormNACdata.Freq', 'overallNACLL', 'overallNACLLdata.Freq', 'overallNACRL', 'overallNACRLdata.Freq');
    end
end

