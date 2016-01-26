function [ data  ] = readRawFiles( badChannels, varargin )
%READRAWFILES Reads in .mcd and MedPC and returns data struct
%   Detailed explanation goes here

if nargin >= 3
    filenameMCD = varargin{1};
    filenameMedPC = varargin{2};
    if nargin == 4 && strcmp (varargin{3}, 'raw')
        method='rawOnly';
    end
    
elseif nargin == 1
    [filename,pathname] = uigetfile('*.mcd');
    filenameMCD = [ pathname filename];
    [filename,pathname] = uigetfile('*.txt');
    filenameMedPC = [ pathname filename];
else
    display ('You did not provide a correct number of inputs');
    
    
    
end

%this should end up with a data file which looks like this...
%{
data.Animal
data.Task
data.FirstLever
data.StartDate
data.EndDate
data.StartTime
data.EndTime
data.L_Press
data.L_Ref
data.R_Press
data.R_Ref
data.MagEntry
data.Nex
data.Nex. (there should be a template somewhere.
%}

%Get FileNames...
%MCD Read
%MedPC Read
Hd = load('strong bandpass coeff.mat');
channelList = [1:1:16];

data.ChannelList =  channelList(~ismember(channelList, badChannels));
goodChannels = data.ChannelList;
data.MCD = filenameMCD;
data.MedPCFile = filenameMedPC;
readMedPC();
%Merge MCD and Med
readMCD();
MCDMedPCSync();
PreProcessFields();

%LeverPressesSpectrum();
%OverallSpectrum();
save([filenameMCD(1:end-4) '_v01_6_01.mat' ], 'data');


%Preprocessing
%Saving to Nex
%filter data
%These are nested functions so they can get variables from the parent
    function readMedPC ()
        %UNTITLED3 Summary of this function goes here
        %   This would be much smarter if it looked for the name of the
        %   variable, instead of having it predefined.
        
        %current version only accepts a single file name
        fid = fopen(filenameMedPC);
        
        %  val = isspace(test); gets logical array for spaces. if first char is a
        %  space, then the line is a continuation of previous variable
        
        % need to make a struct for each file that contains...
        % Start Time (3)
        % End Time (3)
        % File run (MSN) (1) i think
        % Single variables (A B D F and sometimes J(licks))
        % skip to arrays (C
        
        %Open File
        
        
        %%%%%%%%%%%%start parsing
        
        %%%Text Variables
        %Start Date
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
                    data.MSN = strtrim(curLineSplit{2});
                    success = 1;
                    %%% Extract here the actual protocol for future deviations%%%
                end
            end
        end
        
        %%%Char Variables
        %A (int) # of Left Reinforcers
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'A'));
                    data.MedPC.A = str2double(curLineSplit{2});
                    success = 1;
                end
            end
        end
        %B
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'B'));
                    data.MedPC.B = str2double(curLineSplit{2});
                    success = 1;
                end
            end
        end
        %D
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'D'));
                    data.MedPC.D = str2double(curLineSplit{2});
                    success = 1;
                end
            end
        end
        
        %F
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'F'));
                    data.MedPC.F = str2double(curLineSplit{2});
                    success = 1;
                end
            end
        end
        
        %C first ARRAY!!!!!!
        success = 0;
        stop = 0;
        while ~stop && ~success
            curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'C'));
                    tempdataArray = zeros(20000,1);
                    arrayX = 1;
                    while (~stop && ~success && arrayX < 20001)
                        curLine = fgetl(fid);
                        if curLine == -1
                            stop;
                        elseif curLine(1) == ' ';
                            curLineSplit = strsplit(curLine);  %white space this time
                            for x = 1:length(curLineSplit) - 2
                                tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                arrayX = arrayX +1;
                            end
                        else
                            data.MedPC.C = nonzeros (tempdataArray);
                            success = 1;
                        end
                    end
                end
            end
        end
        
        %E array
        success = 0;
        stop = 0;
        while ~stop && ~success
            %    curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'E'));
                    tempdataArray = zeros(20000,1);
                    arrayX = 1;
                    while (~stop && ~success && arrayX < 20001)
                        curLine = fgetl(fid);
                        if curLine == -1
                            stop;
                        elseif curLine(1) == ' ';
                            curLineSplit = strsplit(curLine);  %white space this time
                            for x = 1:length(curLineSplit) - 2
                                tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                arrayX = arrayX +1;
                            end
                        else
                            data.MedPC.E = nonzeros (tempdataArray);
                            success = 1;
                        end
                    end
                end
            end
        end
        
        %G array
        success = 0;
        stop = 0;
        while ~stop && ~success
            %    curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'G'));
                    tempdataArray = zeros(20000,1);
                    arrayX = 1;
                    while (~stop && ~success && arrayX < 20001)
                        curLine = fgetl(fid);
                        if curLine == -1
                            stop;
                        elseif curLine(1) == ' ';
                            curLineSplit = strsplit(curLine);  %white space this time
                            for x = 1:length(curLineSplit) - 2
                                tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                arrayX = arrayX +1;
                            end
                        else
                            data.MedPC.G = nonzeros (tempdataArray);
                            success = 1;
                        end
                    end
                end
            end
        end
        %I array
        success = 0;
        stop = 0;
        while ~stop && ~success
            %    curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'I'));
                    tempdataArray = zeros(20000,1);
                    arrayX = 1;
                    while (~stop && ~success && arrayX < 20001)
                        curLine = fgetl(fid);
                        if curLine == -1
                            stop;
                        elseif curLine(1) == ' ';
                            curLineSplit = strsplit(curLine);  %white space this time
                            for x = 1:length(curLineSplit) - 2
                                tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                arrayX = arrayX +1;
                            end
                        else
                            data.MedPC.I = nonzeros (tempdataArray);
                            success = 1;
                        end
                    end
                end
            end
        end
        %sometimes I K U and missing the time interval data...
        
        
        %K array
        success = 0;
        stop = 0;
        while ~stop && ~success
            %    curLine = fgetl(fid);
            if curLine == -1;
                stop = 1;
            else
                curLineSplit = strsplit(curLine, ':');
                if (strcmp (curLineSplit{1}, 'K'));
                    tempdataArray = zeros(20000,1);
                    arrayX = 1;
                    while (~stop && ~success && arrayX < 20001)
                        curLine = fgetl(fid);
                        if curLine == -1
                            stop;
                        elseif curLine(1) == ' ';
                            curLineSplit = strsplit(curLine);  %white space this time
                            for x = 1:length(curLineSplit) - 2
                                tempdataArray(arrayX,1) = str2double(curLineSplit{x+2});
                                arrayX = arrayX +1;
                            end
                        else
                            data.MedPC.K = nonzeros (tempdataArray);
                            success = 1;
                        end
                    end
                end
            end
        end
        
        
        
        
        
        
        %close file
        fclose(fid);
        %clear success arrayX curLine curLineSplot tempdataArray curLineSplit fid stop x;
        
    end
    function readMCD
        tempNex = open ('BasicNexTemplate02.mat');
        data.Nex = tempNex.nexFile;
        clear tempNex;
        
        [pathname, fileName, ext]=fileparts(which('nsMCDLibrary64.dll'));
        ns_SetLibrary([pathname filesep fileName ext]);
        [nsresult, hfile] = ns_OpenFile(filenameMCD);
        [nsresult, FileInfo] = ns_GetFileInfo(hfile);
        [nsresult, EntityInfo] = ns_GetEntityInfo(hfile, [1 : 1 : FileInfo.EntityCount]);
        freq = 1/FileInfo.TimeStampResolution;
        
        [nsresult,data.EventData,~,~] = ns_GetEventData(hfile,2,[1:2:EntityInfo(2).ItemCount]);
        startIDX = data.EventData(1) * 40000;
        
        data.EventData = data.EventData (2:length(data.EventData),:) - data.EventData(1);
        channelList = goodChannels;
        %reference acquisition
        x = 1;
        tempArray = zeros (72000000,1);
        sumArray = zeros(72000000,1);
        
        for chanNumber = channelList
            [nsresult,~,tempArray]= ns_GetAnalogData(hfile,chanNumber+17,startIDX,72000000);
            sumArray = sumArray + tempArray;
            %data.Nex = nexAddContinuous(data.Nex, 0, freq, tempdata, 'chanNameHere');
            x = x+1;
        end
        data.refChan = sumArray/length(channelList);
        clear ('sumArray');
        data.chandata = zeros (72000000, length(channelList));
        x=1;
        %         if (strcmp(method, 'rawOnly'))
        for chanNumber = channelList
            [nsresult,~,data.chandata(:,x)]= ns_GetAnalogData(hfile,chanNumber+17,startIDX,72000000);
            %data.Nex = nexAddContinuous(data.Nex, 0, freq, tempdata, 'chanNameHere');
            
            x = x+1;
        end
        %         else
        %             for chanNumber = channelList
        %                 [nsresult,~,tempArray]= ns_GetAnalogdata(hfile,chanNumber+17,startIDX,72000000);
        %                 %data.Nex = nexAddContinuous(data.Nex, 0, freq, tempdata, 'chanNameHere');
        %                 tempArray = tempArray - data.refChan;
        %                 data.chandata(:,x) = filtfilt(Hd.SOS,Hd.G, tempArray);
        %                 x = x+1;
        %             end
        %         end
        clear tempArray;
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
        
        %%
        %find out who is first, left lever (C) or right lever (G)
        %
        if data.MedPC.C(2) < data.MedPC.G(2); firstArray = data.MedPC.C; data.FirstArray = 'L';
        else firstArray = data.MedPC.G; data.FirstArray = 'R'; end;
        
        %calculate linear offset
        lastIDX = length(firstArray);
        data.MedPC.offsetP = polyfit (firstArray, data.EventData(1:lastIDX,1),1);
        
        if data.FirstArray == 'L'
            data.LL_TS = data.EventData(1:lastIDX);
            data.RL_TS = data.EventData(lastIDX+1:length(data.EventData),1);
        else
            data.LL_TS = data.EventData(lastIDX+1:length(data.EventData),1);
            data.RL_TS = data.EventData(1:lastIDX);
        end
        %         data.LeftReinfTS = offsetP(1)*data.I + offsetP(2);
        %         data.RightReinfTS = offsetP(1)*data.E + offsetP(2);
        %         data.BeamBreakTS = offsetP(1)*data.K + offsetP(2);
        %         idx = ismember(round(data.LL_TS*100)/100,round(data.LeftReinfTS*100)/100) ;
        %         data.LeftUnRefTS = data.LL_TS(~idx);
        %         idx = ismember(round(data.RL_TS*100)/100,round(data.RightReinfTS*100)/100) ;
        %         data.RightUnRefTS = data.RL_TS(~idx);
        
        LeftRefTS =  data.MedPC.offsetP(1)*data.MedPC.I +  data.MedPC.offsetP(2);
        RightRefTS =  data.MedPC.offsetP(1)*data.MedPC.E +  data.MedPC.offsetP(2);
        data.BeamBreakTS =  data.MedPC.offsetP(1)*data.MedPC.K +  data.MedPC.offsetP(2);
        data.LL_IDX = round(data.LL_TS .* 400);
        data.RL_IDX = round(data.RL_TS .* 400);
        
        %if these are done at 40000, they don't come out correctly with
        %RF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        %This datype is being retired. it was an index to the ts. its now
        %an idx to the data itself
        
        %data.LL_Rf_IDX = ismember(round(data.LL_TS*100)/100,round(LeftRefTS*100)/100) ;
        %data.LL_URf_IDX = ~data.LL_Rf_IDX;
        %data.RL_Rf_IDX = ismember(round(data.RL_TS*100)/100,round(RightRefTS*100)/100) ;
        %data.RL_URf_IDX = ~data.RL_Rf_IDX;
        
        data.LL_Rf_IDX = round(LeftRefTS .* 400);
        data.RL_Rf_IDX = round(RightRefTS .* 400);
        data.LL_Rf_Mask =ismember(round(data.LL_IDX) ,round(data.LL_Rf_IDX));
        data.RL_Rf_Mask=ismember(round(data.RL_IDX) ,round(data.RL_Rf_IDX));
        
        
        data.Version = 'bulk03';
        
        %clear lastIDX offsetP;
    end
    function PreProcessFields()
        %data.Field = downsample(data.chandata,downRate);
        %data=rmfield (data, 'chandata');
        
        data.nacidx = (data.ChannelList == 2 | data.ChannelList == 3);
        data.pfcidx = ~(data.nacidx);
        Fs = 40000;  % Sampling Frequency
        
        Fstop1 = 0.01;    % First Stopband Frequency
        Fpass1 = 0.1;     % First Passband Frequency
        Fpass2 = 200;     % Second Passband Frequency
        Fstop2 = 250;     % Second Stopband Frequency
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
        Fs = 400;  % Sampling Frequency
        
        Fstop1 = 2;       % First Stopband Frequency
        Fpass1 = 5;       % First Passband Frequency
        Fpass2 = 100;     % Second Passband Frequency
        Fstop2 = 110;     % Second Stopband Frequency
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
        if data.FirstArray == 'L'
            LLrange = logical([data.goodDetectionMask(1:360000); zeros(360000,1)]); RLrange = logical( [zeros(360000,1) ; data.goodDetectionMask(360001:720000) ]);
        else RLrange = logical([zeros(360000,1) ; data.goodDetectionMask(360001:720000) ]); LLrange = logical([data.goodDetectionMask(1:360000); zeros(360000,1)]);
        end
        %% Needs vectorization - improve speed, remove errors
        for y = 1:size(data.ChannelList,2) %x = channel
            %for x=4
            %for y = 1:size(tempDataRefed{fnCount},3) %y = event
            %        %for z = 1:6
            %         %   params.fpass = freqs{z};
            if data.pfcidx (y)
                %normval = bandpower (data.Field(find(data.goodDetectionMask),y), Fs, freqs{z});
                %                overallPFCRI {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RIrange)),y), Fs, freqs{z});%/normval;
                %               overallPFCRR {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RRrange)),y), Fs, freqs{z});%/normval;
                %                %[S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RIrange)),y), params);%/normval;
                %                %overallPFCRI {x} (y,z) = mean(S);
                %                %[S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RRrange)),y), params);%/normval;
                %                %overallPFCRR {x} (y,z) = mean(S);
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask,y), params);%/normval;
                NormPFC (y,:) = S./data.SimpleNormArray(1,1,y);
                NormPFCFreq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(LLrange),y), params);%/normval;
                data.overallPFCLL  (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallPFCLLFreq  (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(RLrange),y), params);%/normval;
                data.overallPFCRL (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallPFCRLFreq (y,:) = f;
            else
                %normval = bandpower (data.Field(find(data.goodDetectionMask),y), Fs, freqs{z});
                %               overallNACRI {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RIrange)),y), Fs, freqs{z});%/normval;
                %               overallNACRR {x} (y,z) = bandpower (data.Field(find(data.goodDetectionMask(RRrange)),y), Fs, freqs{z});%//normval;
                %                [S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RIrange)),y), params);%/normval;
                %                overallNACRI {x} (y,z) = mean(S);
                %                [S,f]= mtspectrumc (data.Field(find(data.goodDetectionMask(RRrange)),y), params);%/normval;
                %                overallNACRR {x} (y,z) = mean(S);
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(),y), params);%/normval;
                NormNAC (y,:) = S./data.SimpleNormArray(1,1,y);
                NormNACFreq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(LLrange),y), params);%/normval;
                data.overallNACLL (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallNACLLFreq (y,:) = f;
                [S,f]= mtspectrumc (data.Field(data.goodDetectionMask(RLrange),y), params);%/normval;
                data.overallNACRL  (y,:) = S./data.SimpleNormArray(1,1,y);
                data.overallNACRLFreq (y,:) = f;
            end
            %tempChanDataBasicNAC {fnCount} (x,y,z) = bandpower (tempDataRefed{fnCount}(:,nacidx,y), Fs, freqs{z});
            %[S, t, f] = mtspecgramc (tempDataRefed{fnCount}(:,x,y), [1 .1], params);
            %tempChanDataBasic{fnCount} (:,:,y,x) = S;
        end
    end


%% Lever Press Basic Spectrum (and some of the detailed spectrum)
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
        %save ('compiled.mat', 'NormPFC', 'fileList', 'NormPFCFreq', 'overallPFCLL', 'overallPFCLLFreq', 'overallPFCRL', 'overallPFCRLFreq', 'NormNAC', 'NormNACFreq', 'overallNACLL', 'overallNACLLFreq', 'overallNACRL', 'overallNACRLFreq');
    end
end

