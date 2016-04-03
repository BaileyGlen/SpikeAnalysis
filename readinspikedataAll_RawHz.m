for animalNumber = 1:4
    for dayNumber = 1:2
        
        load (fileNames13141618D04_D10{1,dayNumber,animalNumber});
  %      filename = 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\25\bintest03_2014-04-19 m25rrri-01_sorted.txt';
        delimiter = ',';
        
        %% Format string for each line of text:
        %   column1: double (%f)
        %	column2: double (%f)
        %   column3: double (%f)
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%f%f%f%[^\n\r]';
        
        %% Open the text file.
        fileID = fopen(fileNames13141618D04_D10{2,dayNumber,animalNumber},'r');
        
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
        goodLL_TS = ones(length(data.LL_TS),1);
        goodRL_TS = ones(length(data.RL_TS),1);
        % Previously Missing Initialization %
        celltest = cell (1);
        % 
        if (animalNumber==3 && dayNumber==1)
            goodLL_TS(1:end-3) = arrayfun (@(x) min(mask(x-(400*0):x+(400*1.5)))  , data.LL_IDX(1:end-3), 'UniformOutput', true);
            goodLL_TS(end-3:end) = 0;
        else
            goodLL_TS = arrayfun (@(x) min(mask(x-(400*0):x+(400*1.5)))  , data.LL_IDX, 'UniformOutput', true);
        end
        %goodLL_TS = arrayfun (@(x) min(mask(x-(400*.5):x+(400*.5)))  , data.LL_IDX, 'UniformOutput', true);
        goodRL_TS = arrayfun (@(x) min(mask(x-(400*0):x+(400*1.5))) , data.RL_IDX, 'UniformOutput', true);
        totcount = 1;
        for numchan = unique(tempdata(:,1)).'
            for numcell=  unique(tempdata(tempdata(:,1)==numchan,2)).'
                shorttemp = tempdata(tempdata(:,1)==numchan & tempdata(:,2)==numcell,:);
                celltest{totcount,1} = shorttemp(:,3).';
                totcount = totcount +1;
            end
        end
        %Previously Missing Initialization
        LLSum1 = cell(1);
        RLSum1 = cell(1);
        
        
        for y = 1:length(celltest)
            
            LLSum1{y,1} = arrayfun (@(x) celltest{y}(celltest{y}> (x) & celltest{y}< (x+1.5))-x,data.LL_TS,'UniformOutput', false);
            RLSum1{y,1} = arrayfun (@(x) celltest{y}(celltest{y}> (x) & celltest{y}< (x+1.5))-x,data.RL_TS,'UniformOutput', false);
            %    LLSum2{y,1} = arrayfun (@(x) celltest{y}(celltest{y}> (x+2) & celltest{y}< (x+4))-x+2,data.LL_TS,'UniformOutput', false);
            %    RLSum2{y,1} = arrayfun (@(x) celltest{y}(celltest{y}> (x+2) & celltest{y}< (x+4))-x+2,data.RL_TS,'UniformOutput', false);
            %val = cellfun(@(x)  x(x>data.LL_TS((y))-3 & x<data.LL_TS((y))+3)-data.LL_TS((y)),celltest,'UniformOutput', false);
        end
        LLSumArray{y,1}= [];
        RLSumArray{y,1} = [];
        for y = 1:length(celltest)
            LLSumArray1{y,1} = arrayfun (@(x) cat(2, LLSumArray{y,1}, celltest{y}(celltest{y}> (x) & celltest{y}< (x+1.5))-x),data.LL_TS,'UniformOutput', false);
            %    LLSumArray2{y,1} = arrayfun (@(x) cat(2, LLSumArray{y,1}, celltest{y}(celltest{y}> (x) & celltest{y}< (x+2))-x),data.LL_TS,'UniformOutput', false);
            
            RLSumArray1{y,1} = arrayfun (@(x) cat(2, RLSumArray{y,1}, celltest{y}(celltest{y}> (x) & celltest{y}< (x+1.5))-x),data.RL_TS ,'UniformOutput', false);
            %    RLSumArray2{y,1} = arrayfun (@(x) cat(2, RLSumArray{y,1}, celltest{y}(celltest{y}> (x+2) & celltest{y}< (x+4))-x+2),data.RL_TS ,'UniformOutput', false);
            
            %val = cellfun(@(x)  x(x>data.LL_TS((y))-3 & x<data.LL_TS((y))+3)-data.LL_TS((y)),celltest,'UniformOutput', false);
        end
        clearvars numtest totcount numchan numcell;
        A_LLRL_changecell_RfDel_Raw_Baseline {animalNumber,dayNumber} = cellfun(@(y) length(y)/((length(mask)/400)),celltest);
        A_LL_changecell_RfDel_Raw{animalNumber,dayNumber}(:,1) = cellfun (@(x) length(cat(2,x{goodLL_TS & data.LL_Rf_Mask}))/length(find(goodLL_TS & data.LL_Rf_Mask)), LLSum1)/1.5;
        A_LL_changecell_RfDel_Raw{animalNumber,dayNumber}(:,2) = cellfun (@(x) length(cat(2,x{goodLL_TS & ~data.LL_Rf_Mask & data.burstInfo.LL_Last}))/length(find(goodLL_TS & ~data.LL_Rf_Mask & data.burstInfo.LL_Last)), LLSum1)/1.5;
        
        A_RL_changecell_RfDel_Raw{animalNumber,dayNumber}(:,1) = cellfun (@(x) length(cat(2,x{goodRL_TS & data.RL_Rf_Mask}))/length(find(goodRL_TS & data.RL_Rf_Mask)), RLSum1)/1.5;
        A_RL_changecell_RfDel_Raw{animalNumber,dayNumber}(:,2) = cellfun (@(x) length(cat(2,x{goodRL_TS & ~data.RL_Rf_Mask & data.burstInfo.RL_Last}))/length(find(goodRL_TS & ~data.RL_Rf_Mask & data.burstInfo.RL_Last)), RLSum1)/1.5;
        %
        % %%
    end
end
