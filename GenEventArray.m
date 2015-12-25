function [newEventTS,varargout] = GenEventArray( eventCellArray,var1, varargin  )
%UNTITLED This gets a variety of event timings. Provide the event train,
%the component(s) you are looking for, and the relationship between them
%   GenEventArray is expecting the event array provided by code written
%   with Chris Lapisch. !This should be formalized in the future!. For now,
%   base the existing event array in as a cell array, containing the time
%   text in cell 1 and the timestamps in cell 2. Then, p events to be
%   sorted by, and optionally, the relationship between them. Currently
%   implemented operations include 'opt' -'iso','firstAfter','lastBefore'.
%   The naming of these options could use some serious work. See code for
%   exact details on options


%% Handle the variable number of incoming arguments
opt='none';
optRange=[];
var2='';
ParseInput;
switch opt
    case 'none'
        % Passes all var1 back. Be sure to process as cell or char
        var1TS=getvarTS(var1);
        %         if ~isempty(var2)
        %             var2TS=getvarTS(var2);
        %         end
        %         newEventTS=sort([var1TS var2TS]);
        newEventTS=var1TS;
    case 'iso'
        if isempty(optRange)
            optRange=[-3 3];
        end
        % Gets all of vari one. then makes sure its isolated from var2 by range
        var1TS=getvarTS(var1);
        if ~isempty(var2)
            var2TS=getvarTS(var2);
        else error('var2 is required for ISO to work');
        end
        newEventTS=var1TS(cellfun(@isempty, ...
            arrayfun(@(x) find(var2TS>(x+optRange(1)) & var2TS<(x+optRange(2)) & ...
            abs(var2TS-x)>.0001), var1TS, 'UniformOutput', false)));
    case 'firstAfter'
        % Gets all of vari one, then looks for any of first within range after
        if isempty(optRange)
            optRange=[0 8];
        end
        % Gets all of vari one. then makes sure its isolated from var2 by range
        var1TS=getvarTS(var1);
        if ~isempty(var2)
            var2TS=getvarTS(var2);
        else error('var2 is required for ISO to work');
        end
        
        % Get the first var2 following each var1. then make sure there is
        % only 1 of each var 2
        newEventTS=var2TS(unique(cell2mat(...
            arrayfun(@(x) find(var2TS>x+optRange(1) & var2TS<x+optRange(2),1), ...
            var1TS, 'UniformOutput', false))));
        %%%%%%%%%%%%%%%%%%%%%
        % check for stats out
        %%%%%%%%%%%%%%%%%%%%%
        
        if nargout==2
        % Get the closes('last') var1, to each newEventTS, allowing for lag
        % measurement
        firstEventTS=var1TS(unique(cell2mat(...
            arrayfun(@(x) find(var1TS<x-optRange(1) & var1TS>x-optRange(2),1,'last'), ...
            newEventTS, 'UniformOutput', false))));
        stats.Lag=newEventTS-firstEventTS;
        %stats.IgnoredVar1=length(var1TS)-length(newEventTS);
        varargout{1}=stats;
        end
    case 'lastBefore'
        % what is this even?
        display ('lastbefore has not been implemented yet');
        
    otherwise
        error([varargin{optVarIdx+1} ' is not a valid option input.']);
end

    function [curVarTS] = getvarTS (curVar)
        if iscell(curVar)
            curVarIdx=any(cell2mat(cellfun(@(x) strcmp(x,eventCellArray{2}),curVar,'UniformOutput',false)),2);
            curVarTS=eventCellArray{1}(curVarIdx);
        else
            curVarIdx=any(strcmp(curVar,eventCellArray{2}),2);
            curVarTS=eventCellArray{1}(curVarIdx);
        end
    end
% p=inputParser;
% p.addRequired('eventCellArray',@(x) length(x)>0 && iscell(x));
% p.addRequired('var1',@ischar);
% p.addOptional('var2','',@ischar);
% p.addOptional('optRange',[3 3], @isnumeric);
% p.addParamValue('opt','none',@ischar);
%
% p.parse(eventCellArray, varargin{:});%varargin or varargin{} won't work
% Results=p.Results;
%% old attempt at parser
%nVarargs = length(varargin);
%display(nargin);

% Look for 'opt'
%optVarIdx=findoptVar
% Use switch-case to check for implemented opts
% impOpts={'none' 'iso' 'firstAfter' 'lastBefore'};
% % need to check for string?
%
%
% % check for opt
% isopt = cellfun(@(x) ischar(x) && strcmp(x, 'opt'), varargin);
% if any(isopt)
%     switch varargin{optVarIdx+1}
%         case 'none'
%         case 'iso'
%         case 'firstAfter'
%         case 'lastBefore'
%         otherwise
%             error([varargin{optVarIdx+1} ' is not a valid option input.']);
%             %is it possible to event get here? seems out of index
%     end
% end
%%
    function ParseInput()
        var2Bool=false;
        optBool=false;
        needRangeBool=false;
        optRangeBool=false;
        %%%%%%%%%%%%%%%
        %eventCellArray
        %%%%%%%%%%%%%%%
        
        if ~iscell(eventCellArray)%validate eventArraystruct is a cell array
            error(['EventCellArray is not a cell and is instead a ' class(eventCellArray)]);
        elseif ~xor(size(eventCellArray,1)~=2, size(eventCellArray,2)~=2) ...
                || isempty(eventCellArray{1}) || isempty(eventCellArray{2})%2 components not empty components
            error(['EventCellArray has too few or empty cells ' class(eventCellArray)]);
        elseif ~isnumeric(eventCellArray{1})%the first of which is number
            error(['EventCellArray{1} is not a numeric vector and is instead a ' class(eventCellArray)]);
        elseif ~iscell(eventCellArray{2}) %and the second of which is a cell
            error(['EventCellArray{2} is not a cell and is instead a ' class(eventCellArray)]);
        elseif ~ischar(eventCellArray{2}{1})  %of strings
            error(['EventCellArray is not a cell of strings and is instead a cell of ' class(eventCellArray)]);
        end
        
        %%%%%%%%%%%%%%%%
        % var1
        %%%%%%%%%%%%%%%%
        if isempty(var1)
            error('var1 can not be empty');
        elseif iscell(var1)
            if ~all(cellfun (@ischar,var1))
                error('var1 should be a cell of strings');
            elseif any(cellfun(@isempty,var1))
                error('var1 can not contain any empty cells');
            end
        elseif ischar(var1)
            if strcmp(var1,'opt')||strcmp(var1,'optRange')
                error('var1 can not be opt or optRange');
            end
        else error(['should be a string or cell of strings, but instead is a ' class(var1)]);
        end
        %%%%%%%%
        %Process varargin
        %%%%%%%%%%%%%%%%
        % var2
        %%%%%%%%%%%%%%%%
        if length(varargin)>=1
            if isempty(varargin{1})
                error('var2/opt can not be empty');
            elseif iscell(varargin{1})
                if ~all(cellfun (@ischar,varargin{1}))
                    error('var1 should be a cell of strings');
                elseif any(cellfun(@isempty,varargin{1}))
                    error('var1 can not contain any empty cells');
                else
                    var2=varargin{1};
                    var2Bool=true;
                end
            elseif ischar(varargin{1})
                if ~(strcmp(varargin{1},'opt')||strcmp(varargin{1},'optRange'))
                    var2=varargin{1};
                    var2Bool=true;
                end
            else error('var2/opt can not be empty');
            end
        end
        %%%%%%%%%%%%%%%%
        % 'opt'
        %%%%%%%%%%%%%%%%
        if length(varargin)>=var2Bool+1
            if ~ischar(varargin{1+var2Bool})||isempty(varargin{1+var2Bool})
                error(['opt can not be empty and should be a char array but is a ' class(varargin{1+var2Bool})]);
            elseif strcmp(varargin{1+var2Bool},'opt')
                if length(varargin)>=var2Bool+2
                    if ~ischar(varargin{2+var2Bool})||isempty(varargin{2+var2Bool})
                        error(['opt2 can not be empty and should be a char array but is a ' class(varargin{2+var2Bool})]);
                    else switch varargin{2+var2Bool}
                            case 'none'
                                optBool=true;
                                opt=varargin{2+var2Bool};
                            case 'iso'
                                optBool=true;
                                opt=varargin{2+var2Bool};
                            case 'firstAfter'
                                optBool=true;
                                opt=varargin{2+var2Bool};
                            case 'lastBefore'
                                optBool=true;
                                opt=varargin{2+var2Bool};
                            otherwise
                                error([varargin{2+var2Bool} ' is not a valid option input.']);
                                %is it possible to event get here? seems out of index
                        end
                    end
                else error('opt requires another parameter to function');
                end
            else error([varargin{1+var2Bool} ' was not expected, opt was.']);
            end
        end
    end
end
%
%
%
%
%
% if strcmp(evtTrigger,'RL_C')
%     Press_raw=strmatch('RL_R',behaveEvt_Raw);
%     Press_ts=behaveEvtTm_Raw(Press_raw);
%     BB_raw=strmatch('BB',behaveEvt_Raw);
%     BB_ts=behaveEvtTm_Raw(BB_raw);
%     k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x+8,1), Press_ts, 'UniformOutput', false)));
%     hld=BB_ts(k);
%
% elseif strcmp(evtTrigger,'LL_C')
%     Press_raw=strmatch('LL_R',behaveEvt_Raw);
%     Press_ts=behaveEvtTm_Raw(Press_raw);
%     BB_raw=strmatch('BB',behaveEvt_Raw);
%     BB_ts=behaveEvtTm_Raw(BB_raw);
%     k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x+8,1), Press_ts, 'UniformOutput', false)));
%     hld=BB_ts(k);
% elseif strcmp(evtTrigger,'RL_I')
%     Press_raw=strmatch('RL_U',behaveEvt_Raw);
%     Press_ts=behaveEvtTm_Raw(Press_raw);
%     Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp({'BB' 'RL_R' 'LL_R' 'RL_U' 'LL_U'},x),behaveEvt_Raw,'UniformOutput', false) );
%     Other_ts=behaveEvtTm_Raw(Other_raw);
%     k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
%     hld=Press_ts(k);
%     postEvt=3;
% elseif strcmp(evtTrigger,'LL_I')
%     Press_raw=strmatch('LL_U',behaveEvt_Raw);
%     Press_ts=behaveEvtTm_Raw(Press_raw);
%     Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp({'BB' 'RL_R' 'LL_R' 'RL_U' 'LL_U'},x),behaveEvt_Raw,'UniformOutput', false) );
%     Other_ts=behaveEvtTm_Raw(Other_raw);
%     k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
%     hld=Press_ts(k);
%     postEvt=3;
% else
%     k=strmatch(evtTrigger,behaveEvt_Raw);
%     hld=behaveEvtTm_Raw(k);
% end
%
% end
%
