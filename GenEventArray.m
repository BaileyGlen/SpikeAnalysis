function [] = GenEventArray( eventCellArray,var1, varargin )
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
ParseInput;
% p=inputParser;
% p.addRequired('eventCellArray',@(x) length(x)>0 && iscell(x));
% p.addRequired('var1',@ischar);
% p.addOptional('var2','',@ischar);
% p.addOptional('range',[3 3], @isnumeric);
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
        var2bool=false;
        optbool=false;
        needRangebool=false;
        rangebool=false;
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
        if ~ischar(var1)||isempty(var1)
            if iscell(
            error(['var1 can not be empty and should be a char array ' class(eventCellArray)]);
        end
        
        %%%%%%%%%%%%%%%%
        % var2
        %%%%%%%%%%%%%%%%
        if ~ischar(varargin{1})||isempty(varargin{1})
            error(['var2/opt can not be empty and should be a char array ' class(eventCellArray)]);
        elseif ~strcmp(varargin{1},'opt')&&~strcmp(varargin{1},'range')
            var2=varargin{1};
            var2Bool=true;
        end
        %%%%%%%%%%%%%%%%
        % 'opt'
        %%%%%%%%%%%%%%%%
        if ~ischar(varargin{1+var2bool})||isempty(varargin{1+var2bool})
            error(['opt can not be empty and should be a char array but is a ' class(eventCellArray)]);
        elseif ~strcmp(varargin{1+var2bool},'opt')
            error([varargin{1+var2bool} ' was not expected, opt was.']);
        elseif ~ischar(varargin{2+var2bool})||isempty(varargin{2+var2bool})
            error(['opt2 can not be empty and should be a char array but is a ' class(eventCellArray)]);
        else switch varargin{2+var2bool}
                case 'none'
                    optbool=true;
                    opt=varargin{2+var2bool};
                case 'iso'
                    optbool=true;
                    opt=varargin{2+var2bool};
                case 'firstAfter'
                    optbool=true;
                    opt=varargin{2+var2bool};
                case 'lastBefore'
                    optbool=true;
                    opt=varargin{2+var2bool};
                otherwise
                    error([varargin{2+var2bool} ' is not a valid option input.']);
                    %is it possible to event get here? seems out of index
            end
        end    var2=varargin{1};
        
        var2Bool=true;
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
