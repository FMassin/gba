classdef traceList < handle
    %TRACELIST A list of wform meta-info; everything but the wform itself. 
    %   Version: i32/
    
    properties
        
        fullName        % Waveform filename including absolute path
        
        eqLat           % [deg]
        eqLon           % [deg]
        eqZ             % [km]
        eqDate     
        t0              % Origin time
        ts              % First value of time vector relative to t0 (t(t0)=0)
        
        m           
        mType       
        eqType      
        eqMech     
        eqIdx           % Event index

        epiDist         % [km]
        hypDist         % [km]
        
        stationName
        networkName
        stationLat      % [deg]
        stationLon      % [deg]
        stationAlt      % [m]
        stationIdx      % [#]
        
        sRate      
        %bandCode		% Sensor codes [H]
        instrCode		% [H,L]
        orntCode		% [E,N,Z]
        
        dataSetName
        
        tppx            % p-wave pick in [sec] after origin time
        tspx            % s-wave pick in [sec] after origin time
        ppxIdx          % Index of p-wave pick in time vector
        spxIdx          % Index of s-wave pick in time vector
        hasPpx          % Flag if trace has already been picked
        hasSpx          % Flag if trace has already been picked
        allPx           % Structure for saving other picks
        
        amax            % Maximum amplitudes in each frequency band
        cav             % Cumulative absolute velocities in each frequency band
        pga             % Peak acclereation value after p-pick in a broad [0.1 49]Hz frequency band
        pgv             % Peak velocity     value after p-pick in a broad [0.1 49]Hz frequency band
        pgd             % Peak displacement value after p-pick in a broad [0.1 49]Hz frequency band
        pgaIdx          % vector index of PGA
        pgvIdx          % vector index of PGV
        pgdIdx          % vector index of PGD
        pgv_up          % Upper bound (2sigma) of pgv from CH2007
        pgv_lo          % Lower bound (2sigma) of pgv from CH2007
        nbpga           % narrow band PGA (= max amplitude over all of record in each frequency band)
        nbpgv           % narrow band PGV (= max amplitude over all of record in each frequency band)
        nbpgd           % narrow band PGD (= max amplitude over all of record in each frequency band)
        
        snr             % Signal to noise ratio (on velocity trace)
        nbnoise         % Subband noise level before p-pick (--> in all <nband> frequency bands, on velocity trace)
        dspNoise       % Noise on vertical displacement trace
        velNoise       % Noise on vertical velocity trace
        accNoise       % Noise on vertical acceleration trace
       
        pd              % peak displacement          
        pv              % peak velocity
        pa              % peak acceleration

        mErr            % Magnitude prediction errors at <r> different times with <s> different methods <r-by-s>
        rErr            % Distance        "                 "            " 
        
        tauP
        tauC

        flt             % Filter information

        trec            % absolute time of record start
        relIdx          % relative line-index for co-recorded traces
        
        comment         % e.g. skipReasons for skipped traces
        var1            % Placeholder for unexpected variables 
        var2
        var3
        var4
        var5
    end
    
    methods
        
        % Invoke object   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
        function obj = traceList(ntr)
            
            obj.fullName    = cell(ntr,1);
            
            obj.eqLat       = zeros(ntr,1,'single');
            obj.eqLon       = zeros(ntr,1,'single');
            obj.eqZ         = zeros(ntr,1,'single');
            obj.eqDate      = cell(ntr,1);
            obj.t0          = cell(ntr,1);
            obj.ts          = zeros(ntr,1,'single');
            
            obj.m           = zeros(ntr,1,'single');
            obj.mType       = cell(ntr,1);
            obj.eqType      = cell(ntr,1);
            obj.eqMech      = cell(ntr,1);
            obj.eqIdx       = zeros(ntr,1,'single');

            obj.epiDist     = zeros(ntr,1,'single');
            obj.hypDist     = zeros(ntr,1,'single');
            
            obj.stationName = cell(ntr,1);
            obj.networkName = cell(ntr,1);
            obj.stationLat  = zeros(ntr,1,'single');
            obj.stationLon  = zeros(ntr,1,'single');
            obj.stationAlt  = zeros(ntr,1,'single');
            obj.stationIdx  = zeros(ntr,1,'single');
            
            obj.sRate      = zeros(ntr,1,'single');
            %obj.bandCode   = cell(ntr,1);
            obj.instrCode  = cell(ntr,1);
            obj.orntCode   = cell(ntr,1);
            
            obj.dataSetName = cell(ntr,1);
            
            obj.tppx        = zeros(ntr,1,'single');
            obj.tspx        = zeros(ntr,1,'single');
            obj.ppxIdx      = zeros(ntr,1,'single');
            obj.spxIdx      = zeros(ntr,1,'single');
            obj.hasPpx      = false(ntr,1);
            obj.hasSpx      = false(ntr,1);
            obj.allPx       = cell(ntr,1);
            
            obj.amax        = cell(ntr,1);
            obj.cav         = cell(ntr,1);
            obj.pga         = zeros(ntr,1,'single');
            obj.pgv         = zeros(ntr,1,'single');
            obj.pgd         = zeros(ntr,1,'single');
            obj.pgaIdx      = zeros(ntr,1,'single');
            obj.pgvIdx      = zeros(ntr,1,'single');
            obj.pgdIdx      = zeros(ntr,1,'single');
            
            obj.pgv_up      = zeros(ntr,1,'single');
            obj.pgv_lo      = zeros(ntr,1,'single');
            obj.nbpga       = cell(ntr,1);
            obj.nbpgv       = cell(ntr,1);
            obj.nbpgd       = cell(ntr,1);
            
            obj.snr         = zeros(ntr,1,'single');
            obj.nbnoise     = cell(ntr,1);
            obj.dspNoise    = cell(ntr,1);
            obj.velNoise    = cell(ntr,1);
            obj.accNoise    = cell(ntr,1);

            obj.pd          = cell(ntr,1);
            obj.pv          = cell(ntr,1);
            obj.pa          = cell(ntr,1);

            obj.mErr        = cell(ntr,1);
            obj.rErr        = cell(ntr,1);
            
            obj.tauC        = cell(ntr,1);
            obj.tauP        = cell(ntr,1);

            obj.flt         = cell(ntr,1);
            
            obj.trec        = cell(ntr,1);          
            obj.relIdx      = cell(ntr,1);          
            
            obj.comment     = cell(ntr,1);
            obj.var1        = cell(ntr,1);
            obj.var2        = cell(ntr,1);
            obj.var3        = cell(ntr,1);
            obj.var4        = cell(ntr,1);
            obj.var5        = cell(ntr,1);
        end
        
        % Clear out empty entries  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        function consolidate(obj)
            idx_lastEntry = find(cellfun(@(x) ~isempty(x), obj.fullName),1,'last');
            names         = fieldnames(obj);
            nf            = size(names,1);
            for i = 1:nf
                obj.(names{i}) = obj.(names{i})(1:idx_lastEntry);
            end
        end
        
        % Remove skipped traces from Global list  -  -  -  -  -  -  -  -  -  
        function removeSkipped(obj,flg_skipTrace)
            
            nflgVals = numel(flg_skipTrace);
            ntraces  = numel(obj.m);     
            if nflgVals~=ntraces
                error('Use logicals/flags for specifing traces to skip, not indices\n')
            end
            
            dontSkip = (flg_skipTrace==0);
            
            names = fieldnames(obj);
            nf    = size(names,1);
            ntr   = numel(obj.(names{1}));
            
            for i = 1:nf
                nentries = numel(obj.(names{i}));
                if nentries==ntr
                    obj.(names{i}) = obj.(names{i})(dontSkip);
                else 
                    fprintf(1,sprintf('Note: Field %s only has %i entries; not removing any entries in this field\n',names{i},nentries))
                end
            end
        end

        % Append lists to existing one   -  -  -  -  -  -  -  -  -  -  -  -  
        function appendList(obj,NewList)
            names = fieldnames(obj);
            nf    = size(names,1);
            for i = 1:nf
                obj.(names{i}) = [obj.(names{i}); NewList.(names{i})];
            end
        end

%         % Delete fields except the ones in <keepFieldNames>     -  -  -  -  -  -  
%         function keepOnly(obj,keepFieldNames)
%             names = fieldnames(obj);
%             nf    = size(names,1);
%             for i = 1:nf
%                 fieldName = names{i};
%                 if ~ismember(fieldName,keepFieldNames)
%                     obj = rmfield(obj,fieldName);
%                     %obj.(names{i}) = [obj.(names{i}); NewList.(names{i})];
%                 end
%             end
%         end
        
        % Select certain entries of GlobalList and output them as a new list
        % note: can also be used to obtain independent clone list to which
        % changes can be made without changing the original list
        function [subList] = selectSubList(obj,selectIndex) 
            if (islogical(selectIndex))
                error('Use indices instead of logical vectors to select entries from traceList.')
            end
            names   = fieldnames(obj);
            nf      = size(names,1);
            ntr     = numel(obj.(names{1}));
            ntr_sub = numel(selectIndex);
            subList = traceList(ntr_sub);      %generate new object
            for i = 1:nf
                
                nentries = numel(obj.(names{i})(:,1));
                %nentries = numel(obj.(names{i}));
                if nentries==ntr
                    subList.(names{i}) = obj.(names{i})(selectIndex,:);
                    %subList.(names{i}) = obj.(names{i})(selectIndex);
                else
                    subList.(names{i}) = 0;
                    %fprintf(1,sprintf('Note: Field %s only has %i entries; set this field to zero in new subList\n',names{i},nentries))
                end
            end
        end
        
        % Change GlobalList to only contain certain entries   -  -  -  -  -  
        function shrinkList(obj,selectFlag) 
            names = fieldnames(obj);
            nf    = size(names,1);
            ntr   = numel(obj.(names{1}));
            for i = 1:nf
                nentries = numel(obj.(names{i}));
                if nentries==ntr
                    obj.(names{i}) = obj.(names{i})(selectFlag);
                else
                    fprintf(1,sprintf('Note: Field %s only has %i entries; not removing any entries in this field\n',names{i},nentries))
                end
            end
        end
        
        % Sort lists  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
        function sortList(obj,sortIndex)
            names = fieldnames(obj);
            nf    = size(names,1);
            for i = 1:nf
                obj.(names{i}) = obj.(names{i})(sortIndex);
            end
        end
        
        % Add single or several entries     -  -  -  -  -  -  -  -  -  -  -  
        function addEntries(obj,sourceList,traceIdx)
            idx_nextline = find(cellfun(@(x) isempty(x), obj.fullName),1,'first');
            nadds = numel(traceIdx);
            names = fieldnames(obj);
            nf    = size(names,1);
            for i = 1:nf
                obj.(names{i})(idx_nextline:idx_nextline+nadds-1) = sourceList.(names{i})(traceIdx);
            end
        end
        
        % Add minimum infor of single entry which is to be skipped  -  -  -
        function addSkipEntry(obj,fullName,reason,hypDist,m)
            idx_nextline = find(cellfun(@(x) isempty(x), obj.fullName),1,'first');
            obj.fullName{idx_nextline}   = fullName;
            obj.comment{idx_nextline}    = reason;
            obj.hypDist(idx_nextline)    = hypDist;
            obj.m(idx_nextline)          = m;
        end
        
        % Return independent clone of the list
        % use obj.selectSublist for this purpose
        
        % Return a clone of the list, but containing only the "fullName" property
        function [nameList] = cloneNameList(obj)
            nameList          = traceNameList(size(obj.fullName,1));	%generate new object
            nameList.fullName = obj.fullName;
            nameList.ppxIdx   = obj.ppxIdx;
            nameList.tppx     = obj.tppx;
        end

        % Print most important info 
        function printSingleTraceSummary(obj,idx)
            fprintf(1,['\n\nFile name:\t\t\t',obj.fullName{idx},'\n'])
            fprintf(1,['EQ-coords (lat/lon/z):\t\t',num2str([obj.eqLat(idx),obj.eqLon(idx),obj.eqZ(idx)]),'km\n'])
            fprintf(1,['Magnitude:\t\t\t',num2str(obj.m(idx)),'\t(of type ',obj.mType{idx},')\n'])
            fprintf(1,['Hypocentral Distance:\t\t',num2str(obj.hypDist(idx),'%4.1f'),'km\n'])
            fprintf(1,['Epicentral  Distance:\t\t',num2str(obj.epiDist(idx),'%4.1f'),'km\n'])
            fprintf(1,['Station-coords (lat/lon/alt):\t',num2str([obj.stationLat(idx),obj.stationLon(idx),obj.stationAlt(idx)]),'m\n'])
            fprintf(1,['Station-specifications:\t\torientation code ',obj.orntCode{idx},'\tinstrument code ',obj.instrCode{idx},'\n'])
            fprintf(1,['Comments:\t\t\t',obj.comment{idx},'\n'])
            fprintf(1,['-------------------------------------------------------------------------\n'])
        end

        % Update relative indices  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        function updateRelIdx(obj)
            fprintf(1,'Not yet written. do it.\n')
        end
        
        
        % Compute and print size of traceList  -  -  -  -  -  -  -  -  -  -
        function [totSize] = printObjectSize(obj)
            props   = properties(obj);
            totSize = 0;
            
            for ii=1:length(props)
                currentProperty = getfield(obj, char(props(ii)));
                s               = whos('currentProperty');
                totSize         = totSize + s.bytes;
            end
            
            if (totSize>1e9)
                fprintf(1, '%3.1f GB\n', totSize*1e-9)
            else
                fprintf(1, '%3.1f MB\n', totSize*1e-6)
            end
        end
        
        
        
        
    end
end

