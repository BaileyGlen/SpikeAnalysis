%% NEWEST VERSION OF SELECTIVITY_ D.N.Linsenbardt (November 2015) dlinsen1@gmail.com
close all;
for t=1:4; 
    clearvars ('-except', 'data', 't');   
    %cd 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\CL';
    if t==1;
        %load('LLD04_pEvt02.mat');
        crtVar = 'LL04';
    elseif t==2;
        %load('LLD10_pEvt02.mat');
        crtVar = 'RL04';
    elseif t==3;
        %load('RLD04_pEvt02.mat');
        crtVar = 'LL10';
    else
        %load('RLD10_pEvt02.mat');
        crtVar = 'RL10';
    end;
        
        % GETTING INTO SAME FORMAT AS MY PEVTMTX
        pEvt=data.(crtVar).pEvt';
        
        for i=1:size(pEvt,2);
            NEWpEvt{i}=pEvt{i}';
        end;
        
        pEvtMtx=NEWpEvt;
        
        % CREATING A BASELINE MATRIX TO ADD TO THE BEGGINING OF PEVMTX
%         pEvt_base=data.(crtVar).pEvt_base';
%         
%         for i=1:size(pEvt_base,2);
%             NEWpEvt_base{i}=pEvt_base{i}';
%         end;
%         
%         pEvtMtx_base=NEWpEvt_base;
%                 
         lLim=30;
                
        %% NEW BINS LISTED BELOW
        
        % if t==1;
        epoch=1:30+1; %(3 SECONDS BEFORE AND AFTER LEVER PRESS)
        %     epoch=1:60+1; %(13SEC BASE VS 3 SEC LEVER PRESSING ETC.)
        % elseif t==2;
        %     epoch=20+1:40+1; %
        % elseif t==3;
        %     epoch=40+1:50+1; %
        % elseif t==4;
        %     epoch=50+1:65+1; %
        % elseif t==5;
        %     epoch=65+1:155+1; %
        % elseif t==6;
        %     epoch=155+1:170+1; %
        % elseif t==7;
        %     epoch=170+1:190+1; %
        % elseif t==8;
        %     epoch=20+1:170+1; %
        % else
        %     epoch=20+1:190+1; %
        % end;
        
        m=1;
        h = waitbar(0,'Initializing waitbar...');  % 
        for i=1:size(pEvtMtx,2); % 
            %     i
            PercentComplete = i/6;  waitbar(PercentComplete,h,sprintf('%f%% Complete...D PRIME + ROC',PercentComplete));  % 
            for j=1:size(pEvtMtx{i},2) % 
                
                CellId{m}=[int2str(i) '_' int2str(j)];m=m+1;
                Cellmean=nanmean(pEvtMtx{i}{j},1); % 
                CellmeanPre=nanmean(Cellmean(:,1:lLim),2); % 
                CellmeanPost=nanmean(Cellmean(:,lLim+1:lLim*2),2); % 
                
                %Cellvariance=nanmean(pEvtMtx{i}{j},1); % 
                %Cellvariance_base=nanmean(pEvtMtx_base{i}{j}); % the variance is the squared stdev, so no need to square it later on.
                CellvariancePre=nanvar(Cellmean(:,1:lLim),[],2); % 
                CellvariancePost=nanvar(Cellmean(:,lLim+1:lLim*2),[],2); % 
                
                %sdPre=nanstd(CellvariancePre); % 
                %sdCs=nanstd(CellvariancePost); % 
                
                CellMeanDiff=abs(CellmeanPre-CellmeanPost); %
                CellSqrvariancePre=CellvariancePre; % the stddev^2 = the variance, no need to do anything
                CellSqrvariancePost=CellvariancePost; % 
                CellSumSqrvariance=sqrt((CellSqrvariancePre+CellSqrvariancePost)/2); %
                dSlctIdx{i}(1,j)=(CellMeanDiff)/(CellSumSqrvariance); % computes d: divides the absolute values of the mean differences by the sum of the squared deviations. 
                %% COMPUTING ROC
%                 mnPre=mean(Cellvariance_base(1:lLim)); A_mnPre{i}(j)=mnPre; % 
%                 mnCs=mean(Cellvariance(epoch));A_mnCs{i}(j)=mnCs; % 
%                 sdPre=std(Cellvariance_base(1:lLim)); A_sdPre{i}(j)=sdPre; % 
%                 sdCs=std(Cellvariance(epoch)); A_sdCs{i}(j)=sdCs;  %
%                 
%                 %% COMPUTING p(HITS) + p(MISS)
%                 zVal=[-10:0.1:10];
%                 for z=1:length(zVal);
%                     x(z)= mnCs + (zVal(z)*sdCs);
%                     zH1(z) = normcdf(x(z),mnCs,sdCs); % HITS
%                     zHo(z) = normcdf(x(z),mnPre,sdPre); % MISS
%                 end;
%                                 
%                 %% MAKE ROC AND INTEGRATE AREA - TRAPZ MIGHT BE MAKING THE ROC FUNKY ON SOME NEURONS...
%                 iRoc{i}(1,j)=trapz(zH1,zHo);
%                 rocVals{i}{j}=[zH1;zHo];
                
            end;
        end;
        close(h); 
        
        %% Create surrogate
        h = waitbar(0,'Initializing waitbar...');  % start waitbar to show progress
        for y=1:500; % this changes the number of times surrogates are run (use 1:1 if you dont need surrogate data).
            %     y
            PercentComplete = y/500;  waitbar(PercentComplete,h,sprintf('%f%% Complete...SURROGATE',PercentComplete));  % 
            for i=1:size(pEvtMtx,2); % 
                for j=1:size(pEvtMtx{i},2) % 
                    
                    %% Scrambling data sets
                    CellmeanOrig=nanmean(pEvtMtx{i}{j},1); 
                    x=CellmeanOrig(:);
                    k=length(x);
                    x=x(randperm(k));
                    surMtx=reshape(x,size(CellmeanOrig));
                    
                    
                    surMtx_pre=surMtx(:,1:lLim);
                    surMtx_post=surMtx(:,lLim+1:lLim*2);
                  
                    
                    CellmeanPre=nanmean(surMtx_pre,2);
                    CellmeanPost=nanmean(surMtx_post,2);
%                     Cellmean=nanmean(surMtx_post); % 
%                     CellmeanPre=nanmean(surMtx_pre(1:lLim)); % 
%                     CellmeanPost=nanmean(Cellmean(epoch)); % 
                    
                    %Cellvariance=nanmean(surMtx_post); % 
                    %Cellvariance_pre=nanmean(surMtx_pre);
                    CellvariancePre=nanvar(surMtx_pre,[],2); % 
                    CellvariancePost=nanvar(surMtx_post,[],2); % 
                    
                    CellMeanDiff=abs(CellmeanPre-CellmeanPost); % 
                    CellSqrvariancePre=CellvariancePre; % 
                    CellSqrvariancePost=CellvariancePost; %
                    CellSumSqrvariance=sqrt((CellSqrvariancePre+CellSqrvariancePost)/2); %
                    surSlctIdx{y,i}(1,j)=(CellMeanDiff)/(CellSumSqrvariance); % 
                                                                                
                end;
            end;
        end;
        close(h);  
        
        %% plotting actual versus surrogate data - NEED TO REMOVE NANS AND INFS HERE; THEN NEED TO USE SURMTX TO GET STD AND MEAN TO GET CI
        
        surMtx=cell2mat(surSlctIdx);        
        surMtx_mean=mean(surMtx);        
        surMtx_std=std(surMtx);
        
        % GETTING THE UPPER CI FOR EACH NEURON AND THEN THE MEAN OF THE UPPER CI
        %CI=((surMtx_mean(:))+(surMtx_std(:)*2.575));% 
        CI=((surMtx_mean(:))+(surMtx_std(:)*1.96));% 
        CI_MEAN=nanmean(CI); % THIS WILL BE USED AS THRESHOLD FOR SELECTIVITY - CONSERVATIVE THRESHOLD
        
        % SORTING NEURONS FROM LOWEST DPRIME TO HIGHEST
        [sortData,Idx]=sort(cell2mat(dSlctIdx)); % 
        
        %% NEED TO TAKE THE EACH NEURON THRESHOLD NOW INSTEAD OF ALL NEURON THRESHOLD AND COUNT THESE AS SIGNIFICNATLY SELECTIVE...
        CI_LOWER=((surMtx_mean(:))-(surMtx_std(:)*1.96));
        DP=cell2mat(dSlctIdx);
        DP=DP';
        CI_vs_DP=DP-CI;
                
        [sortCIvDP,CIDPIdx]=sort(CI_vs_DP);
        CellId=CellId(CIDPIdx); % SORT CELL ID BY SIGNIFICANCE

        % TRYING SOME PLOTS...
        %hold off;
        h=figure;
        plot(sortCIvDP,'ko'); hold on;
        plot(DP(CIDPIdx,:),'ro')
        x=zeros([1 60]); plot(x,'k:','LineWidth',2)
        plot(CI_LOWER(CIDPIdx,:),'b'); plot(CI(CIDPIdx,:),'b')
        xlim([0 51]);
        plot(repmat(CI_MEAN,1,size(surMtx,2)),'g:','LineWidth',2);
        title 'Neuron Selectivity [sorted by significance]'
        
        % SIGNIFICANCE WITH INDIVIDUAL CI
        SIG=find(sortCIvDP>0);
        size(SIG)
        NOTSIG=find(sortCIvDP<=0);
        size(NOTSIG)
        
        (length(SIG))/(length(NOTSIG)+length(SIG))*(100)
        
        % SIGNIFICANCE WITH MEAN CI (OLD WAY)
        OLDSIG=find(DP>CI_MEAN);
        size(OLDSIG)
        OLDNOTSIG=find(DP<=CI_MEAN);
        size(OLDNOTSIG)
        
        (length(OLDSIG))/(length(OLDNOTSIG)+length(OLDSIG))*(100)
        
%         rocMtx=cell2mat(iRoc);
%         rocMtx=rocMtx(Idx);
%         figure;
%         scatter(sortData,rocMtx);
        
        cd 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\CL';
        save(sprintf('MkSel_Res_BAILEY__%02d',t));
%        savefig(sprintf('FIG_NODRINK%d.fig',t),h); % will create FIG1, FIG2,...
        data.(crtVar).DPrime.First3.CellId = CellId;
        data.(crtVar).DPrime.First3.SIG = SIG;        
        data.(crtVar).DPrime.First3.CIDPIdx = CIDPIdx;
        %clear all;
        %close all;
end;
