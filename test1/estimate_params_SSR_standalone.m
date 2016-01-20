function [mle] = estimate_params_SSR_standalone(zTarget,hTarget,train,wti,snippet,bands,comps,nsim,r0,logPrior)

% Find nsim most similar traces, according to the (weighted?) L2-norm
% difference. Approximate the m- and r-values of those traces with a
% bivariate Gaussian, optionally combine it with additional likelihood
% functions L(m) and/or L(r) of the two parameters, extract profile
% likelihoods Lp(m) & Lp(r), maximise them to obtain most likely parameter
% estimates mHat & rHat.
%
% menandrin@gmail.com, 140826

% bands         list of bands used for computing misfit between target and
%               training traces
% components    Z and/or H, list of components used for computing misfit 
%               between target and training traces
%
% r0.mean       mean   of independent r-constraint
% r0.sigma      1sigma of independent r-constraint
%
%
% Outstanding task: EXCLUDE BANDS 1 AND 9 FOR NGA TRACES (not yet done)

                        
global o minimax mm rr snpLength fc fMode fOrder iN

o.ppt  = 1;  % If plot is for ppt (landscape format), use o.ppt = 1;
%o.ppt = 0;  % If plot is for poster (portrait format)

fprintf(1,'TMP: hard-coding processing options. check. now. do it.\n')
optRead.process = 1;
optRead.intMode = 'afterPx';

figNum = 999;

ftSize = 14;

mmin    = min(mm);
mmax    = max(mm);
[RR,MM] = meshgrid(rr,mm); 

m_ctlg = zTarget.m;
r_ctlg = zTarget.hypDist;
%r_ctlg = zTarget.epiDist;

if any(strcmp('Z',comps)); useZ=true; else useZ=false; end
if any(strcmp('H',comps)); useH=true; else useH=false; end

% 2.b Find most similar traces
%   Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z  -  Z
if useZ
    az_obs      = zTarget.amax{1}(bands,snippet);
    az_training = train.az(:,bands); 
    %dlmwrite('az_obs.txt',az_obs)
    %dlmwrite('az_training.txt',az_training)
    nt          = size(az_training,1);
%    misfit      = az_training - repmat(az_obs',nt,1);
    misfit      = log10(az_training) - log10(repmat(az_obs',nt,1));
    misfit_l2   = sum(misfit.^2,2); % Cumulative misfit over all frequency bands
    
    % Extract source params (m,r) of nsim most similar traces
    [~,idx_std] = sort(misfit_l2);
    m           =       train.m(idx_std(1:nsim));
    r           = log10(train.r(idx_std(1:nsim)));
    wt          =           wti(idx_std(1:nsim));
    
    % gcf; clf; plot(misfit_l2(idx_std)); set(gca,'yscale','log');
    
    % Weighting: Replicate each value according to it's integer weight in <wt>
    mz = rude(wt,m)';
    rz = rude(wt,r)';
    if o.plotSSR; mzbak=m; rzbak=r; end
end

%   H  -  H  -  H  -  H  -  H  -  H  -  H  -  H  -  H  -  H  -  H  -  H
if useH
    ah_obs      = hTarget.amax{1}(bands,snippet);
    ah_training = train.ah(:,bands); 
    nt          = size(ah_training,1);
    %misfit      = ah_training - repmat(ah_obs',nt,1);
    misfit      = log10(ah_training) - log10(repmat(ah_obs',nt,1));
    misfit_l2   = sum(misfit.^2,2); % Cumulative misfit over all frequency bands
    
    % L2 norm cumulative misfit
    [~,idx_std] = sort(misfit_l2);
    m           =       train.m(idx_std(1:nsim));
    r           = log10(train.r(idx_std(1:nsim)));
    wt          =           wti(idx_std(1:nsim));
    
    % gcf; clf; plot(misfit_l2(idx_std)); set(gca,'yscale','log');
    
    % Weighting: Replicate each value according to it's integer weight in <wt>
    mh = rude(wt,m)';
    rh = rude(wt,r)';
    
    % Back up unweighted values for plotting
    if o.plotSSR; mhbak=m; rhbak=r; end
end

if (useZ && useH)               % ... use both verticals and horizontals
    mVal = [mz; mh];
    rVal = [rz; rh];
    if o.plotSSR
        aVal = [az_training; ah_training]; 
        mbak = [mzbak; mhbak];
        rbak = [rzbak; rhbak];
    end

elseif (useZ && ~useH)          % ... use only verticals
    mVal = mz;
    rVal = rz;
    if o.plotSSR
        aVal = az_training;
        mbak = mzbak;
        rbak = rzbak;
    end
    
elseif (~useZ && useH)          % ... use only horizontals
    mVal = mh;
    rVal = rh;
    if o.plotSSR
        aVal = ah_training;
        mbak = mhbak;
        rbak = rhbak;
    end
end
%   -     -     -     -     -     -     -     -     -     -     -     -


% Approximate m- and r-values with bivariate Gaussian
covmat = cov([rVal mVal]);
mmean  = mean(mVal);
rmean  = mean(rVal);
l      = log(mvnpdf([RR(:) MM(:)],[rmean mmean],covmat));
l      = reshape(l,length(mm),length(rr));   % log-likelihood function on grid

if o.usePrior
    l = l + logPrior;
    figNum = 998;
end

% Compute mle's before adding additional R-constraint
if o.plotSSR; 
    [mHat_noR,rHat_noR,~,~] = max_llh_function(l);
end

% Add hypothetical R-constraint, i.e. combine l(r,m) with l(r)
if o.useRconstraint
    lr = log(normpdf(10.^RR,r0.mean,r0.sigma));
    l  = l+lr;
end

% Relative loglikelihoods
%lr  = l-max(l(:));

% Profile log-likelihood for magnitude and distance lp(m) & lp(r)
% [mHat,rHat,lpmn,lprn] = max_llh_function(l);
% mle.mHat = mHat;
% mle.rHat = rHat;
% mle.lpmn = lpmn;
% mle.lprn = lprn;


% Marginal log-likelihood for magnitude and distance lp(m) & lp(r)
[mHat,rHat,lmmn,lmrn] = max_llh_function_marginal(l);
mle.mHat = mHat;
mle.rHat = rHat;
mle.lmmn = lmmn;
mle.lmrn = lmrn;











%% Plot
if o.plotSSR
    
    
    % Intitiate figure ....................................................
    figure(figNum); clf; hold on; 
    whitebg('w'); set(gcf,'InvertHardcopy','off')

    px = load_picker_settings('production_short');

    shiftFraction =.1;
    shrinkFactor  = .9;
    %level         = -2.3;           % Make colormap saturate at <level>
    level         = -4.6;
    sWindow       = 1;              % Window lenght for measuring signal amps after p-pick [sec]
    nWindow       = 1;              % Window lenght for measuring noise before p-pick [sec]
    npad          = 100;
    tFrac         = 0.1;
    T             = 12;
    
    
    
    % ----------------------------------------
    % PLOT FIGURE FOR POSTER (PORTRAIT FORMAT)
    % ----------------------------------------
    if ~o.ppt
        
        fprintf(1,'Plotting script was updated for i36 for ppt-mode only. Modify portrait mode if needed.\n')
        
        % Waveform plots
        axWf = zeros(10,1);
        for i=1:10
            ii = (i-1)*7+1;
            axWf(i) = subplot(10,7,ii:ii+2);
        end
        
        % Spectra plot
        axSpect = subplot(10,7,[4:6,11:13,18:20]);
        pos     = get(axSpect,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axSpect,'position',pos)
        
        % Data plot
        axDat   = subplot(10,7,[25:27,32:34,39:41]);
        pos     = get(axDat,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axDat,'position',pos)
        
        % Magnitude profile likelihood function Lpm
        axLpm   = subplot(10,7,[49,56,63]);
        pos     = get(axLpm,'position');
        pos(1)  = pos(1)+shiftFraction/2*pos(1);
        set(axLpm,'position',pos)
        
        % Bivariate likelihood function L(m,r)
        axL     = subplot(10,7,[46:48,53:55,60:62]);
        pos     = get(axL,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axL,'position',pos)
        
        % Distance profile likelihood function Lpr
        axLpr = subplot(10,7,[67:69]);
        pos     = get(axLpr,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axLpr,'position',pos)
        
        
        
        %Plot wforms .........................................................
        zFullName = zTarget.fullName{1};    % Get non-available waveforms from bigstar?
        if ( (~exist(zFullName,'file')) && (o.scp_wforms) )
            slash_idx        = regexp(zFullName,'/');
            pathName         = zFullName(1:slash_idx(end));
            [recordName]     = get_recordName(zFullName);
            threeCompPattern = strcat([pathName,recordName,'*']);
            scp_wform(threeCompPattern)
        end
        
        [s,meta] = read_any_trace(zFullName,zTarget,optRead);
        
        ns     = numel(s.raw);
        sr     = meta.sr;
        t      = (1:ns)/sr+zTarget.ts(1);
        ppxIdx = zTarget.ppxIdx(1);
        tppx   = zTarget.tppx(1);
        %tspx   = zTarget.tspx(1);  % Instead of tspx, give it end of window
        t2     = snippet*snpLength;
        %[~,noise] = get_noise(s.vel,ppxIdx,sr,sWindow,nWindow);
        [~,noiseAmps,~,~] = get_snr(s.vel,ppxIdx,sr,sWindow,nWindow);        
        noise             = prctile(noiseAmps,84.1);
        
        nspersnp                      = round(snpLength*sr);
        [amax,amaxIdx,~,nbNoise,sout] = fbank_1trace(s.vel,ppxIdx,sr,0,20,nspersnp,fc,fMode,fOrder);

        sIdx = round(ppxIdx - T*sr*tFrac);      % t-index of window-start
        eIdx = round(ppxIdx + T*sr*(1-tFrac));  % t-index of window-end
        plot_nbWforms_scec(s.vel,sout,amaxIdx,t,tppx,tppx+t2,noise,nbNoise,sIdx,eIdx,fc,'',tFrac,T,axWf,snippet,ftSize);
        
        
        
        % Spectra plot ........................................................
        subplot(axSpect); hold on; grid on
        atmp = aVal(idx_std(1:nsim),:)';
        
        p3 = plot(atmp);
        plot(az_obs,'-w','lineWidth',8);
        plot(az_obs,'-k','lineWidth',6);
        %p4 = plot(az_obs,'-r','lineWidth',3);
        p4 = plot(az_obs,'xy','lineWidth',2,'markerSize',9);
        set(gca,'yscale','log','ylim',[1e-7 1e0],'fontSize',ftSize)
        ylabel('PGV^{nb} [m/s]','fontSize',ftSize)
        xlabel('Frequency bands [Hz]','fontSize',ftSize)
        l1=legend([p3(2);p4],['PGV^{nb} of ',num2str(nsim),' most similar traces'],'PGV^{nb} of target trace');set(l1,'fontSize',ftSize)
        set(l1,'location','southWest')
        
        % Customise annotation
        xT     = cell(9,1);
        for iband = 1:9
            xT{iband} = strcat([num2str(fc(iband,1),'%3.1f'),'-',num2str(fc(iband,2),'%3.1f')]);
        end
        x  = 1:9;
        xT([2,4,6,8]) = {''};
        set(gca,'xtick',x,'xtickLabel',xT,'xlim',[1,9])
        
        % Make upper end lie where wform plot has upper end
        posSpect    = get(axSpect,'position');
        posSpect(4) = posSpect(4)*shrinkFactor;
        posWform    = get(axWf(1),'position');
        y           = posWform(2)+posWform(4)-posSpect(4);
        posSpect(2) = y;
        set(axSpect,'position',posSpect)
        
        
        
        % Plot data ...........................................................
        subplot(axDat); hold on; grid on
        scatter(rVal,mVal)
        l1=line([0 2],                        [m_ctlg m_ctlg],'color','w');
        line([log10(r_ctlg) log10(r_ctlg)],[mmin mmax],    'color','w');
        p2=plot(rHat_noR,mHat_noR,'pw','markerSize',18,'markerFaceColor','b','lineWidth',2);
        p1=plot(rHat,mHat  ,'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
        lgnd=legend([p1;p2;l1],'MLE','MLE wo/ r0','Catalog values');
        if (m_ctlg>5.5); set(lgnd,'location','southeast','fontSize',ftSize); end
        whitebg(gcf,'k'); set(gcf,'InvertHardcopy','off')
        %title('Magnitude and distance values of traces with most similar "spectra"','fontSize',ftSize)
        xlabel(['Hypocentral Distance [km]'],'fontSize',ftSize)
        ylabel(['Magnitude'],'fontSize',ftSize)
        set(gca,'ylim',[mmin mmax],'xlim',[0 2],'fontSize',ftSize)
        rlab = [5,10,20,50,100];
        set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
        posDat = get(axDat,'position');
        posDat(4) = posDat(4)*shrinkFactor;
        posDat(2) = posDat(2)+.1*posDat(4);
        set(axDat,'position',posDat)
        
        
        
        
        % Distance profile likelihood function Lpr ............................
        subplot(axLpr); hold on
        plot(rr,Lprn);
        
        if o.useRconstraint
            Lr0 = 10.^lr(1,:);
            Lr0 = Lr0./sum(Lr0);
            plot(rr,Lr0,':r')
        end
        
        ylabel('Profile likelihood','fontSize',ftSize)
        set(axLpr,'fontSize',ftSize)
        set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
        xlabel(['Hypocentral Distance [km]'],'fontSize',ftSize)
        
        ylm = get(axLpr,'ylim');
        line([log10(r_ctlg) log10(r_ctlg)],[0 ylm(2)],'color','w')
        plot(rHat,max(Lprn),'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
        
        posWf     = get(axWf(end),'position');
        posy      = posWf(2);
        posLpr    = get(axLpr,'position');
        yshift    = posy-posLpr(2);
        posLpr(2) = posy;
        set(axLpr,'position',posLpr)
        
        
        
        % Magnitude profile likelihood function Lpm ...........................
        subplot(axLpm); hold on
        plot(mm,Lpmn);
        ylabel('Profile likelihood','fontSize',ftSize)
        %plot(mm,lpmn);
        view([-90 90])
        set(axLpm,'XAxisLocation','top','Ydir','reverse','fontSize',ftSize)
        xtl = get(axLpm,'xTickLabel');
        xt  = get(axLpm,'xTick');
        xt  = xt(1:end-1);
        xtl = xtl(1:end-1);
        
        ylm = get(axLpm,'ylim');
        line([m_ctlg m_ctlg],[0 ylm(2)],'color','w')
        plot(mHat,max(Lpmn),'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
        
        posLpm = get(axLpm,'position');
        posLpm(2) =posLpm(2)+yshift;
        posLpm(3) =posLpr(4);
        posLpm(4) =posLpm(4)*shrinkFactor;
        set(axLpm,'position',posLpm,'yAxisLocation','Right','xTick',xt,'xTickLabel',xtl)
        
        
        
        
        % Plot bivariate likelihood function L(m,r) ...........................
        subplot(axL); hold on;
        posbak = get(axL,'position');
        %scatter(rVal,mVal,'.w')
        lrel            = l-max(l(:));
        clr             = double(lrel);
        clr(lrel<level) = level;
        surf(rr,mm,double(lrel),clr); view([0 90]);
        hc = colorbar('location','northoutside');
        %set(axL,'position',posbak,'xtickLabel','','fontSize',ftSize)
        title('Relative log-posterior','fontSize',ftSize)
        ylabel(['Magnitude'],'fontSize',ftSize)
        
        l1=line([0 2],                     [m_ctlg m_ctlg],'color','w');
        line([log10(r_ctlg) log10(r_ctlg)],[mmin mmax],    'color','w');
        
        %p1 = plot(rr(colIdx),mm        ,'-ow','markerFaceColor','k');
        %p2 = plot(rr        ,mm(rowIdx),'-ow','markerFaceColor','w','markerEdgeColor','k');
        %lgnd = legend([p1,p2],'P_m(m)','P_m(r)');
        %set(lgnd,'location','southeast','fontSize',ftSize);
        set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
        
        posbak(2) = posbak(2)+yshift;
        posbak(4) = posbak(4)*shrinkFactor;
        set(axL,'position',posbak,'fontSize',ftSize)
        
        
        
        
        
    else
        
        % ----------------------------------------------
        % PLOT FIGURE FOR POWER POINT (LANDSCAPE FORMAT)
        % ----------------------------------------------
        
        % Waveform plots
        axWf = zeros(3,1);
        for i=1:3
            ii = (i-1)*7+1;
            axWf(i) = subplot(7,7,ii:ii+2);
        end
        
        % Spectra plot
        axSpect = subplot(7,7,[22:24,29:31,36:38,43:45]);
        %axSpect = subplot(7,7,[22:24,29:31,36:38]);
        pos     = get(axSpect,'position');
        %pos(1)  = pos(1)+shiftFraction*pos(1);
        %set(axSpect,'position',pos)
        
        % Data plot
        axDat   = subplot(7,7,[4:6,11:13,18:20]);
        pos     = get(axDat,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axDat,'position',pos)
        
        % Magnitude profile likelihood function Lpm
        axLpm   = subplot(7,7,[28,35,42]);
        pos     = get(axLpm,'position');
        pos(1)  = pos(1)+shiftFraction/2*pos(1);
        set(axLpm,'position',pos)
        
        % Bivariate likelihood function L(m,r)
        axL     = subplot(7,7,[25:27,32:34,39:41]);
        pos     = get(axL,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axL,'position',pos)
        
        % Distance profile likelihood function Lpr
        axLpr = subplot(7,7,[46:48]);
        pos     = get(axLpr,'position');
        pos(1)  = pos(1)+shiftFraction*pos(1);
        set(axLpr,'position',pos)
        
        
        
        %Plot wforms .........................................................
        zFullName = zTarget.fullName{1};    % Get non-available waveforms from bigstar?
        if ( (~exist(zFullName,'file')) && (o.scp_wforms) )
            slash_idx        = regexp(zFullName,'/');
            pathName         = zFullName(1:slash_idx(end));
            [recordName]     = get_recordName(zFullName);
            threeCompPattern = strcat([pathName,recordName,'*']);
            scp_wform(threeCompPattern)
        end
        
        [s,meta] = read_any_trace(zFullName,zTarget,optRead);
        
        ns     = numel(s.raw);
        sr     = meta.sr;
        t      = (1:ns)/sr+zTarget.ts(1);
        ppxIdx = zTarget.ppxIdx(1);
        tppx   = zTarget.tppx(1);
        %tspx   = zTarget.tspx(1);  % Instead of tspx, give it end of window
        t2     = snippet*snpLength;
        fprintf(1,'noiseAmps dont make sense. Modify script to use trList.noiseAmps\n')
        [~,noiseAmps,~,~] = get_snr(s.vel,ppxIdx,sr,px.Param.signalWindow,px.Param.gapWindow,px.Param.noiseWindow);
        noise             = prctile(noiseAmps,84.1);
        
%         nspersnp                      = round(snpLength*sr);
%         [amax,amaxIdx,~,nbNoise,sout] = fbank_1trace(s.vel,ppxIdx,sr,0,20,nspersnp,fc,fMode,fOrder);
        [amax,~,amaxIdx,sout,ts] = measure_tdpa(s.vel,ppxIdx,sr,0,10000,snpLength,fc);

        sIdx = round(ppxIdx - T*sr*tFrac);      % t-index of window-start
        eIdx = round(ppxIdx + T*sr*(1-tFrac));  % t-index of window-end
        plot_nbWforms_eewUCB(s.vel,sout,amaxIdx,ts,tppx,tppx+t2,noise,sIdx,eIdx,fc,'',tFrac,T,axWf,snippet,ftSize);
        %plot_nbWforms_eewUCB(s.vel,sout,amaxIdx,t,tppx,tppx+t2,noise,nbNoise,sIdx,eIdx,fc,'',tFrac,T,axWf,snippet,ftSize);
        xlabel('Time since origin time [sec]','fontSize',ftSize)
        posWf  = get(axWf(3),'position');
        % plot_wform(s.vel,meta.t,tppx,zTarget.tspx,99919)
        
        % Spectra plot ........................................................
        subplot(axSpect); hold on; grid on;  box on;
        atmp = aVal(idx_std(1:nsim),:)';
        
        p3 = plot(atmp,'lineWidth',1);
        %plot(az_obs,'-w','lineWidth',7);
        %plot(az_obs,'-k','lineWidth',10);
        %p4 = plot(az_obs,'-r','lineWidth',3);
%        p4 = plot(az_obs,'xr','lineWidth',3,'markerSize',12);
        p4 = plot(az_obs,'xw','lineWidth',6,'markerSize',18);
        p4 = plot(az_obs,'xr','lineWidth',2,'markerSize',12);
        set(gca,'yscale','log','ylim',[1e-7 1e0],'fontSize',ftSize)
        ylabel('PGV^{nb} [m/s]','fontSize',ftSize)
        xlabel('Frequency bands [Hz]','fontSize',ftSize)
        l1=legend([p3(2);p4],['PGV^{nb} of ',num2str(nsim),' most similar traces'],'PGV^{nb} of target trace');set(l1,'fontSize',ftSize)
        set(l1,'location','southWest')
        
        % Customise annotation
        xT     = cell(9,1);
        for iband = 1:9
            xT{iband} = strcat([num2str(fc(iband,1),'%3.1f'),'-',num2str(fc(iband,2),'%3.1f')]);
        end
        x  = 1:9;
        xT([2,4,6,8]) = {''};
        set(gca,'xtick',x,'xtickLabel',xT,'xlim',[1,9])
        
        % Make upper end be somewhat lower
        posSpect    = get(axSpect,'position');
        posSpect(1) = posWf(1);
        posSpect(4) = posSpect(4)*shrinkFactor;
        posSpect(2) = posSpect(2)-0.05;
        set(axSpect,'position',posSpect)
        
        
        
        % Scatter plot of data ............................................
        subplot(axDat); hold on; grid on;
        scatter(rVal,mVal,'markerEdgeColor','k')
        l1=line([0 2],                     [m_ctlg m_ctlg],'color','b','lineWidth',2);
        line([log10(r_ctlg) log10(r_ctlg)],[mmin mmax],    'color','b','lineWidth',2);
        p2=plot(rHat_noR,mHat_noR,'dr','markerSize',12,'markerFaceColor','w','lineWidth',2);
        %    p1=plot(rHat,mHat  ,'pw','markerSize',18,'markerFaceColor','r','lineWidth',2);
        p1=plot(rHat,mHat  ,'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
        lgnd=legend([p1;p2;l1],'MLE','MLE wo/ distance constraint','Catalog values');
        if (m_ctlg>5.5); set(lgnd,'location','southeast','fontSize',ftSize); end
        %whitebg(gcf,'k'); set(gcf,'InvertHardcopy','off')
        %title('Magnitude and distance values of traces with most similar "spectra"','fontSize',ftSize)
        xlabel(['Hypocentral Distance [km]'],'fontSize',ftSize)
        ylabel(['Magnitude'],'fontSize',ftSize)
        set(gca,'ylim',[mmin mmax],'xlim',[0 2],'fontSize',ftSize)
        rlab = [5,10,20,50,100];
        set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
        posDat = get(axDat,'position');
        posDat(2) = posWf(2);
        set(axDat,'position',posDat)
        
        % Distance profile likelihood function Lpr ............................
        subplot(axLpr); hold on;  box on;
        %plot(rr,Lprn,'color','k');
        plot(rr,exp(lmrn),'color','k');
        
        if o.useRconstraint
            Lr0 = 10.^lr(1,:);
            Lr0 = Lr0./sum(Lr0);
            plot(rr,Lr0,'-.r')
        end
        
        ylabel('P_m(R)','fontSize',ftSize)
        %plot(rr,lprn);
        set(axLpr,'fontSize',ftSize,'ylim',[0 .4])
        set(gca,'xtick',log10(rlab),'xtickLabel',rlab,'fontSize',ftSize)
        xlabel(['Hypocentral Distance [km]'],'fontSize',ftSize)
        
        ylm = get(axLpr,'ylim');
        line([log10(r_ctlg) log10(r_ctlg)],[0 ylm(2)],'color','b','lineWidth',2)
        %plot(rHat_noR,max(Lprn),'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
        plot(rHat,max(Lprn),'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);
%         posWf     = get(axWf(end),'position');
%         posy      = posWf(2);
%         posLpr    = get(axLpr,'position');
%         yshift    = posy-posLpr(2);
%         posLpr(2) = posy;
%         set(axLpr,'position',posLpr)
%         
        posLpr    = get(axLpr,'position');
        posLpr(2) = posSpect(2); 
        set(axLpr,'position',posLpr)
        
        
        % Magnitude profile likelihood function Lpm ...........................
        subplot(axLpm); hold on; box on;
        %plot(mm,Lpmn,'color','k');
        plot(mm,Lmmn,'color','k');
        ylabel('P_m(M)','fontSize',ftSize)
        xlabel('Magnitude','fontSize',ftSize)
        %plot(mm,lpmn);
        view([-90 90])
        set(axLpm,'XAxisLocation','top','ylim',[0 .4],'Ydir','reverse','fontSize',ftSize)
        xtl = get(axLpm,'xTickLabel');
        xt  = get(axLpm,'xTick');
        xt  = xt(1:end-1);
        xtl = xtl(1:end-1);
        
        ylm = get(axLpm,'ylim');
        line([m_ctlg m_ctlg],[0 ylm(2)],'color','b','lineWidth',2)
        plot(mHat,max(Lpmn),'pr','markerSize',18,'markerFaceColor','w','lineWidth',2);

        % Make it end where spectra plot ends
        dp        = 0.03;                                                  % Distance between profile and bivariate llh functions
        posLpm    = get(axLpm,'position');
        y_Lpm     = posLpr(2)+posLpr(4)+dp;
        h_Lpm     = posSpect(2)+posSpect(4)-y_Lpm;
        posLpm    = [posLpm(1) y_Lpm 0.8*posLpm(3) h_Lpm];
        set(axLpm,'position',posLpm)
       
        
        
        
        % Plot bivariate likelihood function L(m,r) ...........................
        subplot(axL); hold on; box on;
        posbak = get(axL,'position');
        %scatter(rVal,mVal,'.w')
        lrel            = l-max(l(:));
        clr             = double(lrel);
        clr(lrel<level) = level;
        surf(rr,mm,double(lrel),clr); view([0 90]);
        %hc = colorbar('location','northoutside');
        %set(axL,'position',posbak,'xtickLabel','','fontSize',ftSize)
        %title('Relative log-likelihood','fontSize',ftSize)
        ylabel(['Magnitude'],'fontSize',ftSize)
        
        l1=line([0 2],                     [m_ctlg m_ctlg],'color','w','lineWidth',2);
        line([log10(r_ctlg) log10(r_ctlg)],[mmin mmax],    'color','w','lineWidth',2);
        
        %p1 = plot(rr(colIdx),mm        ,'-ok','markerFaceColor','m');
        %p2 = plot(rr        ,mm(rowIdx),'-ow','markerFaceColor','w','markerEdgeColor','k');
        %lgnd = legend([p1,p2],'P_m(M)','P_m(R)');
        %set(lgnd,'location','southeast','fontSize',ftSize);
        set(gca,'xtick',log10(rlab(1:end-1)),'xtickLabel',rlab(1:end-1),'fontSize',ftSize)
        
        % Make it end where spectra plot ends
        posbak(2) = posLpm(2);
        posbak(4) = posLpm(4);
%        posbak(4) = posSpect(4)-posSpect(2);
        set(axL,'position',posbak)
        %posbak(2) = posbak(2)+yshift;
        %posbak(4) = posbak(4)*shrinkFactor;
        %set(axL,'position',posbak,'fontSize',ftSize)
        
    
        % COLORBAR
        hc = colorbar('northOutside');
        %posScatter = get(axDat,'position');
        posC = get(hc,'position');
        cy   = posbak(2)+posbak(4)+2*dp;
        %cw   = 0.2;
        set(hc,'position',[posC(1) cy posC(3) .015 ],'fontSize',ftSize, ...
            'xAxisLocation','bottom')
        xlabel(hc,'Relative log-posterior','fontSize',ftSize)
        
    end
     
    
    
    if o.printSSR
        figDir = strcat(['~/programs/filterBank/fig/i',num2str(iN),'/methodExamples/new/']);
        if ~o.ppt
            figName=strcat([figDir,'m',num2str(zTarget.m,'%3.1f'),'_',num2str(zTarget.hypDist,'%3.0f'), ...
                'km_',num2str(snippet,'%3.1f'),'snp']);
        else
            figName=strcat([figDir,'m',num2str(zTarget.m,'%3.1f'),'_',num2str(zTarget.hypDist,'%3.0f'), ...
                'km_',num2str(snippet,'%3.1f'),'snp_ppt']);
        end
        
        set(gcf,'PaperPositionMode','auto')
        %set(gcf,'PaperPositionMode','manual')
        %set(gcf,'DefaultFigureColor','remove')
        %set(gcf,'DefaultFigureColor',[1 1 1])
        
        print('-depsc2',[figName,'.eps'])
        print('-dpng',[figName,'.png'])
    end
    1+1;
end
    





%     % Read and plot waveform
%     subplot(4,1,1); hold on; grid on
%     [zS,~,meta] = read_any_trace_proc(zFullName,zTarget);
%     zPxIdx      = meta.ppxIdx;
%     tz          = meta.t;
%     ztpx        = tz(zPxIdx);
%     plot(tz,zS.vel,'lineWidth',1)
%     ymax = 1.1*max(abs(zS.vel));
%     line([ztpx ztpx],[-ymax ymax],'Color','r','lineWidth',2)
%     xlabel('Time [sec]','fontSize',ftSize)
%     ylabel('Velocity [m/s]','fontSize',ftSize)
%     set(gca,'fontSize',ftSize)
%     set(gca,'ylim',[-ymax ymax])
%     set(gca,'xlim',[min(tz) max(tz)])



% .........................................................................
    function [mHat,rHat,lpmn,lprn] = max_llh_function(l)
        
        % Profile likelihood for magnitude Lp(m)
        [lpm,colIdx] = max(l,[],2);
        Lpm          = exp(lpm);
        Lpmn         = Lpm/sum(Lpm);
        [~,idx]      = max(Lpmn);
        lpmn         = log(Lpmn);
        if numel(idx)>1
            fprintf(1,'8ung: non-unique maximum, using middle value. Not sure if ok.\n')
            idx = floor(median(idx));
        end
        mHat         = mm(idx);
        
        
        % Profile likelihood for distance Lp(r)
        [lpr,rowIdx] = max(l,[],1);
        Lpr          = exp(lpr);
        Lprn         = Lpr/sum(Lpr);
        [~,idx]      = max(Lprn);
        lprn         = log(Lprn);
        if numel(idx)>1
            fprintf(1,'8ung: non-unique maximum, using middle value. Not sure if ok.\n')
            idx = floor(median(idx));
        end
        rHat         = rr(idx);
    end



% .........................................................................
    function [mHat,rHat,lmmn,lmrn] = max_llh_function_marginal(l)
        
        L       = exp(l);
        
        % Marginal pdf for magnitude Lm(m)
        dr      = rr(2)-rr(1);          
        Ltmp    = cumsum(L,2)*dr;   % Intrate magnitude out
        Lmm     = Ltmp(:,end);
        Lmmn    = Lmm/sum(Lmm);
        [~,idx] = max(Lmmn);
        lmmn    = log(Lmmn);
        if numel(idx)>1
            fprintf(1,'8ung: non-unique maximum, using middle value. Not sure if ok.\n')
            idx = floor(median(idx));
        end
        mHat    = mm(idx);
        
        
        % Marginal pdf for distance Lm(r)
        dm      = mm(2)-mm(1);          % Intrate magnitude out
        Ltmp    = cumsum(L,1)*dm;
        Lmr     = Ltmp(end,:);
        Lmrn    = Lmr/sum(Lmr);
        [~,idx] = max(Lmrn);
        lmrn    = log(Lmrn);
        if numel(idx)>1
            fprintf(1,'8ung: non-unique maximum, using middle value. Not sure if ok.\n')
            idx = floor(median(idx));
        end
        rHat    = rr(idx);

    end
end