function infant_analysis_m2D3D(projectedData,rcaType,doFreq)
    %% Add Paths
    close all;
    setenv('DYLD_LIBRARY_PATH','')

    
    if nargin < 1
        % PROJECTED DATA OR NOT?
        projectedData = true;
    else
    end
    if nargin < 2
        rcaType = 'freq';
    else
    end
    if nargin < 3
        doFreq = [1,2];
    else
    end

    if numel(doFreq) > 1
        for z = 1:length(doFreq)
            infant_analysis_m2D3D(projectedData,rcaType,doFreq(z));
        end
        return;
    else
    end
    if strcmp(rcaType,'all');
        plotComp = 1;   % plot only 1st RC.
        flipVal = -1;
    else
        if doFreq == 2
            plotComp = 1;   % plot only 1st RC.
            flipVal = 1;
        elseif doFreq == 1
            plotComp = 5;   % plot only 3rd RC.
            flipVal = 1;
        else
            error('freq %0.0f, not sure which RC to plot!',doFreq);
        end
    end

    mainPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
    %mainPath = '/Users/kohler/Desktop';
    figureFolder = sprintf('%s/figures/infant_exp',mainPath);
    load(sprintf('%s/BabyDataOutput.mat', figureFolder),'babyRCA'); %file name from prep workspace export

    %% PLOT RCs

    if any(~cell2mat(arrayfun(@(x) any(ismember([1,2],babyRCA(x).settings.rcaFreqs)), 1:length(babyRCA),'uni',false)))
        error('different frequencies used in different RCA iterations');
    else
        freqsToUse = 1:2;
    end

    if any(any( diff(cell2mat(arrayfun(@(x) babyRCA(x).settings.rcaConds, 1:length(babyRCA),'uni',false)'))>0 ,2))
        error('different conditions used in different RCA iterations');
    else
        condsToUse = babyRCA(1).settings.rcaConds;
    end

    binVals = cellfun(@(x) str2num(x), babyRCA(1).settings.binLabels);
    logStep = diff(reallog(binVals(1:2))); % step size
    extraBins = arrayfun(@(x) exp(reallog(binVals(end))+x), [logStep,logStep*2]);
    xMin = exp(reallog(binVals(1))-logStep*.5);
    xMax = exp(reallog(binVals(end))+logStep*2.5); % add 2.5 steps
    figHeight = 8; % set figure size in the beginning
    figWidth = 16;
    lWidth = 1.5;
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    subColors = repmat([color1; color2],2,1);
    subColors = subColors(condsToUse,:);
    fSize = 12;
    gcaOpts = {'tickdir','out','ticklength',[0.0250,0.0250],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
    condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
    condLabels = condLabels(condsToUse);

    if strcmp(rcaType,'all');
        curFreq = length(freqsToUse)+doFreq;
    else
        curFreq = doFreq;
    end

    figure;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = figHeight;
    figPos(3) = figWidth;
    set(gcf,'pos',figPos);

    % EGI PLOT
    egiH = subplot(1,3,3);
    hold on
    rcaColorBar = [min(babyRCA(curFreq).A(:,plotComp)),max(babyRCA(curFreq).A(:,plotComp))];
    newExtreme = round(max(abs(rcaColorBar))*10)./10;
    rcaColorBar = [-newExtreme,newExtreme*1.001];
    [figH,cH] = mrC.plotOnEgi(babyRCA(curFreq).A(:,plotComp).*flipVal,rcaColorBar,true);

    set(cH,'YTickMode','manual','location','southoutside','units','centimeters')
    set(cH,'fontsize',fSize,'fontname','Helvetica','XTickMode','manual','Xtick',linspace(min(rcaColorBar),min(rcaColorBar)*-1,5));
    xlabel(cH,'weights','fontsize',fSize,'fontname','Helvetica')
    hold off

    eegH = subplot(1,3,1:2);
    curConds = find(ismember(condsToUse,1:4));
    titleStr = sprintf('horizontal: %s',babyRCA(3).settings.freqLabels{doFreq});
    hold on

    NRvals = babyRCA(curFreq).stats.NR_Params;
    NRerrs = babyRCA(curFreq).stats.NR_JKSE;
    NRmodel = babyRCA(curFreq).stats.hModel;

    % compute new signal values, averaged over bins
    [rcaDataReal,rcaDataImag] = getRealImag(babyRCA(curFreq).data);
    rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,plotComp,:),3)),rcaDataReal,'uni',false);
    rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
    rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,plotComp,:),3)),rcaDataImag,'uni',false);
    rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
    realBinMean = squeeze(nanmean(rcaDataReal));
    imagBinMean = squeeze(nanmean(rcaDataImag));
    % compute new noise values, averaged over bins
    [noiseLoReal,noiseLoImag] = getRealImag(babyRCA(curFreq).noiseData.lowerSideBand(curConds,:));
    [noiseHiReal,noiseHiImag] = getRealImag(babyRCA(curFreq).noiseData.lowerSideBand(curConds,:));
    noiseReal = cellfun(@(x,y) (x+y)./2, noiseLoReal,noiseHiReal, 'uni',false);
    noiseImag = cellfun(@(x,y) (x+y)./2, noiseLoImag,noiseHiImag, 'uni',false);
    noiseReal = cellfun(@(x) squeeze(nanmean(x(:,plotComp,:),3)),noiseReal,'uni',false);
    noiseReal = cell2mat(permute(noiseReal,[3,2,1]));
    noiseImag = cellfun(@(x) squeeze(nanmean(x(:,plotComp,:),3)),noiseImag,'uni',false);
    noiseImag = cell2mat(permute(noiseImag,[3,2,1]));
    realBinMeanNoise = squeeze(nanmean(noiseReal));
    imagBinMeanNoise = squeeze(nanmean(noiseImag));

    if ~projectedData
        % grab values from data structure
        valSet = squeeze(babyRCA(curFreq).stats.Amp(:,plotComp,:));
        errSet1 = squeeze(babyRCA(curFreq).stats.ErrLB(:,plotComp,:));
        errSet2 = squeeze(babyRCA(curFreq).stats.ErrUB(:,plotComp,:));
        % compute vector means of mean bins, and errors
        valSet(length(binVals)+1,:) = sqrt(nanmean(realBinMean,1).^2+nanmean(imagBinMean,1).^2);

        % store t-values for later
        rc_tSqrdP = permute(squeeze(babyRCA(curFreq).stats.tSqrdP(:,plotComp,curConds)),[2,1]);
        rc_tSqrdVal = permute(squeeze(babyRCA(curFreq).stats.tSqrdVal(:,plotComp,curConds)),[2,1]);

        % compute elliptical error and do Hotelling's T2 against zero
        for c=1:length(curConds)
            xyBinMean = cat(2,realBinMean(:,c),imagBinMean(:,c));
            nanVals = sum(isnan(xyBinMean),2)>0;
            binErrs = fitErrorEllipse(xyBinMean(~nanVals,:),'SEM');
            errSet1(length(binVals)+1,c) = binErrs(1);
            errSet2(length(binVals)+1,c) = binErrs(2);
            % compute t-values
            tStruct = tSquaredFourierCoefs(xyBinMean(~nanVals,:));
            rc_tSqrdP(c,length(binVals)+1) = tStruct.pVal;
            rc_tSqrdVal(c,length(binVals)+1) = tStruct.tSqrd;
            % note, just using mean df for all values, would need
            % to be fixed if multi data was ever used seriously
            rc_tSqrdDF(c,1:length(binVals)+1) = tStruct.df2;
        end
    else
        % compute projected vector mean
        % move subjects to first dim
        realVector = permute(rcaDataReal,[2,1,3]);
        imagVector = permute(rcaDataImag,[2,1,3]);
        project_amps = vectorProjection(realVector,imagVector);

        % compute mean of projected vector amplitude, over bins
        project_amps(:,end+1,:) = nanmean(project_amps,2);
        % if all values are NaN
        if any(sum(squeeze(all(isnan(project_amps),2)),2)>0)
            nan_subs = find(sum(squeeze(all(isnan(project_amps),2)),2)>0);
            for z = 1:length(nan_subs)
                msg = ...
                    sprintf('Subject %d has no values in one or more conditions', ...
                    nan_subs(z));
                warning(msg);
                %project_amps(nan_subs(z),:,:) = NaN;
            end
        else
        end
        % make new valset and error set
        valSet = squeeze(nanmean(project_amps,1));
        temp_SE = squeeze(nanstd(project_amps,0,1)./sqrt(size(project_amps,1)));
        errSet1 = temp_SE;
        errSet2 = temp_SE;
        % compute t-values
        [~,temp_p,~,temp_stats] = ttest(project_amps,0,'alpha',0.05,'dim',1,'tail','right');
        rc_tSqrdP = permute(temp_p,[3,2,1]);
        rc_tSqrdVal = permute(temp_stats.tstat,[3,2,1]);
        rc_tSqrdDF = permute(temp_stats.df,[3,2,1]);
        clear temp_*;
    end
    % maximum noise across all four conditions being plotted
    % since this is just means, we can compute it the same way for
    % projected and not-projected
    noiseSet = max(babyRCA(curFreq).stats.NoiseAmp(:,plotComp,curConds),[],3);
    % max of mean over bins
    noiseSet(length(binVals)+1) = max(sqrt(nanmean(realBinMeanNoise,1).^2+nanmean(imagBinMeanNoise,1).^2));

    % do paired tests for Naka-Rushton
    testVal = squeeze(babyRCA(curFreq).stats.NR_JKParams(:,:,plotComp,curConds));
    testVal = permute(testVal,[2,1,3]); % move subjects to first dim
    tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp
    paramIdx = [1,2,3,4];% only look at c50 and rMax
    jkDf = size(testVal,1)-1;
    for pT=1:length(tIdx) % do four paired tests
        diffErr = jackKnifeErr(testVal(:,paramIdx,tIdx(pT,1))-testVal(:,paramIdx,tIdx(pT,2)));
        grandDiff = NRvals(paramIdx,plotComp,tIdx(pT,1)) - NRvals(paramIdx,plotComp,tIdx(pT,2));
        paramPairedT(:,pT) = grandDiff'./diffErr;
        paramPairedP(:,pT) = 2*tcdf( -abs(paramPairedT(:,pT)) , jkDf);
    end

    tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp

    for pT=1:length(tIdx) % do four paired tests
        for b = 1:(length(binVals)+1)
            if ~projectedData
                if b <= (length(binVals))
                    xyData = permute(cat(1,rcaDataReal(b,:,tIdx(pT,:)),rcaDataImag(b,:,tIdx(pT,:))),[2,1,3]);
                else
                    xyData = cat(3,realBinMean(:,tIdx(pT,:)), imagBinMean(:,tIdx(pT,:)));
                end
                tempStrct = tSquaredFourierCoefs(xyData);
                tempP(:,b) = tempStrct.pVal;
                tempT(:,b) = tempStrct.tSqrd;
                tempD(:,b) = tempStrct.mahalanobisD;
                tempU(:,b) = tempStrct.cohenNonOverlap;
                tempMu1(:,b) = valSet(b,tIdx(pT,1));
                tempMu2(:,b) = valSet(b,tIdx(pT,2));
                tempDF(:,b) = tempStrct.df2;
            else
                curData = squeeze(project_amps(:,b,tIdx(pT,:)));
                not_nan = ~any(isnan(curData),2);
                [~,tempP(:,b),~,tempStrct] = ttest(curData(not_nan,1),curData(not_nan,2),'alpha',0.05,'dim',1,'tail','both');
                tempT(:,b) = tempStrct.tstat;
                tempD(:,b) = tempT(:,b)./sqrt(tempStrct.df+1); % Cohen's D
                OVL = 2*normcdf(-abs(tempD(:,b))/2);
                tempU(:,b) = 1-OVL/(2-OVL);
                tempMu1(:,b) = mean(curData(not_nan,1));
                tempMu2(:,b) = mean(curData(not_nan,2));
                tempDF(:,b) = tempStrct.df;
                %tempMu1(:,b) = project_valSet(b,tIdx(pT,1));
                %tempMu2(:,b) = project_valSet(b,tIdx(pT,2));
            end
        end
        rc_tSqrdP = cat(1,rc_tSqrdP,tempP);
        rc_tSqrdVal = cat(1,rc_tSqrdVal,tempT);
        rc_tSqrdDF = cat(1,rc_tSqrdDF,tempDF);
        if pT == 1
            rc_tSqrdD = tempD;
            rc_tSqrdU = tempU;
            rc_tSqrdMu1 = tempMu1;
            rc_tSqrdMu2 = tempMu2;
        else
            rc_tSqrdD = cat(1,rc_tSqrdD,tempD);
            rc_tSqrdU = cat(1,rc_tSqrdU,tempU);
            rc_tSqrdMu1 = cat(1,rc_tSqrdMu1,tempMu1);
            rc_tSqrdMu2 = cat(1,rc_tSqrdMu2,tempMu2);
        end
        clear temp*
    end
    
    % do the actual plotting
    for c=1:length(curConds)
        ampH(c) =plot([binVals;extraBins((c>2)+1)],valSet(:,curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
        hE = ErrorBars([binVals;extraBins((c>2)+1)],valSet(:,curConds(c)),[errSet1(:,curConds(c)),errSet1(:,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
        uistack(ampH(c),'bottom')
        cellfun(@(x) uistack(x,'bottom'), hE);
        hold on

        if ~isnan(NRvals(1,plotComp,curConds(c))) % plot Naka-Rushton
            nFine = 1e2;
            nrX = linspace( min(binVals), max(binVals), nFine )';
            nrVals = NRmodel( nrX, NRvals(:,plotComp,curConds(c)));
            nR(c) = plot( nrX, nrVals, '-','color',subColors(curConds(c),:),'LineWidth',lWidth);

        else
        end
    end
    arrayfun(@(x) uistack(x,'bottom'),nR);
    if freqsToUse(doFreq) == 2
        yUnit = 2;
        yMax = 12.0;
    else
        yUnit = .5;
        yMax = 3;
    end
    set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[2,4,8,16,32],'ytick',0:yUnit:yMax,'Layer','top','clipping','off', 'color','none');
    set(gca,'XMinorTick','off');
    xlim([xMin,xMax]);
    ylim([0,yMax])

    % split point between bin and average vals
    xSplit = exp(reallog(binVals(end))+logStep*.5);
    plot(ones(2,1)*xSplit,[0,yMax],'k','LineWidth',lWidth)
    % plot noise patch,
    % add mean bin noise twice to the end with extra bin
    % we are plotting the mean values over two x-axis points
    xNoiseVals = [xMin,xMin,binVals',xSplit,xSplit];
    yNoiseVals = [0,noiseSet(1),noiseSet(1:end-1)',noiseSet(end-1),0]; % start and end points just repeats of first and last
    % additional vals for averages
    xNoiseVals = [xNoiseVals,xSplit,xSplit,extraBins,xMax,xMax];
    yNoiseVals = [yNoiseVals,0,repmat(noiseSet(end),1,4),0];
    pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
    uistack(pH,'bottom')

    text(binVals(1),yMax*0.95,titleStr,'fontsize',fSize,'fontname','Helvetica');
    ylabel('amplitude (\muV)');
    xlabel('displacement (arcmins)');

    text(extraBins(1),yMax*0.95,sprintf('n = %0d',size(babyRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Helvetica'); %(60,11,...) for first 2 args in 2F1 plot; (60,4...) for 1F1

    hold off;
    drawnow;
    addX = 1.5;
    addY = 1.5;
    set(egiH,'units','centimeters');
    if doFreq == 1
        newPos = get(egiH,'position');
        newPos(1) = newPos(1) - newPos(3)*addX*0.48;
        newPos(2) = newPos(2) - newPos(4)*addY*0.52;
        newPos(3) = newPos(3)+newPos(3)*addX;
        newPos(4) = newPos(4)+newPos(4)*addY;
        set(egiH,'position',newPos);
        save(sprintf('%s/EEGpos.mat',figureFolder),'newPos');
    else
        load(sprintf('%s/EEGpos.mat',figureFolder),'newPos');
        set(egiH,'position',newPos);
    end
    cBarPos = get(cH,'position');
    
    xCenter = newPos(1)+newPos(3)/2;
    yCenter = newPos(2)+newPos(4)/2;
    
    cBarPos(3) = newPos(3)/2;%ceil(cBarPos(3));
    cBarPos(4) = cBarPos(3)/10;
    cBarPos(1) = xCenter-cBarPos(3)/2;
    %cBarPos(2) = newPos(2)-cBarPos(4)/2;
    
    set(cH,'position',cBarPos);
    
    plotPos = get(eegH,'position');
    plotPos(4) = plotPos(4)-.1;
    plotPos(2) = plotPos(2)+.1;
    set(eegH,'position',plotPos);

    if projectedData
        export_fig(sprintf('%s/2DInfant_%dF1_proj_%s.pdf',figureFolder,freqsToUse(doFreq),rcaType),'-nocrop','-pdf','-transparent',gcf);
    else
        export_fig(sprintf('%s/2DInfant_%dF1_multi_%s.pdf',figureFolder,freqsToUse(doFreq),rcaType),'-nocrop','-pdf','-transparent',gcf);
    end
    %% PLOT NR
    close all;
    figure;
    for z = 1:length(paramIdx)
        subplot(1,4,z);
        hold on
        plot(1:4,squeeze(NRvals(paramIdx(z),plotComp,:)),'-o')
        errorb(1:4,squeeze(NRvals(paramIdx(z),plotComp,:)),squeeze(NRerrs(paramIdx(z),plotComp,:)));
        hold off
        xlim([0.5,4.5]);
    end
    
    if ~projectedData || doFreq ~= 2
        return
    else
    end
    %% MAKE BIN-BY-BIN TABLE
    
    dfTest = repmat(rc_tSqrdDF(:,1),1,size(rc_tSqrdDF,2))==rc_tSqrdDF;
    if ~all(dfTest(:));
        msg = '\n df differs among bins \n';
        error(msg);
    else
    end
    
    testList =  {'In_RefVsNo', 'Anti_RefVsNo', 'Ref_InVsAnti', 'NoRef_InVsAnti'};
    % make table of paired results
    freqNum = 2;
    finishedArray = [];
    for z = 1:length(testList)
        switch testList{z}
            case 'Ref_InVsAnti'
                testIdx = 3;
                testName = 'Ref: In vs Anti';
            case 'NoRef_InVsAnti'
                testIdx = 4;
                testName = 'UnRef: In vs Anti';
            case 'In_RefVsNo'
                testIdx = 1;
                testName = 'In-phase: Ref';
            case 'Anti_RefVsNo'
                testIdx = 2;
                testName = 'Anti-phase: Ref';
            otherwise
        end 
        stringBinVals = cat(1,arrayfun(@(x) num2str(x,'%0.2f'),binVals,'uni',false),{'n/a'});
        
        numBinP = rc_tSqrdP(4+testIdx,:);
        stringBinP = arrayfun(@(x) num2str(x,'%0.4f'), numBinP,'uni',false);
        sigIdx = cell2mat(arrayfun(@(x) x < 0.0001,numBinP,'uni',false));
        stringBinP(sigIdx) = {'<0.0001'};
        numBinT = rc_tSqrdVal(4+testIdx,:);    
        stringBinT = arrayfun(@(x) num2str(x,'%0.4f'), numBinT,'uni',false);
        numBinD = rc_tSqrdD(testIdx,:);    
        stringBinD = arrayfun(@(x) num2str(x,'%0.4f'), numBinD,'uni',false);
            
        numBinMu1 = rc_tSqrdMu1(testIdx,:);  
        stringBinMu1 = arrayfun(@(x) num2str(x,'%0.4f'), numBinMu1,'uni',false);
        numBinMu2 = rc_tSqrdMu2(testIdx,:); 
        stringBinMu2 = arrayfun(@(x) num2str(x,'%0.4f'), numBinMu2,'uni',false);            
        rowLabels = {'disp (arcmins)'};
        readyArray = [stringBinVals'; stringBinP; stringBinT; stringBinD];
        
        rowLabels = {rowLabels{:},'p','t-statistic','Cohen''s D'};
        readyArray = cat(2,rowLabels',readyArray);
        varLabels = cat(2,{sprintf('%s (df=%0.0f)',testName,rc_tSqrdDF(4+testIdx,end))},...
              arrayfun(@(x) num2str(x,'bin%0.0f'),1:length(binVals),'uni',false),...
              {'ave'});
        readyArray = cat(1,varLabels,readyArray);
        finishedArray = cat(1,finishedArray,readyArray);
    end
    
    infant_file = sprintf('%s/2DInfant_%dF1_proj_%s.csv',figureFolder,freqsToUse(doFreq),rcaType);
    if exist(infant_file,'file')
        delete(infant_file);
    else
    end
    binTable = array2table(finishedArray);
    writetable(binTable,infant_file,'WriteRowNames',false,'WriteVariableNames',false);
end
