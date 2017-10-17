function analyze_pphys(expNum,taskFig,logConvert)
    if nargin < 1
        expNum = 1;
    else
    end
    if nargin < 2
        taskFig = true;
    else
    end
    if nargin < 3
        logConvert = true;
    else
    end
        
    
    %% add paths
    close all;
    codeFolder = '/Users/kohler/code';
    rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
    setenv('DYLD_LIBRARY_PATH','')

    %% PRELIMINARY
    % folder names
    topPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
    dataPath = sprintf('%s/pphys%0.0f_data',topPath,expNum);
    subPaths = subfolders(sprintf('%s/*20*',dataPath),1);
    % NOTE! should be path corresponding to blank
    rcaPath = sprintf('%s/figures/exp1',topPath);
    savePath = sprintf('%s/figures/pphys%0.0f',topPath,expNum);
    % string for saving the data
    saveStr = datestr(clock,26);
    saveStr(strfind(saveStr,'/')) ='';
    if ~exist(savePath,'dir')
        mkdir(savePath);
    else
    end    
    saveFileName = sprintf('%s/rcaData_%s.mat',savePath, saveStr); % include the date as a string;
    
    % figure labels
    bothLabels = {'Hori2D','Hori3D','Vert2D','Vert3D'};
    condLabels = {'HoriA2D','HoriD2D','HoriA3D','HoriD3D','VertA2D','VertD2D','VertA3D','VertD3D'};
    if expNum == 1
        sweepValues =  [.1600, .217726, .296280, .403175, .548636, .746579, 1.015937, 1.382476, 1.881260, 2.5600];
        falseOrder = {'HoriA2D','HoriA3D','VertA2D','VertA3D', ...
                        'HoriD2D','HoriD3D','VertD2D','VertD3D', ...
                        'VertD3D','VertD2D','HoriD3D','HoriD2D', ...
                        'VertA3D','VertA2D','HoriA3D','HoriA2D'};
        rcaLabels = {'Asc. rel-Mot','Desc. rel-Mot','Asc. rel-Disp','Desc. rel-Disp'};
        pphysLabels = {'rel-Mot','rel-Disp'};
    elseif expNum == 2
        sweepValues = [0.5000,0.7349,1.0801,1.5874,2.3331,3.4290,5.0397,7.4070,10.8863,16.0000];
        % reordering for experiment 2
        falseOrder = {'HoriA2D','HoriA3D','VertA2D','VertA3D', ...
                        'catch','catch','catch','catch', ...
                        'VertD3D','VertD2D','HoriD3D','HoriD2D', ...
                        'catch','catch','catch','catch'};
        rcaLabels = {'Asc. abs-Mot','Desc. abs-Mot','Asc. abs-Disp','Desc. abs-Disp'};
        pphysLabels = {'abs-Mot','abs-Disp'};
    end
    [lia,reorderIdx]=ismember(falseOrder,[condLabels,'catch']); % capture catch trials from experiment 2 as 9th condition

    % rca options
    dataType = 'RLS';
    rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
    binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
    freqsToUse= 1:4; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
    trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
    nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
    nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
    trialReps = [40,24]; % number of trials per condition, 40 in exp1, 24 in exp2.
    %% LOAD IN PPHYS DATA
    
    numSubs = length(subPaths);
    
    for s = 1:numSubs;
        tempPath = subfolders(subPaths{s},1);
        tempPath = tempPath(cellfun(@(x) isempty(strfind(x,'mff')), tempPath));
        matFiles = subfiles(sprintf('%s/Exp_MATL_HCN_128_Avg/RT*',tempPath{end}),1);
        blockNum = 0;
        for m = 1:length(matFiles)
            tmpData = load(matFiles{m});
            if size(tmpData.TimeLine,1) > 1;
                blockNum = blockNum + 1;
                if s == 1 && blockNum == 1
                    conditions = unique(cat(1,tmpData.TimeLine.cndNmb));
                    timingInfo = tmpData.CndTiming(1);
                    coreDuration = timingInfo.stepDurSec * timingInfo.nmbTrialSteps;
                end
                
                numTrials = size(tmpData.TimeLine,1); % trials per block, should not vary from block to block
                if numTrials ~= trialReps(expNum)
                    %  descending catch trial sneaks in among ascending
                    condIdx = cat(1,tmpData.TimeLine.cndNmb);
                    if min(condIdx) == 1 % if ascending
                        tmpData.TimeLine = tmpData.TimeLine(condIdx <= 8); 
                    else
                        tmpData.TimeLine = tmpData.TimeLine(condIdx >= 9);
                    end 
                    numTrials = size(tmpData.TimeLine,1);
                else
                end
                if blockNum == 1
                    trialIdx = 1:numTrials;
                else
                    trialIdx = (1:numTrials)+size(subData,1);
                end
                
                subData(trialIdx,1) = cat(1,tmpData.TimeLine.cndNmb); % condition label
                respIdx = cell2mat(cellfun(@(x) find(ismember({'Mis','Ra','La'},x)),{tmpData.TimeLine.respString},'uni',false))-1; % response (0 = mis, 1 = Ra, 2 = La )
                if ~isempty(respIdx)
                    subData(trialIdx,2) = respIdx;
                else
                    subData(trialIdx,2) = 0;
                end
                subData(trialIdx,3) = cat(1,tmpData.TimeLine.respTimeSec); % response time
                blockTrials(blockNum,s) = numTrials;
                clear tmpData;
            else
            end
        end
        
        % relabel trials
        % in exp1, some n ran 16 conditions, but the last 8 were repeats of the first 8
        % in exp2, all n ran 16, where the other 8 were catch trials
        if max(subData(:,1)) > 8
            subData(:,1) = reorderIdx(subData(:,1));
        else
        end
        
        % pre-generate respData matrix
        subTrials(s) = size(subData,1);
        if ~exist('respData','var')
            respData = NaN(subTrials(s),4,numSubs);
        else
        end
        % and populate with subData
        respData(1:size(subData,1),1:size(subData,2),s) = subData;
        clear subData;
        totalDur = 12; % 10 secs, 2 sec post/pre
        misIdx(:,s) = respData(:,2,s)==0;
        for c=1:length(conditions)
            curIdx = respData(:,1,s) == conditions(c);
            percMis(c,s) = length(find(misIdx(curIdx,s)))./length(find(curIdx));
            respData(curIdx,3,s) = respData(curIdx,3,s)-timingInfo.preludeDurSec; % subtract prelude
            if ~mod(c,2)
                % if ascending
                respData(curIdx,4,s) = sweepConvert( coreDuration,[max(sweepValues),min(sweepValues)],respData(curIdx,3,s) );
            else
                % if descending
                respData(curIdx,4,s) = sweepConvert( coreDuration,[min(sweepValues),max(sweepValues)],respData(curIdx,3,s) );
            end
            if logConvert
                respData(curIdx,4,s) = reallog(respData(curIdx,4,s));
            else
            end
            aveTime(c,s) = mean(respData(~misIdx(:,s) & curIdx,3,s));
            aveStim(c,s) = mean(respData(~misIdx(:,s) & curIdx,4,s));
        end
        if expNum == 2
            catchIdx = respData(:,1,s) == 9;
            catchMis(s) = length(find(misIdx(catchIdx,s)))./length(find(catchIdx));
        else
        end
    end
    
    % exclude subs with more than 15% misses
    includeIdx = mean(percMis) < .15;
    includeSubs = length(find(includeIdx));
   
    aveBoth = (aveStim(1:2:end,:)+aveStim(2:2:end,:))./2;
    bothGrand = nanmean(aveBoth(:,includeIdx),2);
    bothStderr = nanstd(aveBoth(:,includeIdx),0,2)./sqrt(includeSubs);
    stimGrand = nanmean(aveStim(:,includeIdx),2);
    stimStderr = nanstd(aveStim(:,includeIdx),0,2)./sqrt(includeSubs);

    subPaths = subPaths(includeIdx);
    
    %% PLOT PPHYS
    lWidth = 1.5;
    fSize = 12;
    if logConvert
        baseVal = reallog(.1);
        yMin = reallog(.1); yMax = reallog(10);
        logOpts = {'ytick',reallog([.1,.25,.5,1,2,10]),'yticklabel',[.1,.25,.5,1,2,10]};
    else
        baseVal = 0;
        yMin = 0; yMax = 8;
        logOpts ={'ytick',yMin:1:yMax};
    end
    gcaOpts = {'tickdir','out','ticklength',[0.025,0.025],'box','off','fontsize',fSize,'fontname','Arial','linewidth',lWidth};
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    close all
    pphysFig = figure;
    subplot(1,2,1);
    hold on
    for z = 1:8
        if isempty(strfind(condLabels{z},'3D'));
            h2D = bar(z,stimGrand(z),'facecolor',color1(1,:),'edgecolor','none','basevalue',baseVal);
        else
            h3D = bar(z,stimGrand(z),'facecolor',color2(1,:),'edgecolor','none','basevalue',baseVal);
        end
    end
    errorb(1:8,stimGrand,stimStderr);
    set(gca,gcaOpts{:},'xtick',[2.5,6.5],'xticklabel',{'Horizontal','Vertical'},logOpts{:});
    
    %legend([h2D,h3D],pphysLabels,'location','northwest','fontsize',fSize,'fontname','Arial');
    %legend boxoff;
    plot(1:2:length(conditions),ones(1,length(conditions)/2).*-1.5,'k^','LineWidth',lWidth,'markerfacecolor',[1 1 1],'markeredgecolor','none','MarkerSize',10);
    plot(2:2:length(conditions),ones(1,length(conditions)/2).*-1.5,'kv','LineWidth',lWidth,'markerfacecolor',[1 1 1],'markeredgecolor','none','MarkerSize',10);

    hold off
    xlim([.5,8.5]);
    ylim([yMin,yMax]);
    ylabel('Displacement (arcmins)','fontsize',fSize,'fontname','Arial')
    subplot(1,2,2);
    hold on
    for z = 1:4
        if isempty(strfind(bothLabels{z},'3D'));
            h2D = bar(z,bothGrand(z),'facecolor',color1(1,:),'edgecolor','none','basevalue',baseVal);
        else
            h3D = bar(z,bothGrand(z),'facecolor',color2(1,:),'edgecolor','none','basevalue',baseVal);
        end
    end
    errorb(1:4,bothGrand,bothStderr);
    xlim([.5,4.5]);
    ylim([yMin,yMax]);
    set(gca,gcaOpts{:},'xtick',[1.5,3.5],'xticklabel',{'Horizontal','Vertical'},logOpts{:});
    legend([h2D,h3D],pphysLabels,'location','northwest','fontsize',fSize,'fontname','Arial');
    legend boxoff;
    hold off
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = 10;
    figPos(3) = 20;
    set(gcf,'pos',figPos);
    ylabel('Displacement (arcmins)','fontsize',fSize,'fontname','Arial');

    export_fig(sprintf('%s/pphys%0.0f_beh.pdf',savePath,expNum),'-pdf','-transparent',pphysFig);
    %close(pphysFig);

    %% LOAD IN EEG DATA

    for s = 1:length(subPaths)
        tempPath = subfolders(subPaths{s},1);
        tempPath = tempPath(cellfun(@(x) isempty(strfind(x,'mff')), tempPath));
        sourceDataFileName = sprintf('%s/Exp_TEXT_HCN_128_Avg/sourceData_%s.mat',tempPath{end},dataType);
        if isempty(dir(sourceDataFileName))
            createSourceDataMat(sprintf('%s/Exp_TEXT_HCN_128_Avg',tempPath{end}));
        end
        [signalData,noise1,noise2,subFreqIdx{s},subBinIdx{s},subFreqLabels{s},subBinLabels{s}] = selectDataForTraining(sourceDataFileName,binsToUse,freqsToUse,[],trialsToUse);
        if size(signalData,1) > 8
            for c = 1:length(conditions)
                if expNum == 1
                    % reorder the conditions and trials
                    trialOrder = [1:5,11:15,6:10]; % conditions 9 to 16 were run in the middle
                    sensorData{c,s}(:,:,trialOrder) = cat(3,signalData{reorderIdx==c});
                    cellNoiseData1{c,s}(:,:,trialOrder) = cat(3,noise1{reorderIdx==c});
                    cellNoiseData2{c,s}(:,:,trialOrder) = cat(3,noise2{reorderIdx==c});
                else
                    % reorder the conditions, ignore catch trials
                    sensorData{c,s} = cat(3,signalData{reorderIdx==c});
                    cellNoiseData1{c,s} = cat(3,noise1{reorderIdx==c});
                    cellNoiseData2{c,s} = cat(3,noise2{reorderIdx==c});
                end
            end
        else
            sensorData(:,s)=signalData;
            cellNoiseData1(:,s)=noise1;
            cellNoiseData2(:,s)=noise2;
        end
    end
    %% CHECK FREQUENCY AND BIN INDICES FOR CONSISTENCY ACROSS SUBS
    nSubs = size(sensorData,3);
    for s=1:nSubs
        for c = 1:length(conditions)
            if sum(abs(subFreqIdx{s}{c}-subFreqIdx{1}{c}))~=0 && sum(abs(subBinIdx{s}{c}-subBinIdx{1}{c}))~=0
                error('Frequency and bin indices vary across subjects: check consistency of DFT/RLS exports\n.');
            end
        end
    end

    rcaFiles = subfiles([rcaPath,'/rcaData*.mat'],1);
    load(rcaFiles{end},'allRCA');

    %% repopulate RCA struct with pphys data
    warning('off','all')
    rcaData = rcaProject(sensorData,allRCA.W);
    noiseData.lowerSideBand=rcaProject(cellNoiseData1,allRCA.W); 
    noiseData.higherSideBand=rcaProject(cellNoiseData2,allRCA.W);

    %% create a "component" of just one channel for performance evaluation if requested
    chanToCompare = 75;
    nChannels = size(sensorData{1},2);
    wComparison=zeros(nChannels,1); wComparison(chanToCompare)=1; 
    comparisonData=rcaProject(sensorData,wComparison); 
    comparisonNoiseData.lowerSideBand =rcaProject(cellNoiseData1,wComparison); 
    comparisonNoiseData.higherSideBand =rcaProject(cellNoiseData2,wComparison); 

    %% generate final output struct
    allRCA.data = rcaData;
    allRCA.noiseData = noiseData;
    allRCA.comparisonData = comparisonData;
    allRCA.comparisonNoiseData = comparisonNoiseData;
    allRCA.inputData = sensorData;
    allRCA.settings.binLevels = subBinLabels{1};

    %% rca replaces NaNs with zeroes, correct this
    nanDims = [1,2]; % if all time points are zero, or all channels are zero
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};
    for z=1:length(structVars)
        if strfind(lower(structVars{z}),'noise')
            for n = 1:length(noiseVars)
                allRCA.(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}).(noiseVars{n}),'uni',false);
            end
        else   
            allRCA.(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}),'uni',false);
        end
    end
    warning('on','all')

    %% COMPUTE VALUES FOR PLOTTING
    % AE RCA
    keepConditions = true;
    errorType = 'SEM';
    trialError = true;
    tempDataStrct = aggregateData(allRCA.data,allRCA.settings,keepConditions,errorType,trialError);
    ampVals(:,:,1:5,:) = tempDataStrct.ampBins;
    errLB(:,:,1:5,:)   =tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,:,1:5,:)   =tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(allRCA.noiseData.lowerSideBand,allRCA.settings,keepConditions,errorType,trialError);
    tempNoiseStrct2 = aggregateData(allRCA.noiseData.higherSideBand,allRCA.settings,keepConditions,errorType,trialError);
    [snrVals(:,:,1:5,:),noiseVals(:,:,1:5,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);

    % COMPARISON
    tempDataStrct = aggregateData(allRCA.comparisonData,allRCA.settings,keepConditions,errorType,trialError);
    ampVals(:,:,6,:) = tempDataStrct.ampBins;
    errLB(:,:,6,:)=tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,:,6,:)=tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(allRCA.comparisonNoiseData.lowerSideBand,allRCA.settings,keepConditions,errorType,trialError);
    tempNoiseStrct2 = aggregateData(allRCA.comparisonNoiseData.higherSideBand,allRCA.settings,keepConditions,errorType,trialError);
    [snrVals(:,:,6,:),noiseVals(:,:,6,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);

    %% PLOT RCs
    close all;

    condsToUse = [1,1,3,3,1,1,3,3]; % you are plotting 2D rel and 3D rel x 2
    rcaType = 'all'; % full or freq?
    plotSNR = false;

    % plot settings
    % set figure size in the beginning
    figHeight = 30;
    figWidth = 20;

    binVals = cellfun(@(x) str2num(x), allRCA.settings.binLevels{1});

    lWidth = 1.5;
    %red =   [152 53 49; 231 184 182]./255;
    %green = [117 148 54; 216 229 186]./255;
    %blue =  [52 95 148; 183 204 229]./255;
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    subColors = repmat([color1; color2],2,1);
    subColors = subColors(condsToUse,:);
    fSize = 12;
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Arial','linewidth',lWidth};
    for f=1:length(freqsToUse)
        figure;
        curFreq = freqsToUse(f);
        for r = 1:(nComp+1)
            if r<6
                egiH(r) = subplot(nComp+1,3,3+(r-1)*3);
                hold on
                rcaColorBar = [min(allRCA.A(:,r)),max(allRCA.A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(allRCA.A(:,r),rcaColorBar);
                hold off
            else
            end
            for s=1:2 % vertical or horizontal motion
                if s==1
                    subplot(nComp+1,3,1+(r-1)*3);
                    curConds = find(ismember(allRCA.settings.condsToUse,1:4));
                    titleStr = sprintf('horizontal: %s',allRCA.settings.freqLabels{1}{curFreq});
                else
                    subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(allRCA.settings.condsToUse,5:8));
                    titleStr = sprintf('vertical: %s',allRCA.settings.freqLabels{1}{curFreq});
                end
                hold on
                for c=1:length(curConds)
                    if isempty(strfind(condLabels{curConds(c)},'A'))
                        % descending, download triangle
                        ampMarkerStyle = {'-v','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7}; % 'MarkerEdgeColor','none'
                    else
                        % ascending, upward triangle
                        ampMarkerStyle = {'-^','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7};
                    end
                    if plotSNR
                        valSet = snrVals(:,curFreq,r,:);
                        ampH(c)=plot(binVals,snrVals(:,curFreq,r,curConds(c)),ampMarkerStyle{:});
                        plot(binVals,noiseVals(:,curFreq,r,curConds(c)),'sq','Color',subColors(curConds(c),:),'MarkerSize',5);
                    else
                        valSet = ampVals(:,curFreq,r,:);
                        ampH(c)=plot(binVals,ampVals(:,curFreq,r,curConds(c)),ampMarkerStyle{:});
                        %plot(binVals,noiseVals(:,curFreq,r,curConds(c)),'sq','Color',subColors(curConds(c),:),'MarkerSize',5);
                        %eH = errorbar(binVals,ampVals(:,curFreq,r,curConds(c)),errLB(:,curFreq,r,curConds(c)),errUB(:,curFreq,r,curConds(c)),'Color',subColors(curConds(c),:),'LineWidth',lWidth);
                        hE = ErrorBars(binVals,ampVals(:,curFreq,r,curConds(c)),[errLB(:,curFreq,r,curConds(c)),errUB(:,curFreq,r,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false);
                        cellfun(@(x) uistack(x,'bottom'), hE);
                        hold on
                    end
                end
                if f == 2 && r == 1     
                    yUnit = 0.4;
                    yMax = 2.0;
                else
                    yUnit = 0.2;
                    yMax = 1.2;
                end
            	if expNum == 1
                    xMin = .13;
                    xMax = 3;
                else
                    xMin = .4;
                    xMax = 18;
                    yMax = yMax * 2;
                    yUnit = yUnit * 2;
                end
                   
                %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
                %yMax = ceil(max(valSet(:)))+yUnit;
                set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');
          
                xlim([xMin,xMax]);
                ylim([0,yMax])

                % plot noise patch
                meanNoise = max(noiseVals(:,curFreq,r,curConds),[],4)';
                yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
                xNoiseVals = [min(get(gca,'xlim')), min(get(gca,'xlim')),binVals',max(get(gca,'xlim')),max(get(gca,'xlim'))];
                pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
                uistack(pH,'bottom')

                if s==2
                    if r > nComp
                        text(0.0008,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,'OZ','fontsize',fSize,'fontname','Arial');
                    else
                        text(0.0008,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,['RC ',num2str(r)],'fontsize',fSize,'fontname','Arial');
                    end
                else
                end
                if r==1
                    title(titleStr,'fontsize',fSize,'fontname','Arial');
                elseif r==6
                    if s==1
                        ylabel('Amplitude (\muV)')
                        xlabel('Distance (arcmins)');
                    else
                        lH = legend(ampH,rcaLabels,'location','northeast');
                        legend boxoff
                        lPos = get(lH,'position');
                        lPos(1) = lPos(1) + .23;
                        lPos(2) = lPos(2) + .05;
                        set(lH,'position',lPos);
                        tH = text(3.5,-0.1,sprintf('n = %0d',size(allRCA.data,2)),'fontsize',fSize,'fontname','Arial');
                    end
                end
                hold off;
            end
        end
        drawnow;
        for r = 1:5;
            addX = 1.6;
            addY = 1.6;
            newPos = get(egiH(r),'position');
            newPos(1) = newPos(1)-(newPos(3)*addX*.7);
            newPos(2) = newPos(2)-(newPos(4)*addY*.6);
            newPos(3) = newPos(3)*(1+addX);
            newPos(4) = newPos(4)*(1+addY);
            set(egiH(r),'position',newPos);
        end
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = figHeight;
        figPos(3) = figWidth;
        set(gcf,'pos',figPos);
        if plotSNR        
            export_fig(sprintf('%s/pphys%0.0f_rc%d_%s_snr.pdf',savePath,expNum,f,rcaType),'-pdf','-transparent',gcf);
        else
            export_fig(sprintf('%s/pphys%0.0f_rc%d_%s.pdf',savePath,expNum,f,rcaType),'-pdf','-transparent',gcf);
        end
    end
    close all;

%% PLOT TASK ILLUSTRATION
if taskFig
    for z = 1:2
        if z ==1
            plot(sweepValues,'ma^-','markerfacecolor',[1 1 1],'MarkerSize',10,'linewidth',lWidth);
        else
            plot(fliplr(sweepValues),'mav-','markerfacecolor',[1 1 1],'MarkerSize',10,'linewidth',lWidth);
        end
        hold on
        set(gca,gcaOpts{:},'xtick',1:10);
        xlim([0.5,10.5]);
        xlabel('time(s)');
        ylabel('displacement(arcmins)');
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = 6;
        figPos(3) = 8;
        set(gcf,'pos',figPos);
        export_fig(sprintf('%s/pphys%0.0f_task%d.pdf',savePath,expNum,z),'-pdf','-transparent',gcf);
        if z==2
            plot(sweepValues,'ma^-','markerfacecolor',[1 1 1],'MarkerSize',10,'linewidth',lWidth);
            export_fig(sprintf('%s/pphys%0.0f_task_both.pdf',savePath,expNum),'-pdf','-transparent',gcf);
        else
        end
        hold off

        close(gcf);
    end
else
end

