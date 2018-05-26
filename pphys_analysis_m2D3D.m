function pphys_analysis_m2D3D(expNum,taskFig,logConvert,rcaType,projectedData,flipEEG,mergeEEG)
    if nargin < 1
        expNum = 2;
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
    if nargin < 4
        rcaType = 'freq';
    else
    end
    if nargin < 5
        projectedData = true;
    else
    end
    if nargin < 6
        flipEEG = true;
    else
    end
    if nargin < 7
        mergeEEG = true;
    else
    end
    
    if ~flipEEG && mergeEEG
        msg = sprintf('\ERROR: nbad idea to merge w/o also flipping\n');
        error(msg);
    else
    end 
    
    %% add paths
    close all;
    setenv('DYLD_LIBRARY_PATH','')

    %% PRELIMINARY
    % folder names
    topPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
    dataPath = sprintf('%s/pphys_exp%0.0f',topPath,expNum);
    subPaths = subfolders(sprintf('%s/*20*',dataPath),1);
    adultExp = [1,2,3,4,5];
    % NOTE! should be path corresponding to blank (2);
    rcaPath = sprintf('%s/figures/exp%0.0f/rcaData.mat',topPath,adultExp(2));
    savePath = sprintf('%s/figures/paper_figures/figure4',topPath);
    % string for saving the data
    saveStr = datestr(clock,26);
    saveStr(strfind(saveStr,'/')) ='';
    if ~exist(savePath,'dir')
        mkdir(savePath);
    else
    end    
    saveFileName = sprintf('%s/rcaData_exp%0.0f_%s.mat',savePath,expNum,saveStr); % include the date as a string;
    
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
                        'HoriD2D','HoriD3D','VertD2D','VertD3D', ...
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
    
    %% DO PPHYS STATS
    testData = exp(aveBoth(:,includeIdx));
    for z = 1:2
        [pphysTest(z).Sig, pphysTest(z).P,ci,testStrct] = ttest(testData(z+(z-1),:),testData(z+z,:),'dim',2,'tail','both');
        pphysTest(z).df = testStrct.df;
        pphysTest(z).tstat = testStrct.tstat;
        factorDiff(z) = exp(bothGrand(z+z))/exp(bothGrand(z+(z-1)));
    end
   
    
    %% PLOT PPHYS
    lWidth = 2;
    fSize = 12;
    if logConvert
        baseVal = reallog(.1);
        yMin = reallog(.1); yMax = reallog(5);
        logOpts = {'ytick',reallog([.1,.25,.5,1,2,5]),'yticklabel',[.1,.25,.5,1,2,5]};
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
            h2D = bar(z,stimGrand(z),'facecolor',color1(expNum,:),'edgecolor','none','basevalue',baseVal);
        else
            h3D = bar(z,stimGrand(z),'facecolor',color2(expNum,:),'edgecolor','none','basevalue',baseVal);
        end
    end
    errorb(1:8,stimGrand,stimStderr);
    set(gca,gcaOpts{:},'xtick',[2.5,6.5],'xticklabel',{'horizontal','vertical'},logOpts{:});
    
    %legend([h2D,h3D],pphysLabels,'location','northwest','fontsize',fSize,'fontname','Arial');
    %legend boxoff;
    plot(1:2:length(conditions),ones(1,length(conditions)/2).*-1.5,'k^','LineWidth',lWidth,'markerfacecolor',[1 1 1],'markeredgecolor','none','MarkerSize',10);
    plot(2:2:length(conditions),ones(1,length(conditions)/2).*-1.5,'kv','LineWidth',lWidth,'markerfacecolor',[1 1 1],'markeredgecolor','none','MarkerSize',10);

    hold off
    xlim([.5,8.5]);
    ylim([yMin,yMax]);
    ylabel('displacement (arcmins)','fontsize',fSize,'fontname','Arial')
    subplot(1,2,2);
    hold on
    for z = 1:4
        if isempty(strfind(bothLabels{z},'3D'));
            h2D = bar(z,bothGrand(z),'facecolor',color1(expNum,:),'edgecolor','none','basevalue',baseVal);
        else
            h3D = bar(z,bothGrand(z),'facecolor',color2(expNum,:),'edgecolor','none','basevalue',baseVal);
        end
    end
    errorb(1:4,bothGrand,bothStderr,'barwidth',.5,'linewidth',lWidth);
    xlim([.5,4.5]);
    ylim([yMin,yMax]);
    set(gca,gcaOpts{:},'xtick',[1.5,3.5],'xticklabel',{'horizontal','vertical'},logOpts{:});
    %legend([h2D,h3D],pphysLabels,'location','northwest','fontsize',fSize,'fontname','Arial');
    %legend boxoff;
    hold off
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = 6;
    figPos(3) = 18;
    set(gcf,'pos',figPos);
    ylabel('displacement (arcmins)','fontsize',fSize,'fontname','Arial');

    export_fig(sprintf('%s/pphys%0.0f_beh.pdf',savePath,expNum),'-pdf','-transparent',pphysFig);
    savefig(gcf,sprintf('%s/pphys%0.0f_beh.fig',savePath,expNum)); 
    %close(pphysFig);

    %% LOAD IN EEG DATA
    rcNum = 1;
    
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
        
        for c = 1:(length(conditions)/2)
            if flipEEG
                % FLIP CONDITIONS 5-8
                flipIdx = zeros(size(subBinIdx{s}{c}));
                freqCount = 0;
                for q = 1:length(subBinIdx{s}{c})
                    flipIdx(q) = (10+1-subBinIdx{s}{c}(q))+(freqCount*10);
                    if subBinIdx{s}{c}(q) == 10
                        freqCount = freqCount+1;
                    else
                    end
                end
                flipIdx = [flipIdx;flipIdx+length(flipIdx)]; % add imag indices
                sensorData{c+c,s} = sensorData{c+c,s}(flipIdx,:,:);
                cellNoiseData1{c+c,s} = cellNoiseData1{c+c,s}(flipIdx,:,:);
                cellNoiseData2{c+c,s} = cellNoiseData2{c+c,s}(flipIdx,:,:);
            else
            end
            if mergeEEG
                % ADD THEM TO 1-4
                sensorData{c+c-1,s} = cat(3,sensorData{c+c-1,s},sensorData{c+c,s});
                cellNoiseData1{c+c-1,s} = cat(3,cellNoiseData1{c+c-1,s},cellNoiseData1{c+c,s});
                cellNoiseData2{c+c-1,s} = cat(3,cellNoiseData2{c+c-1,s},cellNoiseData2{c+c,s});
            else
            end
        end
    end
    if mergeEEG
        sensorData = sensorData(1:2:length(conditions),:);
        cellNoiseData1 = cellNoiseData1(1:2:length(conditions),:);
        cellNoiseData2 = cellNoiseData2(1:2:length(conditions),:);
        conditions = conditions(1:(length(conditions)/2));
    else
    end
    
    %% CHECK FREQUENCY AND BIN INDICES FOR CONSISTENCY ACROSS SUBS
    nSubs = size(sensorData,2);
    for s=1:nSubs
        for c = 1:length(conditions)
            if sum(abs(subFreqIdx{s}{c}-subFreqIdx{1}{c}))~=0 && sum(abs(subBinIdx{s}{c}-subBinIdx{1}{c}))~=0
                error('Frequency and bin indices vary across subjects: check consistency of DFT/RLS exports\n.');
            end
        end
    end

    %% PROJECT DATA ONTO RCs
    % do this separately for each harmonic
    load(rcaPath,'readyRCA');
    warning('off','all')
    doNR = false(4,8,8); % 4 freqs, 8 RCs (with comparison), 8 conditions
    doNR([1,2,4],1,:) = true; % do fitting for first RC, first, second and fourth harmonic, all conditions
    keepConditions = true;
    errorType = 'SEM';
    trialError = false;
    for f = 1:4
        if strcmp(rcaType,'all'); 
            useFreq = f + 8;
        else
            useFreq = f;
        end

        % repopulate RCA struct with pphys data
        rcaData = rcaProject(sensorData,readyRCA(useFreq).W);
        noiseData.lowerSideBand=rcaProject(cellNoiseData1,readyRCA(useFreq).W); 
        noiseData.higherSideBand=rcaProject(cellNoiseData2,readyRCA(useFreq).W);

        % create a "component" of just one channel for performance evaluation if requested
        chanToCompare = 75;
        nChannels = size(sensorData{1},2);
        wComparison=zeros(nChannels,1); wComparison(chanToCompare)=1; 
        comparisonData=rcaProject(sensorData,wComparison); 
        comparisonNoiseData.lowerSideBand =rcaProject(cellNoiseData1,wComparison); 
        comparisonNoiseData.higherSideBand =rcaProject(cellNoiseData2,wComparison); 

        pphysRCA(f) = readyRCA(useFreq);
        fIdx = repmat(subFreqIdx{3}{1}==f,2,1); % repmat because the first half is real, second half is imag with same ordering
        pphysRCA(f).data = cellfun(@(x) x(fIdx,:,:),rcaData,'uni',false);
        pphysRCA(f).noiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.lowerSideBand,'uni',false);
        pphysRCA(f).noiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),noiseData.higherSideBand,'uni',false);
        pphysRCA(f).comparisonData = cellfun(@(x) x(fIdx,:,:),comparisonData,'uni',false);
        pphysRCA(f).comparisonNoiseData.lowerSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.lowerSideBand,'uni',false);
        pphysRCA(f).comparisonNoiseData.higherSideBand = cellfun(@(x) x(fIdx,:,:),comparisonNoiseData.higherSideBand,'uni',false);
        % put into output struct
        pphysRCA(f).settings.rcaConds = 1:length(conditions);
        pphysRCA(f).settings.binLabels = ( arrayfun(@(x) num2str(x,'%.04f'),sweepValues,'uni',false) )';
        pphysRCA(f).settings = orderfields(pphysRCA(f).settings);
    
        rcStruct = aggregateData(pphysRCA(f),keepConditions,errorType,trialError,doNR);        
        % RC
        pphysRCA(f).stats.Amp = squeeze(rcStruct.ampBins);
        pphysRCA(f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp);
        pphysRCA(f).stats.ErrLB = squeeze(rcStruct.ampErrBins(:,:,:,:,1));
        pphysRCA(f).stats.ErrUB = squeeze(rcStruct.ampErrBins(:,:,:,:,2));
        pphysRCA(f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins);
        pphysRCA(f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise);
        % Naka-Rushton
        pphysRCA(f).stats.NR_Params = squeeze(rcStruct.NakaRushton.Params);
        pphysRCA(f).stats.NR_R2 = squeeze(rcStruct.NakaRushton.R2);
        pphysRCA(f).stats.NR_JKSE = squeeze(rcStruct.NakaRushton.JackKnife.SE);
        pphysRCA(f).stats.NR_JKParams = squeeze(rcStruct.NakaRushton.JackKnife.Params);
        pphysRCA(f).stats.hModel = rcStruct.NakaRushton.hModel;
        pphysRCA(f).stats.tSqrdP = squeeze(rcStruct.tSqrdP);
        pphysRCA(f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig);
        pphysRCA(f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal);
    end
    warning('on','all')
    
    
    %% PLOT RCs
    close all;
    if mergeEEG
         pphysConds{1} = 1:2;
         pphysConds{2} = 3:4;
         if expNum == 1
            pphysColors = [1,3,1,3];
         else
            pphysColors = [2,4,2,4];
         end
    else
        pphysConds{1} = 1:4;
        pphysConds{2} = 5:8; % you are plotting 2D rel and 3D rel x 2
        if expNum == 1
            pphysColors = [1,1,3,3,1,1,3,3];
        else
            pphysColors = [2,2,4,4,2,2,4,4];
         end
    end
    % plot settings
    % set figure size in the beginning
    figHeight = 30;
    figWidth = 20;

    binVals = cellfun(@(x) str2num(x), pphysRCA(1).settings.binLabels);

    lWidth = 1.5;
    %red =   [152 53 49; 231 184 182]./255;
    %green = [117 148 54; 216 229 186]./255;
    %blue =  [52 95 148; 183 204 229]./255;
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    subColors = repmat([color1; color2],2,1);
    subColors = subColors(pphysColors,:);
    fSize = 12;
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
    
    rcNum = 1; freqNum = 2; % plot first RC, second harmonic
    
    figure;
    
    for s = 1:2
        curConds = find(ismember(pphysRCA(freqNum).settings.rcaConds,pphysConds{s}));
        % compute new signal values, averaged over bins
        [rcaDataReal,rcaDataImag] = getRealImag(pphysRCA(freqNum).data);
        rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataReal,'uni',false);
        rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
        rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataImag,'uni',false);
        rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
        realBinMean = squeeze(nanmean(rcaDataReal));
        imagBinMean = squeeze(nanmean(rcaDataImag));
        % compute new noise values, averaged over bins
        [noiseLoReal,noiseLoImag] = getRealImag(pphysRCA(freqNum).noiseData.lowerSideBand);
        [noiseHiReal,noiseHiImag] = getRealImag(pphysRCA(freqNum).noiseData.lowerSideBand);
        noiseReal = cellfun(@(x,y) (x+y)./2, noiseLoReal,noiseHiReal, 'uni',false);
        noiseImag = cellfun(@(x,y) (x+y)./2, noiseLoImag,noiseHiImag, 'uni',false);
        noiseReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),noiseReal,'uni',false);
        noiseReal = cell2mat(permute(noiseReal,[3,2,1]));
        noiseImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),noiseImag,'uni',false);
        noiseImag = cell2mat(permute(noiseImag,[3,2,1]));
        realBinMeanNoise = squeeze(nanmean(noiseReal));
        imagBinMeanNoise = squeeze(nanmean(noiseImag));
        
        if ~projectedData
            % grab values from data structure
            valSet = squeeze(pphysRCA(freqNum).stats.Amp(:,rcNum,:));
            errSet1 = squeeze(pphysRCA(freqNum).stats.ErrLB(:,rcNum,:));
            errSet2 = squeeze(pphysRCA(freqNum).stats.ErrLB(:,rcNum,:));
            % compute vector means of mean bins, and errors
            valSet(length(binVals)+1,:) = sqrt(nanmean(realBinMean,1).^2+nanmean(imagBinMean,1).^2);

            % store p-values for later

            % compute elliptical error and do Hotelling's T2 against zero
            for c=1:length(curConds)
                statIdx = c+(s-1)*(length(curConds)+1);
                rc_tSqrdP(1:length(binVals),statIdx,s) = pphysRCA(freqNum).stats.tSqrdP(:,rcNum,curConds(c));
                rc_tSqrdVal(1:length(binVals),statIdx) = pphysRCA(freqNum).stats.tSqrdVal(:,rcNum,curConds(c));
                xyBinMean = cat(2,realBinMean(:,curConds(c)),imagBinMean(:,curConds(c)));
                nanVals = sum(isnan(xyBinMean),2)>0;      
                [binErrs,~,~,errorEllipse] = fitErrorEllipse(xyBinMean(~nanVals,:),'SEM');
                errSet1(length(binVals)+1,c) = binErrs(1);
                errSet2(length(binVals)+1,c) = binErrs(2);
                % compute t-values
                tStruct = tSquaredFourierCoefs(xyBinMean(~nanVals,:));
                rc_tSqrdP(length(binVals)+1,statIdx) = tStruct.pVal;
                rc_tSqrdVal(length(binVals)+1,statIdx) = tStruct.tSqrd;
                % note, just using mean df for all values, would need
                % to be fixed if multi data was ever used seriously
                rc_tSqrdDF(1:length(binVals)+1,statIdx) = tStruct.df2;
            end
        else
            statIdx = (1:length(curConds))+(s-1)*(length(curConds)+1);
            % compute projected vector mean
            % move subjects to first dim
            realVector = permute(rcaDataReal(:,:,curConds),[2,1,3]);
            imagVector = permute(rcaDataImag(:,:,curConds),[2,1,3]);
            project_amps = vectorProjection(realVector,imagVector);
            % compute mean of projected vector amplitude, over bins
            project_amps(:,end+1,:) = nanmean(project_amps,2);
            % if all values are NaN
            if any(sum(squeeze(all(isnan(project_amps),2)),2)>0)
                nan_subs = find(sum(squeeze(all(isnan(project_amps),2)),2)>0);
                for z = 1:length(nan_subs)
                    msg = ...
                        sprintf('Subject %d has no values in one or more conditions, setting all values to NaN', ...
                        nan_subs(z));
                    warning(msg);
                    project_amps(nan_subs(z),:,:) = NaN;
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
            rc_tSqrdP(:,statIdx) = permute(temp_p,[2,3,1]);
            rc_tSqrdVal(:,statIdx) = permute(temp_stats.tstat,[2,3,1]);
            rc_tSqrdDF(:,statIdx) = permute(temp_stats.df,[2,3,1]);
            clear temp_*;
        end

        % maximum noise across all four conditions being plotted
        % since this is just means, we can compute it the same way for
        % projected and not-projected
        noiseSet = max(pphysRCA(freqNum).stats.NoiseAmp(:,rcNum,:),[],3);
        % max of mean over bins
        noiseSet(length(binVals)+1) = max(sqrt(nanmean(realBinMeanNoise,1).^2+nanmean(imagBinMeanNoise,1).^2));
        
        % note, the paired tests below do not really works unless data are merged
        statIdx = length(curConds)+1+(s-1)*(length(curConds)+1);
        for b = 1:(length(binVals)+1)
            if ~projectedData
                if b <= (length(binVals))
                    xyData = permute(cat(1,rcaDataReal(b,:,curConds),rcaDataImag(b,:,curConds)),[2,1,3]);
                else
                    xyData = cat(3,realBinMean(:,curConds), imagBinMean(:,curConds));
                end
                tempStrct = tSquaredFourierCoefs(xyData);
                rc_tSqrdP(b,statIdx) = tempStrct.pVal;
                rc_tSqrdSig(b,statIdx) = tempStrct.pVal < 0.05;
                rc_tSqrdVal(b,statIdx) = tempStrct.tSqrd;
                rc_tSqrdDF(b,statIdx) = tempStrct.df2;
                rc_tSqrdD(b,s) = tempStrct.mahalanobisD;
                rc_tSqrdU(b,s) = tempStrct.cohenNonOverlap;
                rc_tSqrdMu1(b,s) = valSet(b,1);
                rc_tSqrdMu2(b,s) = valSet(b,2);
            else
                curData = squeeze(project_amps(:,b,:));
                not_nan = ~any(isnan(curData),2);
                [rc_tSqrdSig(b,statIdx),rc_tSqrdP(b,statIdx),~,tempStrct] = ttest(curData(not_nan,1),curData(not_nan,2),'alpha',0.05,'dim',1,'tail','both');
                rc_tSqrdVal(b,statIdx) = tempStrct.tstat;
                rc_tSqrdDF(b,statIdx) = tempStrct.df;
                rc_tSqrdD(b,s) = rc_tSqrdVal(b,statIdx)./sqrt(tempStrct.df+1); % Cohen's D
                OVL = 2*normcdf(-abs(rc_tSqrdD(b,s)/2));
                rc_tSqrdU(b,s) = 1-OVL/(2-OVL); 
                rc_tSqrdMu1(b,s) = mean(curData(not_nan,1));
                rc_tSqrdMu2(b,s) = mean(curData(not_nan,2));
            end
        end
                            
        % make Naka-Rushton values
        NRset = pphysRCA(freqNum).stats.NR_Params;          
        NRerrs = pphysRCA(freqNum).stats.NR_JKSE;
        NRmodel = pphysRCA(freqNum).stats.hModel;

        % do paired tests
        testVal = squeeze(pphysRCA(freqNum).stats.NR_JKParams(:,:,rcNum,curConds));
        testVal = permute(testVal,[2,1,3]); % move subjects to first dim    
        paramIdx = [1,2,3,4];% only look at c50 and rMax
        jkDf = size(testVal,1)-1;
        diffErr = jackKnifeErr(testVal(:,paramIdx,1)-testVal(:,paramIdx,2));
        grandDiff = NRset(paramIdx,rcNum,curConds(1)) - NRset(paramIdx,rcNum,curConds(2));
        paramPairedT(:,s) = grandDiff'./diffErr;
        paramPairedP(:,s) = 2*tcdf( -abs(paramPairedT(:,s)) , jkDf);
        
        % DO PLOTTING
        binVals = cellfun(@(x) str2num(x), pphysRCA(freqNum).settings.binLabels);
        logStep = diff(reallog(binVals(1:2))); % step size
        xMin = exp(reallog(binVals(1))-logStep*.5);
        xMax = exp(reallog(binVals(end))+logStep*2.5); % add 2.5 steps
        extraBins = arrayfun(@(x) exp(reallog(binVals(end))+x), [logStep,logStep*2]);
        
        spH = subplot(1,2,s);
        
        hold on
        for c = 1:length(curConds)
            if mergeEEG
                ampMarkerStyle = {'o','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',5}; % 'MarkerEdgeColor','none'
            else
                if isempty(strfind(condLabels{curConds(c)},'A'))
                    % descending, download triangle
                    ampMarkerStyle = {'^','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7}; % 'MarkerEdgeColor','none'
                else
                    % ascending, upward triangle
                    ampMarkerStyle = {'v','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7};
                end
            end
            valH(c)=plot([binVals;extraBins((curConds(c)==max(curConds))+1)],valSet(:,c),ampMarkerStyle{:});
            hE = ErrorBars([binVals;extraBins((curConds(c)==max(curConds))+1)],valSet(:,c),[errSet1(:,c),errSet2(:,c)],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
            uistack(valH(c),'bottom')
            cellfun(@(x) uistack(x,'bottom'), hE);
            hold on

            % plot Naka-Rushton
            nFine = 1e2;
            nrX = linspace( min(binVals), max(binVals), nFine )';
            nrVals = NRmodel( nrX, NRset(:,rcNum,curConds(c)));
            hNR{c} = plot( nrX, nrVals, '-','color',subColors(curConds(c),:), 'LineWidth',lWidth);
        end
        cellfun(@(x) uistack(x,'bottom'), hNR);   
        yUnit = 1;
        yMin = 0;
        yMax = 5;
        legend(valH,{'in-phase','anti-phase'},'fontsize',fSize,'fontname','Helvetica','location','northwest');
        legend boxoff
        if s == 1
            title('horizontal','fontsize',fSize,'fontname','Helvetica')
        else
            title('vertical','fontsize',fSize,'fontname','Helvetica')
        end
        
        set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');

        xlim([xMin,xMax]);
        ylim([yMin,yMax])
%         text(xMin-logStep*.75,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
%             plotLabel(f),'fontsize',fSize*2,'fontname','Helvetica');
%         text(xMax+logStep*1.5,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
%             plotLabel(f+2),'fontsize',fSize*2,'fontname','Helvetica');

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
        
        if projectedData && mergeEEG
            curFig = gcf;
            tmpFig = figure;
            tmpAx = get(spH,'children');
            copyobj(allchild(spH),gca(tmpFig));
            set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');
            xlim([xMin,xMax]);
            ylim([0,yMax])
            if s == 2
                savefig(sprintf('%s/subplot%0.0f_vert.fig',savePath,expNum));
            else
                savefig(sprintf('%s/subplot%0.0f_hori.fig',savePath,expNum));
            end
            close(tmpFig);
            figure(curFig);
        else
        end
    end
    
    if ~projectedData && ~mergeEEG
        return
    else
    end
    
    %% MAKE TABLE
    if ~all(rc_tSqrdDF(:) == rc_tSqrdDF(1,1))
        error('different dfs for different tests');
    else
    end

    stringPairedP = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdP,'uni',false);
    sigIdx = cell2mat(arrayfun(@(x) x < 0.0001,rc_tSqrdP,'uni',false));
    stringPairedP(sigIdx) = {'<0.0001'};
    stringPairedSig = cell(size(rc_tSqrdSig));
    stringPairedSig(rc_tSqrdSig==1) = {'*'};
    stringPairedSig(rc_tSqrdSig==0) = {'ns'};
    stringPairedT = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdVal,'uni',false);
    stringBinVals = cat(1,arrayfun(@(x) num2str(x,'%0.2f'),binVals,'uni',false),{'n/a'});
    stringPairedD = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdD,'uni',false);
    stringPairedDegFree = arrayfun(@(x) num2str(x,'%0.0f'),min(rc_tSqrdDF'),'uni',false)';
    stringPairedU = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdU,'uni',false);
    stringPairedMu1 = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdMu1,'uni',false);
    stringPairedMu2 = arrayfun(@(x) num2str(x,'%0.4f'),rc_tSqrdMu2,'uni',false);
    finishedArray = [];
    
    for s = 1:2
        if s == 1
            varLabels = cat(2,{num2str(rc_tSqrdDF(1),'horizontal (df=%0.0f)')},...
                  arrayfun(@(x) num2str(x,'bin%0.0f'),1:length(binVals),'uni',false),...
                  {'ave'});
            readyArray = [stringBinVals, stringPairedP(:,s*3),stringPairedT(:,s*3),stringPairedD(:,s)]';
            rowLabels = {'disp (arcmins)','p','t-statistic','Cohen''s D'}';
        else
            varLabels = cat(2,{num2str(rc_tSqrdDF(1),'vertical (df=%0.0f)')},...
                  arrayfun(@(x) num2str(x,'bin%0.0f'),1:length(binVals),'uni',false),...
                  {'ave'});
            readyArray = [stringBinVals, stringPairedP(:,s*3),stringPairedT(:,s*3),stringPairedD(:,s)]';
            rowLabels = {'disp (arcmins)','p','t-statistic','Cohen''s D'}';
        end
        readyArray = cat(2,rowLabels,readyArray);
        readyArray = cat(1,varLabels,readyArray);
        finishedArray = cat(1,finishedArray,readyArray);
    end
    if expNum == 1
        titleStr = 'EEG recorded during psychophysics: referenced conditions';
    else
        titleStr = 'EEG recorded during psychophysics: unreferenced conditions';
    end
    
    titleLabels = cell(size(varLabels));
    titleLabels{1} = titleStr;

    finishedArray = cat(1,titleLabels,finishedArray);
    table_file = sprintf('%s/pphys_table%0.0f.csv',savePath,expNum);
    pphysTable = array2table(finishedArray);
    writetable(pphysTable,table_file,'WriteRowNames',false,'WriteVariableNames',false);

    if exist(sprintf('%s/pphys_table1.csv',savePath)) && exist(sprintf('%s/pphys_table2.csv',savePath))
        combined_file = sprintf('%s/pphys_combined.csv',savePath);
        if exist(combined_file,'file')
            delete(combined_file);
        else
        end
        catCmd{1} = sprintf('cd %s',savePath);
        catCmd{2} = 'cat pphys_table1.csv pphys_table2.csv > pphys_combined.csv';
        system(strjoin(catCmd, '; '));
    else
    end
    
    %% PLOT COMBINED FIGURE
    close all;
    % experiment 1, horizontal
    rcaH(1) = openfig(sprintf('%s/subplot%0.0f_hori.fig',savePath,1)); 
    rcaAx(1) = gca;
    % experiment 1, vertical
    rcaH(2) = openfig(sprintf('%s/subplot%0.0f_vert.fig',savePath,1)); 
    rcaAx(2) = gca;
    % experiment 2, horizontal
    rcaH(3) = openfig(sprintf('%s/subplot%0.0f_hori.fig',savePath,2)); 
    rcaAx(3) = gca;
    % experiment 2, vertical
    rcaH(4) = openfig(sprintf('%s/subplot%0.0f_vert.fig',savePath,2)); 
    rcaAx(4) = gca;

    behH(1) = openfig(sprintf('%s/pphys%0.0f_beh.fig',savePath,1)); 
    behAx(1) = gca;
    behH(2) = openfig(sprintf('%s/pphys%0.0f_beh.fig',savePath,2)); 
    behAx(2) = gca;

    figure; %create new figure
    for z = 1:4
        spH = subplot(2,3,z+(z>2)); %create and get handle to the subplot axes
        figH = get(rcaAx(z),'children'); %get handle to all the children in the figure
        copyobj(figH,spH); %copy children to new parent axes i.e. the subplot axes
        delete(findall(spH,'type','text'));
        xMin = min(get(rcaAx(z),'xlim'));
        xMax = max(get(rcaAx(z),'xlim'));
        if z < 3
            xTick = [0,0.2,0.5,1,2];
        else
            xTick = [0,0.2,0.5,1,2,4,8,16];
        end
        yUnit = 1;
        yMax = 3;
        gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
        set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',xTick,'ytick',0:yUnit:yMax,'Layer','top');
        xlim([xMin,xMax]);
        ylim([0,yMax])
        
        newAx = axes('position',get(gca,'position'));
        set(newAx,'visible','off');
        xText(1) = min(get(newAx,'xlim'))-diff(get(newAx,'xlim'))*.15;
        xText(2) = min(get(newAx,'xlim'))+diff(get(newAx,'xlim'))*.15;
        yText(1) = max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1;
        yText(2) = max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.1;
        yText(3) = max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.2;
        
        switch z
            case 1
                plotLabel = 'A';
                title('horizontal','fontsize',fSize,'fontname','Helvetica')
                text(xText(2),yText(2),'reference','fontsize',fSize,'fontname','Helvetica')
            case 2
                plotLabel = 'B';
                title('vertical','fontsize',fSize,'fontname','Helvetica')
                text(xText(2),yText(2),'reference','fontsize',fSize,'fontname','Helvetica')
            case 3
                plotLabel = 'D';
                title('horizontal','fontsize',fSize,'fontname','Helvetica')
                ylabel('amplitude (\muV)','fontsize',fSize,'fontname','Helvetica')
                xlabel('displacement (arcmins)','fontsize',fSize,'fontname','Helvetica');
                text(xText(2),yText(2),'no reference','fontsize',fSize,'fontname','Helvetica')
            case 4
                title('vertical','fontsize',fSize,'fontname','Helvetica')
                plotLabel = 'E';
                text(xText(2),yText(2),'no reference','fontsize',fSize,'fontname','Helvetica')
            otherwise
        end
        
        text(xText(1),yText(1),plotLabel,'fontsize',fSize*2,'fontname','Helvetica');
        warning('figure assumes same n for both experiments, make sure this is the case');
        text(xText(2),yText(3),sprintf('n = %0.0f',rc_tSqrdDF(1,1)+1),'fontsize',fSize,'fontname','Helvetica')
        if z == 4
            for q = 1:2
                spH = subplot(2,3,3+(q-1)*3); %create and get handle to the subplot axes
                figH = get(behAx(q),'children'); %get handle to all the children in the figure
                copyobj(figH,spH); %copy children to new parent axes i.e. the subplot axes
                delete(findall(spH,'type','text'));
                if logConvert
                    baseVal = reallog(.1);
                    yMin = reallog(.1); yMax = reallog(5);
                    logOpts = {'ytick',reallog([.1,.25,.5,1,2,5]),'yticklabel',[.1,.25,.5,1,2,5]};
                else
                    baseVal = 0;
                    yMin = 0; yMax = 8;
                    logOpts ={'ytick',yMin:1:yMax};
                end
                ylabel('threshold (arcmins)','fontsize',fSize,'fontname','Helvetica')
                xlim([.5,4.5]);
                ylim([yMin,yMax]);
                set(gca,gcaOpts{:},'xtick',[1.5,3.5],'xticklabel',{'horizontal','vertical'},logOpts{:});
                 if q == 1
                    title('reference','fontsize',fSize,'fontname','Helvetica');
                    newAx = axes('position',get(gca,'position'));
                    text(xText(1),yText(1),'C','fontsize',fSize*2,'fontname','Helvetica');
                    set(newAx,'visible','off');
                else
                    title('no reference','fontsize',fSize,'fontname','Helvetica');
                    newAx = axes('position',get(gca,'position'));
                    text(xText(1),yText(1),'F','fontsize',fSize*2,'fontname','Helvetica');
                    set(newAx,'visible','off');
                 end
            end
        else
        end  
    end
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = 12;
    figPos(3) = 24;
    set(gcf,'pos',figPos);
    export_fig(sprintf('%s/figure4_complete.pdf',savePath),'-pdf','-transparent',gcf);
    
    close all
    
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
    
    return;

    %% STUFF BELOW IS OLD
    
    close all;
    if mergeEEG
         pphysConds{1} = 1:2;
         pphysConds{2} = 3:4;
         if expNum == 1
            pphysColors = [1,3,1,3];
         else
            pphysColors = [2,4,2,4];
         end
    else
        pphysConds{1} = 1:4;
        pphysConds{2} = 5:8; % you are plotting 2D rel and 3D rel x 2
        if expNum == 1
            pphysColors = [1,1,3,3,1,1,3,3];
        else
            pphysColors = [2,2,4,4,2,2,4,4];
         end
    end
    % plot settings
    % set figure size in the beginning
    figHeight = 30;
    figWidth = 20;

    binVals = cellfun(@(x) str2num(x), pphysRCA(1).settings.binLabels);

    lWidth = 1.5;
    %red =   [152 53 49; 231 184 182]./255;
    %green = [117 148 54; 216 229 186]./255;
    %blue =  [52 95 148; 183 204 229]./255;
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    subColors = repmat([color1; color2],2,1);
    subColors = subColors(pphysColors,:);
    fSize = 12;
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
    for f=1:length(freqsToUse)
        freqFig(f) = figure;
        curFreq = freqsToUse(f);
        
        % make Naka-Rushton values
        NRmodel = pphysRCA(curFreq).stats.hModel;
        NRset = squeeze(pphysRCA(curFreq).stats.NR_Params(:,rcNum,:));          
        NRerrs = pphysRCA(curFreq).stats.NR_JKSE;
        NRmodel = pphysRCA(curFreq).stats.hModel;

        % compute new signal values, averaged over bins
        [rcaDataReal,rcaDataImag] = getRealImag(pphysRCA(curFreq).data);
        rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataReal,'uni',false);
        rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
        rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataImag,'uni',false);
        rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
        realBinMean = squeeze(nanmean(rcaDataReal));
        imagBinMean = squeeze(nanmean(rcaDataImag));
        % compute new noise values, averaged over bins
        [noiseLoReal,noiseLoImag] = getRealImag(pphysRCA(curFreq).noiseData.lowerSideBand);
        [noiseHiReal,noiseHiImag] = getRealImag(pphysRCA(curFreq).noiseData.lowerSideBand);
        noiseReal = cellfun(@(x,y) (x+y)./2, noiseLoReal,noiseHiReal, 'uni',false);
        noiseImag = cellfun(@(x,y) (x+y)./2, noiseLoImag,noiseHiImag, 'uni',false);
        noiseReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),noiseReal,'uni',false);
        noiseReal = cell2mat(permute(noiseReal,[3,2,1]));
        noiseImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),noiseImag,'uni',false);
        noiseImag = cell2mat(permute(noiseImag,[3,2,1]));
        realBinMeanNoise = squeeze(nanmean(noiseReal));
        imagBinMeanNoise = squeeze(nanmean(noiseImag));
        
        if ~projectedData
            % grab values from data structure
            valSet = squeeze(pphysRCA(curFreq).stats.Amp(:,rcNum,:));
            errSet1 = squeeze(pphysRCA(curFreq).stats.ErrLB(:,rcNum,:));
            errSet2 = squeeze(pphysRCA(curFreq).stats.ErrLB(:,rcNum,:));
            % compute vector means of mean bins, and errors
            valSet(length(binVals)+1,:) = sqrt(nanmean(realBinMean,1).^2+nanmean(imagBinMean,1).^2);

            % store p-values for later

            % compute elliptical error and do Hotelling's T2 against zero
            for c=1:length(condsToUse)
                rc_tSqrdP(1:length(binVals),c+(f-1)*3) = pphysRCA(curFreq).stats.tSqrdP(:,rcNum,c);
                rc_tSqrdVal(1:length(binVals),c+(f-1)*3) = pphysRCA(curFreq).stats.tSqrdVal(:,rcNum,c);
                xyBinMean = cat(2,realBinMean(:,c),imagBinMean(:,c));
                nanVals = sum(isnan(xyBinMean),2)>0;      
                [binErrs,~,~,errorEllipse] = fitErrorEllipse(xyBinMean(~nanVals,:),'SEM');
                errSet1(length(binVals)+1,c) = binErrs(1);
                errSet2(length(binVals)+1,c) = binErrs(2);
                % compute t-values
                tStruct = tSquaredFourierCoefs(xyBinMean(~nanVals,:));
                rc_tSqrdP(length(binVals)+1,c+(f-1)*3) = tStruct.pVal;
                rc_tSqrdVal(length(binVals)+1,c+(f-1)*3) = tStruct.tSqrd;
                % note, just using mean df for all values, would need
                % to be fixed if multi data was ever used seriously
                rc_tSqrdDF(1:length(binVals)+1,c+(f-1)*3) = tStruct.df2;
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
                        sprintf('Subject %d has no values in one or more conditions, setting all values to NaN', ...
                        nan_subs(z));
                    warning(msg);
                    project_amps(nan_subs(z),:,:) = NaN;
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
            rc_tSqrdP(:,(1:2)+(f-1)*3) = permute(temp_p,[2,3,1]);
            rc_tSqrdVal(:,(1:2)+(f-1)*3) = permute(temp_stats.tstat,[2,3,1]);
            rc_tSqrdDF(:,(1:2)+(f-1)*3) = permute(temp_stats.df,[2,3,1]);
            clear temp_*;
        end

        % maximum noise across all four conditions being plotted
        % since this is just means, we can compute it the same way for
        % projected and not-projected
        noiseSet = max(pphysRCA(curFreq).stats.NoiseAmp(:,rcNum,:),[],3);
        % max of mean over bins
        noiseSet(length(binVals)+1) = max(sqrt(nanmean(realBinMeanNoise,1).^2+nanmean(imagBinMeanNoise,1).^2));

        for b = 1:(length(binVals)+1)
            if ~projectedData
                if b <= (length(binVals))
                    xyData = permute(cat(1,rcaDataReal(b,:,:),rcaDataImag(b,:,:)),[2,1,3]);
                else
                    xyData = cat(3,realBinMean, imagBinMean);
                end
                tempStrct = tSquaredFourierCoefs(xyData);
                rc_tSqrdP(b,3+(f-1)*3) = tempStrct.pVal;
                rc_tSqrdSig(b,3+(f-1)*3) = tempStrct.pVal < 0.05;
                rc_tSqrdVal(b,3+(f-1)*3) = tempStrct.tSqrd;
                rc_tSqrdDF(b,3+(f-1)*3) = tempStrct.df2;
                rc_tSqrdD(b,f) = tempStrct.mahalanobisD;
                rc_tSqrdU(b,f) = tempStrct.cohenNonOverlap;
                rc_tSqrdMu1(b,f) = valSet(b,1);
                rc_tSqrdMu2(b,f) = valSet(b,2);
            else
                curData = squeeze(project_amps(:,b,:));
                not_nan = ~any(isnan(curData),2);
                [rc_tSqrdSig(b,3+(f-1)*3),rc_tSqrdP(b,3+(f-1)*3),~,tempStrct] = ttest(curData(not_nan,1),curData(not_nan,2),'alpha',0.05,'dim',1,'tail','both');
                rc_tSqrdVal(b,3+(f-1)*3) = tempStrct.tstat;
                rc_tSqrdDF(b,3+(f-1)*3) = tempStrct.df;
                rc_tSqrdD(b,f) = rc_tSqrdVal(b,3+(f-1)*3)./sqrt(tempStrct.df+1); % Cohen's D
                OVL = 2*normcdf(-abs(rc_tSqrdD(b,f)/2));
                rc_tSqrdU(b,f) = 1-OVL/(2-OVL); 
                rc_tSqrdMu1(b,f) = mean(curData(not_nan,1));
                rc_tSqrdMu2(b,f) = mean(curData(not_nan,2));
            end
        end
        
        
        for r = 1:(nComp+1)
            if r<6
                egiH(r) = subplot(nComp+1,3,3+(r-1)*3);
                hold on
                rcaColorBar = [min(pphysRCA(curFreq).A(:,r)),max(pphysRCA(curFreq).A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(pphysRCA(curFreq).A(:,r)*-1,rcaColorBar);
                hold off
            else
            end
            for s=1:2 % vertical or horizontal motion
                if s==1
                    spH = subplot(nComp+1,3,1+(r-1)*3);
                    curConds = find(ismember(pphysRCA(curFreq).settings.rcaConds,pphysConds{s}));
                    titleStr = sprintf('horizontal: %s',pphysRCA(curFreq).settings.freqLabels{1});
                else
                    spH = subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(pphysRCA(curFreq).settings.rcaConds,pphysConds{s}));
                    titleStr = sprintf('vertical: %s',pphysRCA(curFreq).settings.freqLabels{1});
                end
                hold on
                for c=1:length(curConds)
                    if mergeEEG
                        ampMarkerStyle = {'o','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',5}; % 'MarkerEdgeColor','none'
                    else
                        if isempty(strfind(condLabels{curConds(c)},'A'))
                            % descending, download triangle
                            ampMarkerStyle = {'^','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7}; % 'MarkerEdgeColor','none'
                        else
                            % ascending, upward triangle
                            ampMarkerStyle = {'v','LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1 1 1],'MarkerSize',7};
                        end
                    end
                    
                    valH(c)=plot(binVals,valSet(:,r,curConds(c)),ampMarkerStyle{:});
                    hE = ErrorBars(binVals,valSet(:,r,curConds(c)),[errSet1(:,r,curConds(c)),errSet2(:,r,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
                    uistack(valH(c),'bottom')
                    cellfun(@(x) uistack(x,'bottom'), hE);
                    hold on
                    if ~isnan(NRset(1,r,curConds(c)))
                        % plot Naka-Rushton
                        nFine = 1e2;
                        nrX = linspace( min(binVals), max(binVals), nFine )';
                        nrVals = NRmodel( nrX, NRset(:,r,curConds(c)));
                        hNR{c} = plot( nrX, nrVals, '-','color',subColors(curConds(c),:), 'LineWidth',lWidth);
                    else
                    end
                end
                cellfun(@(x) uistack(x,'bottom'), hNR);
                
                if curFreq == 2 && r == 1     
                    yUnit = 0.5;
                    yMax = 3.5;
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
                end
                   
                %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
                %yMax = ceil(max(valSet(:)))+yUnit;
                set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');
          
                xlim([xMin,xMax]);
                ylim([0,yMax])

                % plot noise patch
                meanNoise = max(pphysRCA(curFreq).stats.NoiseAmp(:,r,curConds),[],3)';
                yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
                xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
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
                        ylabel('amplitude (\muV)')
                        xlabel('displacement (arcmins)');
                    else
                        if mergeEEG
                            lH = legend(valH,pphysLabels,'location','northeast');
                        else
                            lH = legend(valH,rcaLabels,'location','northeast');
                        end
                        legend boxoff
                        lPos = get(lH,'position');
                        lPos(1) = lPos(1) + .23;
                        lPos(2) = lPos(2) + .05;
                        set(lH,'position',lPos);
                        tH = text(3.5,-0.1,sprintf('n = %0d',size(pphysRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Arial');
                    end
                end
                if curFreq == 2 && r == 1
                    
                    % make new p-values
                    [rcaDataReal,rcaDataImag] = getRealImag(pphysRCA(curFreq).data(curConds,:));
                    freqIdx = pphysRCA(curFreq).settings.freqIndices == curFreq;
                    rcaDataReal = cellfun(@(x) squeeze(nanmean(x(freqIdx,rcNum,:),3)),rcaDataReal,'uni',false);
                    rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
                    rcaDataImag = cellfun(@(x) squeeze(nanmean(x(freqIdx,rcNum,:),3)),rcaDataImag,'uni',false);
                    rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
                    
                    rc_tSqrdP(:,1+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdP(:,rcNum,1);
                    rc_tSqrdVal(:,1+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdVal(:,rcNum,1);
                    rc_tSqrdSig(:,1+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdSig(:,rcNum,1);
                    rc_tSqrdP(:,2+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdP(:,rcNum,2);
                    rc_tSqrdVal(:,2+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdVal(:,rcNum,2);
                    rc_tSqrdSig(:,2+(s-1)*3) = pphysRCA(curFreq).stats.tSqrdSig(:,rcNum,2);

                    for b = 1:length(binVals)
                        xyData = permute(cat(1,rcaDataReal(b,:,:),rcaDataImag(b,:,:)),[2,1,3]);
                        tempStrct = tSquaredFourierCoefs(xyData);
                        rc_tSqrdP(b,3+(s-1)*3) = tempStrct.pVal;
                        rc_tSqrdVal(b,3+(s-1)*3) = tempStrct.tSqrd;
                        rc_tSqrdSig(b,3+(s-1)*3) = tempStrct.H;
                        diffMag(b,s) = tempStrct.testAmp;
                    end
                    
                    % make Naka-Rushton values
                    
                    NRvals = pphysRCA(curFreq).stats.NR_Params;
                    NRerrs = pphysRCA(curFreq).stats.NR_JKSE;
                    NRmodel = pphysRCA(curFreq).stats.hModel;
        
                     % do paired tests
                    testVal = squeeze(pphysRCA(curFreq).stats.NR_JKParams(:,:,rcNum,curConds));
                    testVal = permute(testVal,[2,1,3]); % move subjects to first dim    
                    paramIdx = [1,2,3,4];% only look at c50 and rMax
                    jkDf = size(testVal,1)-1;
                    diffErr = jackKnifeErr(testVal(:,paramIdx,1)-testVal(:,paramIdx,2));
                    grandDiff = NRvals(paramIdx,rcNum,curConds(1)) - NRvals(paramIdx,rcNum,curConds(2));
                    paramPairedT(:,s) = grandDiff'./diffErr;
                    paramPairedP(:,s) = 2*tcdf( -abs(paramPairedT(:,s)) , jkDf);
                else
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
        export_fig(sprintf('%s/pphys%0.0f_rc%d_%s.pdf',savePath,expNum,curFreq,rcaType),'-pdf','-transparent',gcf);
    end
    close all;
    
    for z=1:4
        subplot(1,4,z);
        hold on
        plot(1:4,squeeze(NRvals(z,1,:)),'o')
        errorb(1:4,squeeze(NRvals(z,1,:)),squeeze(NRerrs(z,1,:)))
    end



