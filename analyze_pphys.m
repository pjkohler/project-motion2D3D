function analyze_pphys(expNum,taskFig,logConvert,rcaType)
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
        flipEEG = true;
    else
    end
    if nargin < 6
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
    codeFolder = '/Users/kohler/code';
    rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
    setenv('DYLD_LIBRARY_PATH','')

    %% PRELIMINARY
    % folder names
    topPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
    dataPath = sprintf('%s/pphys_exp%0.0f',topPath,expNum);
    subPaths = subfolders(sprintf('%s/*20*',dataPath),1);
    adultExp = [1,2,3,4,5];
    % NOTE! should be path corresponding to blank (2);
    rcaPath = sprintf('%s/figures/exp%0.0f/plottingData.mat',topPath,adultExp(2));
    savePath = sprintf('%s/figures/paper_figures/pphys_exp',topPath);
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
            h2D = bar(z,bothGrand(z),'facecolor',color1(expNum,:),'edgecolor','none','basevalue',baseVal);
        else
            h3D = bar(z,bothGrand(z),'facecolor',color2(expNum,:),'edgecolor','none','basevalue',baseVal);
        end
    end
    errorb(1:4,bothGrand,bothStderr,'barwidth',.5,'linewidth',lWidth);
    xlim([.5,4.5]);
    ylim([yMin,yMax]);
    set(gca,gcaOpts{:},'xtick',[1.5,3.5],'xticklabel',{'Horizontal','Vertical'},logOpts{:});
    %legend([h2D,h3D],pphysLabels,'location','northwest','fontsize',fSize,'fontname','Arial');
    %legend boxoff;
    hold off
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = 6;
    figPos(3) = 18;
    set(gcf,'pos',figPos);
    ylabel('Displacement (arcmins)','fontsize',fSize,'fontname','Arial');

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
    nSubs = size(sensorData,3);
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
    doNR = false(1,6,8); % 1 freqs, 6 RCs (w/ comparison), 8 conditions
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

        % generate final output struct
        pphysRCA(f).data = rcaData;
        pphysRCA(f).noiseData = noiseData;
        pphysRCA(f).comparisonData = comparisonData;
        pphysRCA(f).comparisonNoiseData = comparisonNoiseData;
        pphysRCA(f).inputData = sensorData;
        
        pphysRCA(f).settings.binLevels = subBinLabels{3};
        pphysRCA(f).settings.freqLabels = subFreqLabels{3};
        pphysRCA(f).settings.freqIndices = subFreqIdx{3};
        pphysRCA(f).settings.binIndices = subBinIdx{3};
        pphysRCA(f).settings.freqsToUse = f;
        pphysRCA(f).settings.binsToUse = 1:length(subBinLabels{3}{1});
        pphysRCA(f).settings.condsToUse = 1:length(conditions);
        pphysRCA(f).A = readyRCA(useFreq).A;
    
        % rca replaces NaNs with zeroes, correct this
        nanDims = [1,2]; % if all time points are zero, or all channels are zero
        structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
        noiseVars = {'lowerSideBand','higherSideBand'};
        for z=1:length(structVars)
            if strfind(lower(structVars{z}),'noise')
                for n = 1:length(noiseVars)
                    pphysRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),pphysRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                end
            else   
                pphysRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),pphysRCA(f).(structVars{z}),'uni',false);
            end
        end

        % COMPUTE VALUES FOR PLOTTING
        % AE RCA
%         keepConditions = true;
%         errorType = 'SEM';
%         trialError = false;
%         tempDataStrct = aggregateData(pphysRCA(f).data,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         ampVals(:,f,1:5,:) = tempDataStrct.ampBins;
%         errLB(:,f,1:5,:)   = tempDataStrct.ampErrBins(:,:,:,:,1);
%         errUB(:,f,1:5,:)   = tempDataStrct.ampErrBins(:,:,:,:,2);
%         tempNoiseStrct1 = aggregateData(pphysRCA(f).noiseData.lowerSideBand,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         tempNoiseStrct2 = aggregateData(pphysRCA(f).noiseData.higherSideBand,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         [snrVals(:,f,1:5,:),noiseVals(:,f,1:5,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
% 
%         % COMPARISON
%         tempDataStrct = aggregateData(pphysRCA(f).comparisonData,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         ampVals(:,f,6,:) = tempDataStrct.ampBins;
%         errLB(:,f,6,:) = tempDataStrct.ampErrBins(:,:,:,:,1);
%         errUB(:,f,6,:) = tempDataStrct.ampErrBins(:,:,:,:,2);
%         tempNoiseStrct1 = aggregateData(pphysRCA(f).comparisonNoiseData.lowerSideBand,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         tempNoiseStrct2 = aggregateData(pphysRCA(f).comparisonNoiseData.higherSideBand,pphysRCA(f).settings,keepConditions,errorType,trialError);
%         [snrVals(:,f,6,:),noiseVals(:,f,6,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
        
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

    binVals = cellfun(@(x) str2num(x), pphysRCA(1).settings.binLevels{1});

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
                        
        valSet = pphysRCA(curFreq).stats.Amp;
        errSet1 = pphysRCA(curFreq).stats.ErrLB;
        errSet2 = pphysRCA(curFreq).stats.ErrUB;
        NRset = pphysRCA(curFreq).stats.NR_Params;
        NRmodel = pphysRCA(curFreq).stats.hModel;
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
                    curConds = find(ismember(pphysRCA(curFreq).settings.condsToUse,pphysConds{s}));
                    titleStr = sprintf('horizontal: %s',pphysRCA(curFreq).settings.freqLabels{1}{curFreq});
                else
                    spH = subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(pphysRCA(curFreq).settings.condsToUse,pphysConds{s}));
                    titleStr = sprintf('vertical: %s',pphysRCA(curFreq).settings.freqLabels{1}{curFreq});
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
                        hNR{c} = plot( nrX, nrVals, '-k', 'LineWidth',lWidth);
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
                    figure(freqFig(curFreq));
                    
                    % make new p-values
                    [rcaDataReal,rcaDataImag] = getRealImag(pphysRCA(curFreq).data(curConds,:));
                    freqIdx = pphysRCA(curFreq).settings.freqIndices{3} == curFreq;
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
    if z < 3
        xMin = .13;
        xMax = 3;
        textPos = 0.07;
    else
        xMin = .4;
        xMax = 18;
        textPos = 0.17;
    end
    yUnit = 0.5;
    yMax = 3.5;
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
    set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');
    xlim([xMin,xMax]);
    ylim([0,yMax])
    switch z
        case 1
            plotLabel = 'A';
            title('Horizontal','fontsize',fSize,'fontname','Helvetica')
            text(.15,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.1,'Reference','fontsize',fSize,'fontname','Helvetica')
        case 2
            plotLabel = 'B';
            title('Vertical','fontsize',fSize,'fontname','Helvetica')
            text(.15,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.1,'Reference','fontsize',fSize,'fontname','Helvetica')
        case 3
            plotLabel = 'D';
            title('Horizontal','fontsize',fSize,'fontname','Helvetica')
            ylabel('amplitude (\muV)','fontsize',fSize,'fontname','Helvetica')
            xlabel('displacement (arcmins)','fontsize',fSize,'fontname','Helvetica');
            text(.5,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.1,'No reference','fontsize',fSize,'fontname','Helvetica')
        case 4
            title('Vertical','fontsize',fSize,'fontname','Helvetica')
            plotLabel = 'E';
            text(.5,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.1,'No reference','fontsize',fSize,'fontname','Helvetica')
        otherwise
    end
    text(textPos,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
        plotLabel,'fontsize',fSize*2,'fontname','Helvetica');
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
            set(gca,gcaOpts{:},'xtick',[1.5,3.5],'xticklabel',{'Horizontal','Vertical'},logOpts{:});
             if q == 1
                title('Reference','fontsize',fSize,'fontname','Helvetica');
                text(-.25,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
                    'C','fontsize',fSize*2,'fontname','Helvetica');
            else
                title('No reference','fontsize',fSize,'fontname','Helvetica');
                text(-.25,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
                    'F','fontsize',fSize*2,'fontname','Helvetica');
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

%% MAKE TABLE
T = table([rc_tSqrdSig(:,1),rc_tSqrdP(:,1),rc_tSqrdVal(:,1)], ...
          [rc_tSqrdSig(:,2),rc_tSqrdP(:,2),rc_tSqrdVal(:,2)], ...
          [rc_tSqrdSig(:,3),rc_tSqrdP(:,3),rc_tSqrdVal(:,3)], ...
          [rc_tSqrdSig(:,4),rc_tSqrdP(:,4),rc_tSqrdVal(:,4)], ...
          [rc_tSqrdSig(:,5),rc_tSqrdP(:,5),rc_tSqrdVal(:,5)], ...
          [rc_tSqrdSig(:,6),rc_tSqrdP(:,6),rc_tSqrdVal(:,6)]);
T.Properties.VariableNames = {'HoriInPh','HoriAntiPh','HoriPaired','VertInPh','VertAntiPh','VertPaired'}

