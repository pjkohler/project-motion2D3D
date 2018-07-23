%% Script for running Condition-Specific RCA on Disparity data 

codeFolder = '/Users/kohler/code';
rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
addpath(genpath(rcaCodePath));
addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
setenv('DYLD_LIBRARY_PATH','')

clear all; close all;

%% Set up inputs 
plotSupplemental = false;
projectedData = true;

% we now hardcode the paths to the two data sets
dateStr = datestr(clock,26);
dateStr(strfind(dateStr,'/')) ='';

% make subject list
topFolder = '/Volumes/svndl/FinishedExperiments/2018_Kohler_NatureCommunications/';
if plotSupplemental
    expList = {'exp4','exp5'};
    condsToUse = [4,8];   
else
    expList = {'exp1','exp2','exp3','exp4'};
    condsToUse = [3,7];
end
idList = [];
for e = 1:length(expList);
    expFolder = subfolders(sprintf('%s/*%s*',topFolder,expList{e}),1);
    if e == 1
        subjList = subfolders(sprintf('%s/*nl-*',expFolder{1}),1);
    else
    end
    tempList = subfolders(sprintf('%s/*nl-*',expFolder{1}),1);
    for s = 1:length(tempList);
        [~,curSub]=fileparts(tempList{s});
        idIdx = strfind(curSub,'nl-');
        curSub = curSub(idIdx:(idIdx+6));
        if e == 1
            idList = [idList;{curSub}];
        else
            if any(~cell2mat(cellfun(@(x) isempty(strfind(x,curSub)),subjList,'uni',false)))
                continue;
            else
                subjList = [subjList;tempList(s)];
                idList = [idList;{curSub}];
            end
        end
    end 
end

% go deeper into the folder names
for f = 1:length(subjList)
    tempFolders = subfolders(subjList{f},1);
    tempFolders = tempFolders(cellfun(@(x) isempty(strfind(x,'mff')),tempFolders));
    folderNames{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',tempFolders{end});
end

saveFilePath = '/Volumes/svndl/FinishedExperiments/2018_Kohler_NatureCommunications/figures/paper_figures/figure5';

binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= [2,4]; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
trialsToUse = []; %1:10; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false;
trialError = false;

%% RC analysis 
condNames = {'HorizontalDisparity', 'VerticalDisparity'};
for f = 1:length(freqsToUse)
    superRCA(f) = rcaSweep(folderNames,binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,1,rcPlotStyle,forceSourceData);
    export_fig(sprintf('%s/SuperSetHoriVert_%0.0f_cov.pdf',saveFilePath,freqsToUse(f)),'-pdf','-transparent',gcf);
    close gcf;
end
zeroData = cellfun(@(x) any(x(:)==0), cat(1,superRCA(:).data));
    
if any(zeroData(:))
    error('data values exactly zero, this should not happen');
else
end
    
%% COMPUTE VALUES FOR PLOTTING

keepConditions = true;
errorType = 'SEM';
doNR = false(2,6,2); % 2 freqs, 5 RCs, 2 conditions
doNR(:,1,:) = true; % do fitting for first RC, second and fourth harmonic, all conditions

for f = 1:length(superRCA)
        fprintf('\n ... rc no. %0.0f ...\n',f);
        rcStruct = aggregateData(superRCA(f),keepConditions,errorType,trialError,doNR);
        % RC
        superRCA(f).stats.Amp = squeeze(rcStruct.ampBins);
        superRCA(f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp);
        superRCA(f).stats.ErrLB = squeeze(rcStruct.ampErrBins(:,:,:,:,1));
        superRCA(f).stats.ErrUB = squeeze(rcStruct.ampErrBins(:,:,:,:,2));
        superRCA(f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins);
        superRCA(f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise);
        % Naka-Rushton
        superRCA(f).stats.NR_Params = squeeze(rcStruct.NakaRushton.Params);
        superRCA(f).stats.NR_R2 = squeeze(rcStruct.NakaRushton.R2);
        superRCA(f).stats.NR_JKSE = squeeze(rcStruct.NakaRushton.JackKnife.SE);
        superRCA(f).stats.NR_JKParams = squeeze(rcStruct.NakaRushton.JackKnife.Params);
        superRCA(f).stats.hModel = rcStruct.NakaRushton.hModel;
        % t-values
        superRCA(f).stats.tSqrdP = squeeze(rcStruct.tSqrdP);
        superRCA(f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig);
        superRCA(f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal);
end
% shut down parallel pool, which was used for fitting Naka-Rushton
delete(gcp('nocreate'));
clc;

%% MAKE FIGURE
close all
rcNum = 1;
lWidth = 1.5;
fSize = 12;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
cBrewer = load('colorBrewer.mat');
mainColors = [cBrewer.rgb20(5,:); cBrewer.rgb20(7,:)];

topoVals = [superRCA(1).A(:,rcNum);superRCA(2).A(:,rcNum)];
rcaColorBar = [min(topoVals),max(topoVals)];
newExtreme = round(max(abs(rcaColorBar(:,f)))*5)./5;
%rcaColorBar = [-newExtreme,newExtreme*1.001];
rcaColorBar = [-0.4,0.4001];
plotLabel = {'A','B','C','D'};
if plotSupplemental
    flipIdx = [1,1];
else
    flipIdx = [1,1];
end
for f=1:length(freqsToUse)
    binVals = cellfun(@(x) str2num(x), superRCA(f).settings.binLabels);
    logStep = diff(reallog(binVals(1:2))); % step size
    xMin = reallog(binVals(1))-logStep*.5;
    xMax = reallog(binVals(end))+logStep*2.5; % add 2.5 steps
    extraBins = arrayfun(@(x) exp(reallog(binVals(end))+x), [logStep,logStep*2]);
    
    % make Naka-Rushton values
    NRmodel = superRCA(f).stats.hModel;
    NRset = squeeze(superRCA(f).stats.NR_Params(:,rcNum,:));          
    NRerrs = superRCA(f).stats.NR_JKSE;
    NRmodel = superRCA(f).stats.hModel;
    
    % compute new signal values, averaged over bins
    [rcaDataReal,rcaDataImag] = getRealImag(superRCA(f).data);
    rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataReal,'uni',false);
    rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
    rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataImag,'uni',false);
    rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
    realBinMean = squeeze(nanmean(rcaDataReal));
    imagBinMean = squeeze(nanmean(rcaDataImag));
    % compute new noise values, averaged over bins
    [noiseLoReal,noiseLoImag] = getRealImag(superRCA(f).noiseData.lowerSideBand);
    [noiseHiReal,noiseHiImag] = getRealImag(superRCA(f).noiseData.lowerSideBand);
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
        valSet = squeeze(superRCA(f).stats.Amp(:,rcNum,:));
        errSet1 = squeeze(superRCA(f).stats.ErrLB(:,rcNum,:));
        errSet2 = squeeze(superRCA(f).stats.ErrLB(:,rcNum,:));
        % compute vector means of mean bins, and errors
        valSet(length(binVals)+1,:) = sqrt(nanmean(realBinMean,1).^2+nanmean(imagBinMean,1).^2);
        
        % store p-values for later
        
        % compute elliptical error and do Hotelling's T2 against zero
        for c=1:length(condsToUse)
            rc_tSqrdP(1:length(binVals),c+(f-1)*3) = superRCA(f).stats.tSqrdP(:,rcNum,c);
            rc_tSqrdVal(1:length(binVals),c+(f-1)*3) = superRCA(f).stats.tSqrdVal(:,rcNum,c);
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
    noiseSet = max(superRCA(f).stats.NoiseAmp(:,rcNum,:),[],3);
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
    
    % do Naka-Rushton, paired tests
    testVal = squeeze(superRCA(f).stats.NR_JKParams(:,:,rcNum,:));
    testVal = permute(testVal,[2,1,3]); % move subjects to first dim    
    paramIdx = [1,2,3,4];% only look at c50 and rMax
    jkDf = size(testVal,1)-1;
    diffErr = jackKnifeErr(testVal(:,paramIdx,1)-testVal(:,paramIdx,2));
    grandDiff = NRset(paramIdx,1) - NRset(paramIdx,2);
    paramPairedT(:,f) = grandDiff'./diffErr;
    paramPairedP(:,f) = 2*tcdf( -abs(paramPairedT(:,f)) , jkDf);
    
    egiH(f) = subplot(2,2,2+(f-1)*2);
    hold on
    if f == 1
        [figH(f),cH(f)] = mrC.plotOnEgi(superRCA(f).A(:,rcNum)*flipIdx(f),rcaColorBar,true);
    else
        [figH(f),cH(f)] = mrC.plotOnEgi(superRCA(f).A(:,rcNum)*flipIdx(f),rcaColorBar,true);
    end
    hold off
    figH(f+2) = subplot(2,2,f+(f-1));
    hold on
    for c = 1:length(condsToUse)
        if c == 1
            ampMarkerStyle = {'o','LineWidth',lWidth,'Color',mainColors(c,:),'markerfacecolor',[1 1 1],'MarkerSize',7}; % 'MarkerEdgeColor','none'
        else
            ampMarkerStyle = {'sq','LineWidth',lWidth,'Color',mainColors(c,:),'markerfacecolor',[1 1 1],'MarkerSize',7}; % 'MarkerEdgeColor','none'
        end
        xVals = reallog([binVals;extraBins((condsToUse(c)>4)+1)]);
        yVals = valSet(:,c);
        errVals = [errSet1(:,c),errSet2(:,c)];
        valH(c)=plot(xVals,yVals,ampMarkerStyle{:});
        hE = ErrorBars(xVals,yVals,errVals,'color',mainColors(c,:),'type','bar','cap',false,'barwidth',lWidth);
        uistack(valH(c),'bottom')
        cellfun(@(x) uistack(x,'bottom'), hE);
        hold on
        
        % plot Naka-Rushton
        nFine = 1e2;
        nrX = linspace( min(binVals), max(binVals), nFine )';
        nrVals = NRmodel( nrX, NRset(:,c));
        hNR{c} = plot( reallog(nrX), nrVals, '-','color',mainColors(c,:), 'LineWidth',lWidth);
    end
    cellfun(@(x) uistack(x,'bottom'), hNR);
    
    if f == 1     
        yUnit = 1;
        yMin = 0;
        yMax = 5;
        legend(valH,{'horizontal','vertical'},'fontsize',fSize,'fontname','Helvetica','location','northwest');
        legend boxoff
        titleStr = '\it\fontname{Helvetica}2F ';
        title(titleStr,'fontsize',fSize,'fontname','Helvetica','interpreter','tex');
    else
        yUnit = 0.5;
        yMin = 0;
        yMax = 1.5;
        ylabel('amplitude (\muV)','fontsize',fSize,'fontname','Helvetica')
        xlabel('displacement (arcmins)','fontsize',fSize,'fontname','Helvetica');
        titleStr = '\it\fontname{Arial}4F ';
        title(titleStr,'fontsize',fSize,'fontname','Helvetica','interpreter','tex');
    end
    set(gca,gcaOpts{:},'XMinorTick','off','xtick',reallog([0.2,0.5,1,2,4,8,16]),'xticklabels',arrayfun(@(x) num2str(x),([0.2,0.5,1,2,4,8,16]),'uni',false),'ytick',0:yUnit:yMax,'Layer','top');
          
    xlim([xMin,xMax]);
    ylim([yMin,yMax])
    text(xMin-logStep*.75,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,...
        plotLabel(f),'fontsize',fSize*2,'fontname','Helvetica');
    text(xMax+logStep*1.5,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,...
        plotLabel(f+2),'fontsize',fSize*2,'fontname','Helvetica');

    % split point between bin and average vals
    xSplit = reallog(binVals(end))+logStep*.5;
    plot(ones(2,1)*xSplit,[0,yMax],'k','LineWidth',lWidth)
    % plot noise patch, 
    % add mean bin noise twice to the end with extra bin
    % we are plotting the mean values over two x-axis points
    xNoiseVals = [xMin,xMin,reallog(binVals'),xSplit,xSplit];
    yNoiseVals = [0,noiseSet(1),noiseSet(1:end-1)',noiseSet(end-1),0]; % start and end points just repeats of first and last
    pH(1) = fill(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
    % additional vals for averages
    xNoiseVals = [xSplit,xSplit,reallog(extraBins),xMax,xMax];
    yNoiseVals = [0,repmat(noiseSet(end),1,4),0];
    pH(2) = fill(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
    arrayfun(@(x) uistack(x,'bottom'), pH);
end



for f = 1:length(freqsToUse)
    addX = 0.38;
    addY = 0.38;
    if f == length(freqsToUse)
        set(cH(f),'fontsize',fSize,'fontname','Helvetica','YTick',linspace(min(rcaColorBar),min(rcaColorBar)*-1,5));
        ylabel(cH(f),'weights','fontsize',fSize,'fontname','Helvetica')
        set(cH(f),'location','eastoutside');
        set(cH(f),'units','centimeters');
        %cBarPos = get(cH(e,f),'position');
        cBarPos = [14,4.5,.4,4];
        set(cH(f),'position',cBarPos);
        oldPos = newPos;
        newPos = get(egiH(f),'position');
        newPos(1) = oldPos(1);
        newPos(2) = newPos(2)-(newPos(4)*addY*.6);
        newPos(3) = oldPos(3);
        newPos(4) = oldPos(4);
    else
        set(cH(f),'visible','off');
        newPos = get(egiH(f),'position');
        newPos(1) = newPos(1)-(newPos(3)*addX);
        newPos(2) = newPos(2)-(newPos(4)*addY*.6);
        newPos(3) = newPos(3)*(1+addX);
        newPos(4) = newPos(4)*(1+addY);
    end
    set(egiH(f),'position',newPos); 
end
drawnow;
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = figPos(4)/figPos(3)*17.8;
figPos(3) = 17.8;
set(gcf,'pos',figPos);
if projectedData
    suffix = 'proj';
else
    suffix = 'multi';
end

if plotSupplemental
    export_fig(sprintf('%s/suppl_figure5_combined_%s.pdf',saveFilePath,suffix),'-pdf','-transparent',gcf);
else
    export_fig(sprintf('%s/figure5_combined_%s.pdf',saveFilePath,suffix),'-pdf','-transparent',gcf);
end

if ~projectedData
    return
else
end

%% MAKE TABLE
%if plotSupplemental
%    paperText = 'In 2nd harmonic data from Experiments 4 and 5, ';
%else
%    paperText = 'In 2nd harmonic data from Experiments 1,2,4 and 8, ';
%end

if ~all(rc_tSqrdDF(1,:) == rc_tSqrdDF(1,1))
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


readyArray = [stringBinVals, stringPairedP(:,3),stringPairedT(:,3),stringPairedD(:,1)...
                stringPairedP(:,6),stringPairedT(:,6),stringPairedD(:,2)]';
rowLabels = {'disp (arcmins)','2F p','2F t-value','2F Cohen''s D','4th p','4F t-value','4F Cohen''s D'}';

readyArray = cat(2,rowLabels,readyArray);

if plotSupplemental
    titleStr = 'Exp. 4-5 (IOVD)';
    tableIdx = 2;
else
    titleStr = 'Exp. 1-4 (CDOT+IOVD)';
    tableIdx = 1;
end

varLabels = cat(2,{num2str(rc_tSqrdDF(1),'df=%0.0f')},...
              arrayfun(@(x) num2str(x,'bin%0.0f'),1:length(binVals),'uni',false),...
              {'ave'});
titleLabels = cell(size(varLabels));
titleLabels{1} = titleStr;

readyArray = cat(1,titleLabels,varLabels,readyArray);
table_file = sprintf('%s/superset_table%0.0f.csv',saveFilePath,tableIdx);
superTable = array2table(readyArray);
writetable(superTable,table_file,'WriteRowNames',false,'WriteVariableNames',false);

if exist(sprintf('%s/superset_table1.csv',saveFilePath)) && exist(sprintf('%s/superset_table2.csv',saveFilePath))
    combined_file = sprintf('%s/superset_combined.csv',saveFilePath);
    if exist(combined_file,'file')
        delete(combined_file);
    else
    end
    catCmd{1} = sprintf('cd %s',saveFilePath);
    catCmd{2} = 'cat superset_table1.csv superset_table2.csv > superset_combined.csv';
    system(strjoin(catCmd, '; '));
else
end


    
