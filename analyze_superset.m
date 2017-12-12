%% Script for running Condition-Specific RCA on Disparity data 

codeFolder = '/Users/kohler/code';
rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
addpath(genpath(rcaCodePath));
addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
setenv('DYLD_LIBRARY_PATH','')

clear all; close all;

%% Set up inputs 

trialError = false;

% we now hardcode the paths to the two data sets
dateStr = datestr(clock,26);
dateStr(strfind(dateStr,'/')) ='';

% make subject list
topFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/';
expList = {'exp1','exp2','exp4','exp8'};
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

saveFilePath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/paper_figures/figure5';

binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= [2,4]; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
condsToUse = [3,7];
trialsToUse = []; %1:10; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false;

%% RC analysis 
condNames = {'HorizontalDisparity', 'VerticalDisparity'};
for f = 1:length(freqsToUse)
    superRCA(f) = rcaSweep(folderNames,binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,1,rcPlotStyle,forceSourceData);
    export_fig(sprintf('%s/SuperSetHoriVert_%0.0f_cov.pdf',saveFilePath,freqsToUse(f)),'-pdf','-transparent',gcf);
    close gcf;
    % RCA replaces NaNs with zero; this is a fix for that
    nanDims = [1,2]; % if all time points are zero, or all channels are zero
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};

    for z = 1:length(structVars)
        if z == 2 || z == 4
            for n = 1:length(noiseVars)
                superRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),superRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
            end
        else
            superRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),superRCA(f).(structVars{z}),'uni',false);
        end

    end
end
    
%% COMPUTE VALUES FOR PLOTTING

keepConditions = true;
errorType = 'SEM';
doNR = false(2,5,2); % 2 freqs, 5 RCs, 2 conditions
doNR(:,1,:) = true; % do fitting for first RC, second and fourth harmonic, all conditions
for f=1:length(freqsToUse)
    if f == 1
        idxList = 1:2;
    else
        idxList = 3:4;
    end
    rcStruct = aggregateData(superRCA(f).data,superRCA(f).settings,keepConditions,errorType,trialError,doNR);    
    compStruct = aggregateData(superRCA(f).comparisonData,superRCA(f).settings,keepConditions,errorType,trialError);
    % note: no NR fitting on noise data
    rcNoiseStrct1 = aggregateData(superRCA(f).noiseData.lowerSideBand,superRCA(f).settings,keepConditions,errorType,trialError,[]);
    rcNoiseStrct2 = aggregateData(superRCA(f).noiseData.higherSideBand,superRCA(f).settings,keepConditions,errorType,trialError,[]);
    compNoiseStrct1 = aggregateData(superRCA(f).comparisonNoiseData.lowerSideBand,superRCA(f).settings,keepConditions,errorType,trialError,[]);
    compNoiseStrct2 = aggregateData(superRCA(f).comparisonNoiseData.higherSideBand,superRCA(f).settings,keepConditions,errorType,trialError,[]);
    % snr
    [rcSNR,rcNoiseVals] = computeSnr(rcStruct,rcNoiseStrct1,rcNoiseStrct2,false);
    [compSNR,compNoiseVals] = computeSnr(compStruct,compNoiseStrct1,compNoiseStrct2,false);
    
    % RC
    superRCA(f).stats.ampVals(:,1:5,:) = rcStruct.ampBins;
    superRCA(f).stats.errLB(:,1:5,:) = rcStruct.ampErrBins(:,:,:,:,1);
    superRCA(f).stats.errUB(:,1:5,:) = rcStruct.ampErrBins(:,:,:,:,2);
    superRCA(f).stats.snrVals(:,1:5,:) = rcSNR;
    superRCA(f).stats.noiseVals(:,1:5,:) = rcNoiseVals;
    superRCA(f).stats.NR_pOpt(:,1:5,:) = rcStruct.NakaRushton.pOpt;
    superRCA(f).stats.NR_JKSE(:,1:5,:) = rcStruct.NakaRushton.JKSE;
    superRCA(f).stats.NR_R2(1:5,:) = rcStruct.NakaRushton.R2;
    superRCA(f).stats.hModel = rcStruct.NakaRushton.hModel;
    superRCA(f).stats.tSqrdP(:,1:5,:) = rcStruct.tSqrdP;
    superRCA(f).stats.tSqrdSig(:,1:5,:) = rcStruct.tSqrdSig;
    superRCA(f).stats.tSqrdVal(:,1:5,:) = rcStruct.tSqrdVal;

    % COMPARISON
    superRCA(f).stats.ampVals(:,6,:) = compStruct.ampBins;
    superRCA(f).stats.errLB(:,6,:) = compStruct.ampErrBins(:,:,:,:,1);
    superRCA(f).stats.errUB(:,6,:) = compStruct.ampErrBins(:,:,:,:,2);
    superRCA(f).stats.snrVals(:,6,:) = compSNR;
    superRCA(f).stats.noiseVals(:,6,:) = compNoiseVals;
    superRCA(f).stats.NR_pOpt(:,6,:) = compStruct.NakaRushton.pOpt;
    superRCA(f).stats.NR_JKSE(:,6,:) = compStruct.NakaRushton.JKSE;
    superRCA(f).stats.NR_R2(6,:) = compStruct.NakaRushton.R2;
    superRCA(f).stats.tSqrdP(:,6,:) = compStruct.tSqrdP;
    superRCA(f).stats.tSqrdSig(:,6,:) = compStruct.tSqrdSig;
    superRCA(f).stats.tSqrdVal(:,6,:) = compStruct.tSqrdVal;
end
delete(gcp('nocreate'));

%% MAKE FIGURE
close all
lWidth = 1.5;
fSize = 12;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
cBrewer = load('colorBrewer.mat');
mainColors = [cBrewer.rgb20(5,:); cBrewer.rgb20(7,:)];

topoVals = [superRCA(1).A(:,1);superRCA(2).A(:,1)];
rcaColorBar = [min(topoVals),max(topoVals)];
newExtreme = round(max(abs(rcaColorBar(:,f)))*5)./5;
rcaColorBar = [-newExtreme,newExtreme*1.001];
plotLabel = {'A','B','C','D'};
for f=1:length(freqsToUse)
    binVals = cellfun(@(x) str2num(x), superRCA(f).settings.binLevels{1});
    valSet = squeeze(superRCA(f).stats.ampVals(:,1,:));
    errSet1 = squeeze(superRCA(f).stats.errLB(:,1,:));
    errSet2 = squeeze(superRCA(f).stats.errUB(:,1,:));
    NRset = squeeze(superRCA(f).stats.NR_pOpt(:,1,:));
    NRmodel = superRCA(f).stats.hModel;
    maxNoise = max(superRCA(f).stats.noiseVals(:,1,:),[],3)';
    
    egiH(f) = subplot(2,2,2+(f-1)*2);
    hold on
    if f == 1
        [figH(f),cH(f)] = mrC.plotOnEgi(superRCA(f).A(:,1)*-1,rcaColorBar,true);
    else
        [figH(f),cH(f)] = mrC.plotOnEgi(superRCA(f).A(:,1),rcaColorBar,true);
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
        valH(c)=plot(binVals,valSet(:,c),ampMarkerStyle{:});
        hE = ErrorBars(binVals,valSet(:,c),[errSet1(:,c),errSet2(:,c)],'color',mainColors(c,:),'type','bar','cap',false,'barwidth',lWidth);
        uistack(valH(c),'bottom')
        cellfun(@(x) uistack(x,'bottom'), hE);
        hold on
        
        % plot Naka-Rushton
        nFine = 1e2;
        nrX = linspace( min(binVals), max(binVals), nFine )';
        nrVals = NRmodel( nrX, NRset(:,c));
        hNR{c} = plot( nrX, nrVals, '-','color',mainColors(c,:), 'LineWidth',lWidth);
    end
    cellfun(@(x) uistack(x,'bottom'), hNR);

    if f == 1     
        yUnit = 1;
        yMax = 5;
        legend(valH,{'horizontal','vertical'},'fontsize',fSize,'fontname','Helvetica','location','northwest');
        legend boxoff
        title('Second Harmonic','fontsize',fSize,'fontname','Helvetica')
    else
        yUnit = 0.5;
        yMax = 1.5;
        ylabel('amplitude (\muV)','fontsize',fSize,'fontname','Helvetica')
        xlabel('displacement (arcmins)','fontsize',fSize,'fontname','Helvetica');
        title('Fourth Harmonic','fontsize',fSize,'fontname','Helvetica')
    end
                
    xMin = .4;
    xMax = 18;
    
    set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0,0.2,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top');
          
    xlim([xMin,xMax]);
    ylim([0,yMax])
    text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
        plotLabel(f),'fontsize',fSize*2,'fontname','Helvetica');
    text(30,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,...
        plotLabel(f+2),'fontsize',fSize*2,'fontname','Helvetica');

    % plot noise patch
    yNoiseVals = [0,maxNoise(1),maxNoise,maxNoise(end),0]; % start and end points just repeats of first and last
    xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
    pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
    uistack(pH,'bottom') 
    hold off
end



for f = 1:length(freqsToUse)
    addX = 0.7;
    addY = 0.7;
    if f == length(freqsToUse)
        set(cH(f),'fontsize',fSize,'fontname','Helvetica','YTick',linspace(min(rcaColorBar),min(rcaColorBar)*-1,5));
        ylabel(cH(f),'weights','fontsize',fSize,'fontname','Helvetica')
        set(cH(f),'location','eastoutside');
        set(cH(f),'units','centimeters');
        %cBarPos = get(cH(e,f),'position');
        cBarPos = [14.5,4.5,.4,4];
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
        newPos(1) = newPos(1)-(newPos(3)*addX*.7);
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
export_fig(sprintf('%s/figure5_combined.pdf',saveFilePath),'-pdf','-transparent',gcf);