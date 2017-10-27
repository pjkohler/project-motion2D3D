%function draft_analyze_babymotion2D3D(doExp,trialError)
%% Add Paths
close all;
if ~exist('codeFolder','var')
    codeFolder = '/Users/labmanager/Desktop/LabManager/MatAnal';
    rcaCodePath = sprintf('%s/rcaBase',codeFolder);
    addpath(genpath(codeFolder));
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/mrC',codeFolder)));
    addpath(genpath(sprintf('%s/matlab_lib',codeFolder)));
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
    addpath(genpath(sprintf('%s/export_fig',codeFolder)));
    
else
end
setenv('DYLD_LIBRARY_PATH','')

%% TAKE CARE OF INPUT ARGUMENTS
%if nargin < 1
    doExp = 1;
%else
%end
%if nargin < 2
    trialError = true;
%else
%end

%saveLocation = sprintf('/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/exp%d',doExp);
saveLocation = '/Users/labmanager/Desktop/LabManager/WM_Data/2017_2DBaby/2DBaby_Data';
saveFileName = sprintf('%s/rcaData.mat',saveLocation);
load(saveFileName);

%% COMPUTE VALUES FOR PLOTTING
keepConditions = true;
errorType = 'SEM';
doNR = true(2,6,4); % 5 freqs, 5 RCs + OZ, 4 conditions

for f=1:length(freqsToUse)+1
    fprintf('%d\n',f);
    if f==length(freqsToUse)+1
        curRCA = allRCA;
        curIdx = 3:4;

    else
        curRCA = freqRCA(f);
        curIdx = f;

    end
    tempDataStrct = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,trialError,doNR);
    ampVals(:,curIdx,1:5,:) = tempDataStrct.ampBins;
    errLB(:,curIdx,1:5,:)   =tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,curIdx,1:5,:)   =tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    [snrVals(:,curIdx,1:5,:),noiseVals(:,curIdx,1:5,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
    NR_pOpt(:,curIdx,1:5,:) = tempDataStrct.NakaRushton.pOpt;
    NR_JKSE(:,curIdx,1:5,:) = tempDataStrct.NakaRushton.JKSE;
    NR_R2(curIdx,1:5,:) = tempDataStrct.NakaRushton.R2;
    if ~isempty(tempDataStrct.NakaRushton.hModel) && ~exist('hModel','var')
        hModel = tempDataStrct.NakaRushton.hModel;
    else
    end
    delete(gcp('nocreate'));
    
    % COMPARISON
    tempDataStrct = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,trialError,doNR);
    ampVals(:,curIdx,6,:) = tempDataStrct.ampBins;
    errLB(:,curIdx,6,:)=tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,curIdx,6,:)=tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    [snrVals(:,curIdx,6,:),noiseVals(:,curIdx,6,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
    NR_pOpt(:,curIdx,6,:) = tempDataStrct.NakaRushton.pOpt;
    NR_JKSE(:,curIdx,6,:) = tempDataStrct.NakaRushton.JKSE;
    NR_R2(curIdx,6,:) = tempDataStrct.NakaRushton.R2;
    delete(gcp('nocreate'));
end

for f= length(allRCA.data(:,1)):-1:1
    if isnan(allRCA.data{f,1}(1,1,1))
        allRCA.data(f,:)=[];
        allRCA.comparisonData(f,:)=[];
    end
end

freqRCA(3:4) = allRCA;

% repeat allRCA, this is only done to preserve topographies in the code below
%freqRCA(9:12) = allRCA; % repeat allRCA, this is only done to preserve topographies in the code below


%% Plot RCs

    %% PLOT RCs
    close all;
    nFreq = length(allRCA.settings.freqsToUse);
    nComp = allRCA.settings.nComp;
    plotSNR = false;
    rcaType = 'all';
    % plot settings
    % set figure size in the beginning
    figHeight = 30;
    figWidth = 20;

    binVals = cellfun(@(x) str2num(x), allRCA.settings.binLevels{1});

    lWidth = .75;
    errorWidth = .3;
    cBrewer = load('colorBrewer.mat');
    color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
    color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
    subColors = repmat([color1; color2],2,1);
    subColors = subColors(condsToUse,:);
    fSize = 12;
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','color','none','fontsize',fSize,'fontname','Arial','linewidth',lWidth};
    condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
    condLabels = condLabels(condsToUse);
    xMin = 1; xMax = 40;
    
    
    % Get parameter standard errors
%     if logScale
%         % put c50 back in original units
%         pA1(1,:) = 10.^pA1(1,:) / fx;
%         pA2(1,:) = 10.^pA2(1,:) / fx;
%         pN( 1,:) = 10.^pN( 1,:) / fx;
%         pA1jk(1,:) = 10.^pA1jk(1,:) / fx;
%         pA2jk(1,:) = 10.^pA2jk(1,:) / fx;
%         pNjk( 1,:) = 10.^pNjk( 1,:) / fx;
%         fprintf( '\nFIT LOG10( %g * OFFSET )\n', fx )
%     else
%         fprintf( '\nFIT LINEAR OFFSET\n' )
%     end

    

    for f= 1:length(freqsToUse)
        if strcmp(rcaType,'all'); 
            curFreq = f+2;
        else
            curFreq = f;
        end
        figure;
        for r = 1:(nComp+1)

            if r<6
                egiH(r) = subplot(nComp+1,3,3+(r-1)*3);
                hold on
                rcaColorBar = [min(freqRCA(curFreq).A(:,r)),max(freqRCA(curFreq).A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(freqRCA(curFreq).A(:,r),rcaColorBar);
                hold off
            else
            end
            %title(['Full RCA ' num2str(r)],'fontsize',fSize,'fontname','Arial');
            for s=1:1 % vertical or horizontal motion
                if s==1
                    sH(f,r,1) = subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(condsToUse,[1,3]));
                    titleStr = sprintf('horizontal: %s',allRCA.settings.freqLabels{1}{f});
%                 else
%                     sH(f,r,2) = subplot(nComp+1,3,2+(r-1)*3);
%                     curConds = find(ismember(condsToUse,5:8));
%                     titleStr = sprintf('vertical: %s',fullRCA.settings.freqLabels{f});
                end
                hold on
                for c=1:length(curConds)
                    if plotSNR
                        valSet = snrVals(:,curFreq,r,:);
                        ampH(c)=plot(binVals,snrVals(:,curFreq,r,curConds(c)),'o','LineWidth',lWidth*1.5,'Color',subColors(curConds(c),:));
                        plot(binVals,noiseVals(:,curFreq,r,curConds(c)),'sq','Color',subColors(curConds(c),:),'MarkerSize',5);
                    else
                        valSet = ampVals(:,curFreq,r,:);
                        ampH(c) = errorbar(binVals,ampVals(:,curFreq,r,curConds(c)), errLB(:,curFreq,r,curConds(c)), errUB(:,curFreq,r,curConds(c)),'.', 'Color', subColors(curConds(c),:),'MarkerSize', 12, 'LineWidth',errorWidth);
                        %cellfun(@(x) uistack(x,'bottom'), ampH(c));
                        hold on
                    end
                     if ~isnan(NR_pOpt(1,curFreq,r,curConds(c)))
                         % plot Naka-Rushton
                         nFine = 1e2;
                         nrX = linspace( min(binVals), max(binVals), nFine )';
                         nrVals = hModel( nrX, NR_pOpt(:,curFreq,r,curConds(c)));
                         plot( nrX, nrVals, '-', 'Color', subColors(curConds(c),:),'LineWidth',lWidth);
                         
                         
                     else
                     end
                end
             

                if freqsToUse(f) == 2 && r == 1
                    yUnit = 3;
                    yMax = 12.0;
                elseif freqsToUse(f) == 2 && r~=1 && r~=6
                    yUnit = .75;
                    yMax = 3;
                elseif freqsToUse(f) == 1 && r ==2
                    yUnit = 1.5;
                    yMax = 6;
                elseif freqsToUse(f) == 2 && r == (nComp +1)
                    yUnit = 1.5;
                    yMax = 6;
                elseif freqsToUse(f) == 3
                    yUnit = .5;
                    yMax = 2;
                elseif freqsToUse(f) == 4 && r == 1
                    yUnit = .75;
                    yMax =3;
                elseif freqsToUse(f) == 4 && r == 2 
                    yUnit = .5;
                    yMax = 2;
                elseif freqsToUse(f) == 4 && r == (nComp +1)
                    yUnit = .5;
                    yMax = 2;
                elseif freqsToUse(f) == 4 
                    yUnit = .25;
                    yMax = 1;
                else
                    yUnit = 1;
                    yMax = 4;
                end

                %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
                %yMax = ceil(max(valSet(:)))+yUnit;
                set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0.1,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top','clipping','off');
                xlim([xMin,xMax]);
                ylim([0,yMax])
                % plot noise patch
                meanNoise = max(noiseVals(:,curFreq,r,curConds),[],4)';
                yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
                xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
                pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
                uistack(pH,'bottom')
%changing to .1 to see if this improves plotting, .0008 x
                if s==1
                    if r > nComp
                        text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,'OZ','fontsize',fSize,'fontname','Arial');
                    else
                        text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,['RC ',num2str(r)],'fontsize',fSize,'fontname','Arial');
                    end
                else
                end
                if r==1
                    title(titleStr,'fontsize',fSize,'fontname','Arial');
                elseif r==6
                    if s==1
                        ylabel('Amplitude (\muV)')
                        xlabel('Distance (arcmins)');
                    %else
                        lH = legend(ampH,condLabels(curConds),'location','northeast');
                        legend boxoff
                        lPos = get(lH,'position');
                        lPos(1) = lPos(1) + .2;
                        lPos(2) = lPos(2) + .05;
                        set(lH,'position',lPos);            
                        tH = text(60,-0.1,sprintf('n = %0d',size(allRCA.data,2)),'fontsize',fSize,'fontname','Arial');
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
            newPos(1) = newPos(1)-(newPos(3)*addX*.5);
            newPos(2) = newPos(2)-(newPos(4)*addY*.4);
            newPos(3) = newPos(3)*(.5+addX);
            newPos(4) = newPos(4)*(.5+addY);
            set(egiH(r),'position',newPos);
        end
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = figHeight;
        figPos(3) = figWidth;
        set(gcf,'pos',figPos);
        if plotSNR        
            export_fig(sprintf('%s/exp%d_%dF1_%s_snr.pdf',saveLocation,doExp,freqsToUse(f),rcaType),'-pdf','-transparent',gcf);
        else
            export_fig(sprintf('%s/exp%dcond13_%dF1_%s',saveLocation,doExp,freqsToUse(f),rcaType),'openegl','-png','-eps','-pdf','-transparent',gcf);
            print(sprintf('2DBabycond13_%dF1',freqsToUse(f)),'-dsvg');
        end
        if f < 3
            for v = 1:1
                for z=1:5;
                    fH = figure;
                    set(fH, 'units', 'centimeters');
                    figPos = get(fH,'pos');
                    figPos(4) = figHeight/4.25;
                    figPos(3) = figWidth*(1/2.75);
                    set(fH,'pos',figPos);
                    aH = copyobj(sH(f,1,v),fH);
                    set(aH,'position',[.2,.2,.7,.6]);
                    dataObjs = get(aH, 'Children'); %handles to low-level graphics objects in axes
                    lineObjs = findobj(dataObjs, 'type', 'line');
                    if z == 1
                         keepIdx = [];
                         delete(lineObjs);
                    elseif z < 5
                        keepIdx = [keepIdx,6-z,(1:10)+10*(z-2)+4];
                        deleteIdx = true(1,length(lineObjs));
                        deleteIdx(keepIdx) = false;
                        delete(lineObjs(deleteIdx));
                    else
                    end
                    ylabel('Amplitude (\muV)');
                    xlabel('Distance (arcmins)');
                    if ~exist(sprintf('%s/exp%d_talkFigures',saveLocation,doExp),'dir')
                        mkdir(sprintf('%s/exp%d_talkFigures',saveLocation,doExp));
                    else
                    end
                    if v == 1
                        %savefig(sprintf('%s/exp%d_talkFigures/exp%d_rc%d_%s_hori%d.fig',saveLocation,doExp,doExp,f,rcaType,z));
                        export_fig(sprintf('%s/exp%d_talkFigures/exp%dcond13_%dF1_%s_hori%d',saveLocation,doExp,doExp,freqsToUse(f),rcaType,z),'-opengl','-png', '-pdf', '-eps', '-transparent', gcf);
                        %print(sprintf('2DBabycond13_talkFigures%dF1_vert%d',freqsToUse(f),z),'-dsvg');

                    end
                    close(fH);
                end
            end
        else
        end
    end    
   
save('BabyDataOutput');

