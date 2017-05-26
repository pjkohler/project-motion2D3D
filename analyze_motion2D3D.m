function analyze_motion2D3D(doExp,trialError)
    %% ADD PATHS
    close all;
    codeFolder = '/Users/kohler/code';
    rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
    setenv('DYLD_LIBRARY_PATH','')
    
    %% TAKE CARE OF INPUT ARGUMENTS
    if nargin < 1
        doExp = 1;
    else
    end
    if nargin < 2
        trialError = true;
    else
    end
    
    saveLocation = sprintf('/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/exp%d',doExp);
    saveFileName = sprintf('%s/rcaData.mat',saveLocation);
    load(saveFileName);

    %% COMPUTE VALUES FOR PLOTTING
    keepConditions = true;
    errorType = 'SEM';
    for f=1:6
        if f==6
            curRCA = allRCA;
            curIdx = 9:12;
        elseif f==5
            curRCA = fullRCA;
            curIdx = 5:8;
        else
            curRCA = freqRCA(f);
            curIdx = f;
        end
        tempDataStrct = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,trialError);
        ampVals(:,curIdx,1:5,:) = tempDataStrct.ampBins;
        errLB(:,curIdx,1:5,:)   =tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
        errUB(:,curIdx,1:5,:)   =tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
        tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError);
        tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError);
        [snrVals(:,curIdx,1:5,:),noiseVals(:,curIdx,1:5,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);

        % COMPARISON
        tempDataStrct = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,trialError);
        ampVals(:,curIdx,6,:) = tempDataStrct.ampBins;
        errLB(:,curIdx,6,:)=tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
        errUB(:,curIdx,6,:)=tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
        tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError);
        tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError);
        [snrVals(:,curIdx,6,:),noiseVals(:,curIdx,6,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
    end

    freqRCA(5:8) = fullRCA; % repeat fullRCA, this is only done to preserve topographies in the code below
    freqRCA(9:12) = allRCA; % repeat allRCA, this is only done to preserve topographies in the code below


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

    binVals = cellfun(@(x) str2num(x), fullRCA.settings.binLevels{1});

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
    condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
    condLabels = condLabels(condsToUse);
    xMin = .4; xMax = 18;

    for f=1:nFreq
        if strcmp(rcaType,'all'); 
            curFreq = f+8;
        elseif strcmp(rcaType,'full'); 
            curFreq = f+4;
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
            for s=1:2 % vertical or horizontal motion
                if s==1
                    sH(f,r,1) = subplot(nComp+1,3,1+(r-1)*3);
                    curConds = find(ismember(condsToUse,1:4));
                    titleStr = sprintf('horizontal: %s',fullRCA.settings.freqLabels{f});
                else
                    sH(f,r,2) = subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(condsToUse,5:8));
                    titleStr = sprintf('vertical: %s',fullRCA.settings.freqLabels{f});
                end
                hold on
                for c=1:length(curConds)
                    if plotSNR
                        valSet = snrVals(:,curFreq,r,:);
                        ampH(c)=plot(binVals,snrVals(:,curFreq,r,curConds(c)),'-','LineWidth',lWidth*1.5,'Color',subColors(curConds(c),:));
                        plot(binVals,noiseVals(:,curFreq,r,curConds(c)),'sq','Color',subColors(curConds(c),:),'MarkerSize',5);
                    else
                        valSet = ampVals(:,curFreq,r,:);
                        ampH(c)=plot(binVals,ampVals(:,curFreq,r,curConds(c)),'-','LineWidth',lWidth*1.5,'Color',subColors(curConds(c),:));
                        %plot(binVals,noiseVals(:,curFreq,r,curConds(c)),'sq','Color',subColors(curConds(c),:),'MarkerSize',5);
                        %errorbar(binVals,ampVals(:,curFreq,r,curConds(c)),errLB(:,curFreq,r,curConds(c)),errUB(:,curFreq,r,curConds(c)),'Color',subColors(curConds(c),:),'LineWidth',lWidth);
                        hE = ErrorBars(binVals',ampVals(:,curFreq,r,curConds(c)),[errLB(:,curFreq,r,curConds(c)),errUB(:,curFreq,r,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth*1.5);
                        cellfun(@(x) uistack(x,'bottom'), hE);
                        hold on
                    end                    
                end

                if f == 2 && r == 1
                    yUnit = 1;
                    yMax = 5.0;
                elseif f == 1 && r == 1
                    yUnit = 1;
                    yMax = 3;
                else
                    yUnit = 0.5;
                    yMax = 2;
                end

                %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
                %yMax = ceil(max(valSet(:)))+yUnit;
                set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0.1,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top','clipping','off');
                xlim([xMin,xMax]);
                ylim([0,yMax])
                % plot noise patch
                meanNoise = max(noiseVals(:,curFreq,r,curConds),[],4)';
                yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
                xNoiseVals = [xMin,xMin,binVals,xMax,xMax];
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
            export_fig(sprintf('%s/exp%d_rc%d_%s_snr.pdf',saveLocation,doExp,f,rcaType),'-pdf','-transparent',gcf);
        else
            export_fig(sprintf('%s/exp%d_rc%d_%s.pdf',saveLocation,doExp,f,rcaType),'-pdf','-transparent',gcf);
        end
        if f < 3
            for v = 1:2
                for z=1:5;
                    fH = figure;
                    set(fH, 'units', 'centimeters');
                    figPos = get(fH,'pos');
                    figPos(4) = figHeight/6;
                    figPos(3) = figWidth*(1/3);
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
                    ylabel('Amplitude (\muV)')
                    %xlabel('Distance (arcmins)');
                    if ~exist(sprintf('%s/exp%d_talkFigures',saveLocation,doExp),'dir')
                        mkdir(sprintf('%s/exp%d_talkFigures',saveLocation,doExp));
                    else
                    end
                    if v == 1
                        export_fig(sprintf('%s/exp%d_talkFigures/exp%d_rc%d_%s_hori%d.pdf',saveLocation,doExp,doExp,f,rcaType,z),'-pdf','-transparent',fH);
                    else
                        export_fig(sprintf('%s/exp%d_talkFigures/exp%d_rc%d_%s_vert%d.pdf',saveLocation,doExp,doExp,f,rcaType,z),'-pdf','-transparent',fH);
                    end
                    close(fH);
                end
            end
        else
        end
    end    
    close all;
end
