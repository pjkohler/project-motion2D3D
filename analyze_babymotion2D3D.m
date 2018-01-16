%% Add Paths
close all;
 codeFolder = '/Users/kohler/code';
 rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
 addpath(genpath(rcaCodePath));
 addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
 addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
 setenv('DYLD_LIBRARY_PATH','')

% codeFolder = '/Users/labmanager/Desktop/LabManager/MatAnal';
% rcaCodePath = sprintf('%s/rcaBase',codeFolder);
% addpath(genpath(codeFolder));
% addpath(genpath(rcaCodePath));
% addpath(genpath(sprintf('%s/mrC',codeFolder)));
% addpath(genpath(sprintf('%s/matlab_lib',codeFolder)));
% addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
% addpath(genpath(sprintf('%s/export_fig',codeFolder)));
% setenv('DYLD_LIBRARY_PATH','')

%mainPath = '/Users/labmanager/Desktop/LabManager/WM_Data/2017_2DBaby';
mainPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
figureFolder = sprintf('%s/figures/exp3',mainPath);

load(sprintf('%s/BabyDataOutput_2F1',figureFolder)); %file name from prep workspace export
%% PLOT RCs
close all;
nFreq = length(allRCA.settings.freqsToUse);
plotSNR = false;
rcaType = 'freq';
figHeight = 5; % set figure size in the beginning
figWidth = 21;
binVals = cellfun(@(x) str2num(x), allRCA.settings.binLevels{1});
lWidth = 1.5;
errorWidth = .3;
cBrewer = load('colorBrewer.mat');
color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
subColors = repmat([color1; color2],2,1);
subColors = subColors(condsToUse,:);
fSize = 12;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
condLabels = condLabels(condsToUse);
xMin = 1.8; xMax = 33.5;

if freqsToUse == 2
    nComp = 1;   %change to reflect desired plotting of only 1st RC.
end

for f= 1:length(freqsToUse)
    if strcmp(rcaType,'all');
        curFreq = f+1;
    else
        curFreq = f;
    end
    figure;
    for r = 1:length(nComp)
        
        if r<nComp+1
            egiH(r) = subplot(length(nComp),3,3+(r-1)*2);
            hold on
            rcaColorBar = [min(freqRCA(curFreq).A(:,nComp(r))),max(freqRCA(curFreq).A(:,nComp(r)))];
            newExtreme = max(abs(rcaColorBar));
            rcaColorBar = [-newExtreme,newExtreme];
            [figH(f),cH(f)] = mrC.plotOnEgi(freqRCA(curFreq).A(:,nComp(r)),rcaColorBar,true);
            hold off
        else
        end
        
        set(cH(f),'fontsize',fSize,'fontname','Helvetica','YTick',linspace(round(min(rcaColorBar(:,f)),1),round(min(rcaColorBar(:,f))*-1,1),5));
        xlabel(cH(f),'weights','fontsize',fSize,'fontname','Helvetica')
        set(cH(f),'location','southoutside');
        set(cH(f),'units','centimeters');
        cBarPos = [14.5,1.2,3,.35];
        set(cH(f),'position',cBarPos);
        
        %title(['Full RCA ' num2str(r)],'fontsize',fSize,'fontname','Arial');
        sH(f,r,1) = subplot(length(nComp),3,2+(r-1)*2);
        curConds = find(ismember(condsToUse,[1:4]));
        titleStr = sprintf('horizontal: %s',allRCA.settings.freqLabels{1}{f});
        hold on
        
        for c=1:length(curConds)
            valSet = ampVals(:,curFreq,nComp(r),:);
            ampH(c) =plot(binVals,ampVals(:,curFreq,nComp(r),curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth*1.5,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
            hE = ErrorBars(binVals,ampVals(:,curFreq,nComp(r),curConds(c)),[errLB(:,curFreq,nComp(r),curConds(c)),errUB(:,curFreq,nComp(r),curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
            uistack(ampH(c),'bottom')
            cellfun(@(x) uistack(x,'bottom'), hE);
            hold on
            
            if ~isnan(NR_pOpt(1,curFreq,nComp(r),curConds(c))) % plot Naka-Rushton
                nFine = 1e2;
                nrX = linspace( min(binVals), max(binVals), nFine )';
                nrVals = hModel( nrX, NR_pOpt(:,curFreq,nComp(r),curConds(c)));
                nR = plot( nrX, nrVals, '-k','LineWidth',lWidth);
                uistack(nR,'bottom')
            else
            end
        end
        
        if freqsToUse(f) == 2 && r == 1
            yUnit = 3;
            yMax = 12.0;
        else
            yUnit = 1;
            yMax = 4;
        end
        %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
        %yMax = ceil(max(valSet(:)))+yUnit;
        get(0,'Factory');
        set(0,'defaultfigurecolor',[1 1 1]);
        set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[2,4,8,16,32],'ytick',0:yUnit:yMax,'Layer','top','clipping','off', 'color','none');
        set(gca,'XMinorTick','off');
        xlim([xMin,xMax]);
        ylim([0,yMax])
        
        % plot noise patch
        meanNoise = max(noiseVals(:,curFreq,nComp(r),curConds),[],4)';
        yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
        xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
        pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
        uistack(pH,'bottom')
        
        %             if s==1
        %                 if r > nComp
        %                     text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,'OZ','fontsize',fSize,'fontname','Arial');
        %                 else
        %                     text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,['RC ',num2str(r)],'fontsize',fSize,'fontname','Arial');
        %                 end
        %             else
        %             end
        if r==length(nComp)
            title(titleStr,'fontsize',fSize,'fontname','Arial');
            ylabel('amplitude (\muV)');
            xlabel('distance (arcmins)');
            %else
            %                 lH = legend(ampH,condLabels(curConds),'location','northeast');
            %                 legend boxoff
            %                 lPos = get(lH,'position');
            %                 lPos(1) = lPos(1) + .13;
            %                 lPos(2) = lPos(2) + -.04; %.05
            %                 set(lH,'position',lPos);
            if freqsToUse ==1
                tH = text(45,4,sprintf('n = %0d',size(allRCA.data,2)),'fontsize',fSize,'fontname','Arial'); %(60,11,...) for first 2 args in 2F1 plot; (60,4...) for 1F1
            else
                tH = text(45,11,sprintf('n = %0d',size(allRCA.data,2)),'fontsize',fSize,'fontname','Arial');
            end
        end
        hold off;
        drawnow;
        
        addX = 1.6;
        addY = 1.6;
        newPos = get(egiH(r),'position');
        newPos(1) = newPos(1)-(newPos(3)*addX*.5);
        %newPos(2) = newPos(2)-(newPos(4)*addY*.4);
        newPos(2) = -0.2;
        newPos(3) = newPos(3)*(.5+addX);
        newPos(4) = newPos(4)*(.5+addY);
        set(egiH(r),'position',newPos);
        
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = figHeight;
        figPos(3) = figWidth;
        set(gcf,'pos',figPos);
        
        export_fig(sprintf('%s/2DInfant_%dF1',figureFolder,freqsToUse),'-painters','-opengl','-png','-eps','-pdf','-transparent',gcf);
    end
end
save(sprintf('%s/BabyDataOutput_%dF1',figureFolder,freqsToUse));