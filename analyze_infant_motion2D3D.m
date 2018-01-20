function analyze_infant_motion2D3D(doFreq)
    %% Add Paths
    close all;
    
    codeFolder = '/Users/kohler/code';
    if exist(codeFolder,'dir')
        rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
        addpath(genpath(rcaCodePath));
        addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
        addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
    else
    end
    setenv('DYLD_LIBRARY_PATH','')

    if nargin < 1
        doFreq = [1,2];
    else
    end

    if numel(doFreq) > 1
        for z = 1:length(doFreq)
            analyze_infant_motion2D3D(doFreq(z));
        end
        return;
    else
    end
    
    if doFreq == 2
        nComp = 1;   % plot only 1st RC.
    elseif doFreq == 1
        nComp = 5;   % plot only 5th RC.
    else
        error('freq %0.0f, not sure which RC to plot!',doFreq);
    end

    mainPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
    figureFolder = sprintf('%s/figures/infant_exp',mainPath);
    load(sprintf('%s/BabyDataOutput_%dF1', figureFolder, doFreq),'babyRCA'); %file name from prep workspace export
    rcaType = 'freq';

    %% PLOT RCs

    if diff(arrayfun(@(x) babyRCA(x).settings.freqsToUse, 1:length(babyRCA))) ~= 0
        error('different frequencies used in different RCA iterations');
    else
        freqsToUse = babyRCA(1).settings.freqsToUse;
        nFreq = length(freqsToUse);
    end

    if any( diff(cell2mat(arrayfun(@(x) babyRCA(x).settings.condsToUse, 1:length(babyRCA),'uni',false)')) ~= 0 )
        error('different conditions used in different RCA iterations');
    else
        condsToUse = babyRCA(1).settings.condsToUse;
    end

    binVals = cellfun(@(x) str2num(x), babyRCA(1).settings.binLevels{1});

    figHeight = 5; % set figure size in the beginning
    figWidth = 21;
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
    
    for f= 1:length(freqsToUse)
        if strcmp(rcaType,'all');
            curFreq = length(freqsToUse)+f;
        else
            curFreq = f;
        end
        
        figure;
        for r = 1:length(nComp)

            if r<nComp+1
                egiH(r) = subplot(length(nComp),3,3+(r-1)*2);
                hold on
                rcaColorBar(:,f) = [min(babyRCA(curFreq).A(:,nComp(r))),max(babyRCA(curFreq).A(:,nComp(r)))];
                newExtreme = round(max(abs(rcaColorBar(:,f)))*10)./10;
                rcaColorBar(:,f) = [-newExtreme,newExtreme*1.001];
                [figH(f),cH(f)] = mrC.plotOnEgi(babyRCA(curFreq).A(:,nComp(r)),rcaColorBar(:,f),true);
                hold off
            else
            end
            set(cH(f),'YTickMode','manual','location','eastoutside','units','centimeters')
            set(cH(f),'fontsize',fSize,'fontname','Helvetica','YTickMode','manual','Ytick',linspace(min(rcaColorBar(:,f)),min(rcaColorBar(:,f))*-1,5));
            ylabel(cH(f),'weights','fontsize',fSize,'fontname','Helvetica')
            cBarPos = [17.1,1.4,.3,3];
            set(cH(f),'position',cBarPos);

            %title(['Full RCA ' num2str(r)],'fontsize',fSize,'fontname','Helvetica');
            sH(f,r,1) = subplot(length(nComp),3,2+(r-1)*2);
            curConds = find(ismember(condsToUse,[1:4]));
            titleStr = sprintf('horizontal: %s',babyRCA(1).settings.freqLabels{1}{f});
            hold on

            valSet = babyRCA(curFreq).stats.Amp;
            errSet1 = babyRCA(curFreq).stats.ErrLB;
            errSet2 = babyRCA(curFreq).stats.ErrUB;
            NRvals = babyRCA(curFreq).stats.NR_Params;
            NRerrs = babyRCA(curFreq).stats.NR_JKSE;
            NRmodel = babyRCA(curFreq).stats.hModel;

            % do paired tests
            testVal = squeeze(babyRCA(curFreq).stats.NR_JKParams(:,:,nComp(r),curConds));
            testVal = permute(testVal,[2,1,3]); % move subjects to first dim
            tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp
            paramIdx = [1,2,3,4];% only look at c50 and rMax
            jkDf = size(testVal,1)-1;
            for pT=1:length(tIdx) % do four paired tests
                diffErr = jackKnifeErr(testVal(:,paramIdx,tIdx(pT,1))-testVal(:,paramIdx,tIdx(pT,2)));
                grandDiff = NRvals(paramIdx,nComp(r),tIdx(pT,1)) - NRvals(paramIdx,nComp(r),tIdx(pT,2));
                paramPairedT(:,pT,f) = grandDiff'./diffErr;
                paramPairedP(:,pT,f) = 2*tcdf( -abs(paramPairedT(:,pT,f)) , jkDf);
            end

            for c=1:length(curConds)
                ampH(c) =plot(binVals,valSet(:,nComp(r),curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
                hE = ErrorBars(binVals,valSet(:,nComp(r),curConds(c)),[errSet1(:,nComp(r),curConds(c)),errSet1(:,nComp(r),curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
                uistack(ampH(c),'bottom')
                cellfun(@(x) uistack(x,'bottom'), hE);
                hold on

                if ~isnan(NRvals(1,nComp(r),curConds(c))) % plot Naka-Rushton
                    nFine = 1e2;
                    nrX = linspace( min(binVals), max(binVals), nFine )';
                    nrVals = NRmodel( nrX, NRvals(:,nComp(r),curConds(c)));
                    nR(c) = plot( nrX, nrVals, '-','color',subColors(curConds(c),:),'LineWidth',lWidth);

                else
                end
            end
            arrayfun(@(x) uistack(x,'bottom'),nR);

            if freqsToUse(f) == 2 && r == 1
                yUnit = 2;
                yMax = 10.0;
            else
                yUnit = 1;
                yMax = 4;
            end
            %         for z = 1:length(yLabels)
            %             if strcmp(yLabels{z}(1),'0')
            %                 yLabels{z}(1) = ' ';
            %             else
            %             end
            %         end

            %yUnit = floor((ceil(max(valSet(:)))/5)*10)/10;
            %yMax = ceil(max(valSet(:)))+yUnit;
            set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[2,4,8,16,32],'ytick',0:yUnit:yMax,'Layer','top','clipping','off', 'color','none');
            set(gca,'XMinorTick','off');
            xlim([xMin,xMax]);
            ylim([0,yMax])

            % plot noise patch
            meanNoise = max(babyRCA(curFreq).stats.NoiseAmp(:,nComp(r),curConds),[],3)';
            yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
            xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
            pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
            uistack(pH,'bottom')

            %             if s==1
            %                 if r > nComp
            %                     text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,'OZ','fontsize',fSize,'fontname','Helvetica');
            %                 else
            %                     text(0.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.1,['RC ',num2str(r)],'fontsize',fSize,'fontname','Helvetica');
            %                 end
            %             else
            %             end
            if r==length(nComp)
                text(2,yMax*1.1,titleStr,'fontsize',fSize,'fontname','Helvetica');
                ylabel('amplitude (\muV)');
                xlabel('distance (arcmins)');
                %else
                %                 lH = legend(ampH,condLabels(curConds),'location','northeast');
                %                 legend boxoff
                %                 lPos = get(lH,'position');
                %                 lPos(1) = lPos(1) + .13;
                %                 lPos(2) = lPos(2) + -.04; %.05
                %                 set(lH,'position',lPos);
                tH = text(35,yMax*1.1,sprintf('n = %0d',size(babyRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Helvetica'); %(60,11,...) for first 2 args in 2F1 plot; (60,4...) for 1F1
            end
            hold off;
            drawnow;
            addX = 0.25;
            addY = 0.25;
            set(egiH(r),'units','centimeters');
            if doFreq == 1
                newPos = get(egiH(r),'position');
                newPos(1) = newPos(1) - newPos(3)*addX*1.5;
                newPos(2) = newPos(2) - newPos(4)*addY*2.2;
                newPos(3) = newPos(3)+newPos(3)*addX;
                newPos(4) = newPos(4)+newPos(4)*addY;
                set(egiH(r),'position',newPos);
                save(sprintf('%s/EEGpos.mat',figureFolder),'newPos');
            else
                load(sprintf('%s/EEGpos.mat',figureFolder),'newPos');
                set(egiH(r),'position',newPos);
            end

            set(gcf, 'units', 'centimeters');
            figPos = get(gcf,'pos');
            figPos(4) = figHeight;
            figPos(3) = figWidth;
            set(gcf,'pos',figPos);

            %export_fig(sprintf('%s/2DInfant_%dF1',figureFolder,freqsToUse),'-painters','-opengl','-png','-eps','-pdf','-transparent',gcf);
            export_fig(sprintf('%s/2DInfant_%dF1',figureFolder,freqsToUse),'-pdf','-transparent',gcf);
        end
    end
    save(sprintf('%s/BabyDataOutput_%dF1',figureFolder,freqsToUse));
    %% PLOT NR
    close all;
    figure;
    for z = 1:length(paramIdx)
        subplot(1,4,z);
        hold on
        plot(1:4,squeeze(NRvals(paramIdx(z),nComp(r),:)),'-o')
        errorb(1:4,squeeze(NRvals(paramIdx(z),nComp(r),:)),squeeze(NRerrs(paramIdx(z),nComp(r),:)));
        hold off
        xlim([0.5,4.5]);
    end
end
