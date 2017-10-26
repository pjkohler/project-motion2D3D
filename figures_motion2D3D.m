%prep_motion2D3D;
%analyze_motion2D3D(1);
%analyze_motion2D3D(2);
%analyze_motion2D3D(4);
%analyze_motion2D3D(5);
%analyze_motion2D3D(8);

adultExp = [1,2,4,5,8];
plotType = 'freq';
figFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures';
%figFolder = '/Users/kohler/Desktop/figures';
%% PLOT RCs
% set figure size in the beginning
figHeight = 25;
figWidth = 20;
condsToUse = 1:8;
% figure params
lWidth = 2;
fSize = 12;
cBrewer = load('colorBrewer.mat');
color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
subColors = repmat([color1; color2],2,1);
subColors = subColors(condsToUse,:);
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
condLabels = condLabels(condsToUse);
xMin = .4; xMax = 18;

nFreq = 4;
numExp = length(adultExp);
close all;
rcNum = 1;
figLabel = {'A','B','C','D','E'};
for e = 1:numExp
    load(sprintf('%s/exp%0.0f/plottingData.mat',figFolder,adultExp(e)));
    for f = 1:nFreq
        figure(f);
        if strcmp(plotType,'all'); 
            curFreq = f+8;
            % grab frequency name from condition 1, 
            % will be the same for all conditions
            freqName = readyRCA(curFreq).settings.freqLabels{1}{f};
        elseif strcmp(plotType,'full'); 
            curFreq = f+4;
            freqName = readyRCA(curFreq).settings.freqLabels{1}{f};
        elseif strcmp(plotType,'freq'); 
            curFreq = f;
            freqName = readyRCA(curFreq).settings.freqLabels{f}{1};
        else
            msg = sprintf('\n unknown plot type: %s\n',plotType);
            error(msg);
        end
        % shared params
        if e == 1 && f == 1
            binVals = cellfun(@(x) str2num(x), readyRCA(curFreq).settings.binLevels{1});
        else
        end
        % plot topography
        egiH(e,f) = subplot(numExp,3,3+(e-1)*3);
        hold on
        rcaColorBar = [min(readyRCA(curFreq).A(:,rcNum)),max(readyRCA(curFreq).A(:,rcNum))];
        newExtreme = max(abs(rcaColorBar));
        rcaColorBar = [-newExtreme,newExtreme];
        mrC.plotOnEgi(readyRCA(curFreq).A(:,rcNum),rcaColorBar);
        hold off
            
        %title(['Full RCA ' num2str(rcNum)],'fontsize',fSize,'fontname','Helvetica');
        for s=1:2 % vertical or horizontal motion
            if s==1
                sH(f,e,1) = subplot(numExp,3,1+(e-1)*3);
                curConds = find(ismember(condsToUse,1:4));
                titleStr = sprintf('horizontal: %s',freqName);
            else
                sH(f,e,2) = subplot(numExp,3,2+(e-1)*3);
                curConds = find(ismember(condsToUse,5:8));
                titleStr = sprintf('vertical: %s',freqName);
            end
            hold on

            valSet = readyRCA(curFreq).stats.ampVals;
            errSet1 = readyRCA(curFreq).stats.errLB;
            errSet2 = readyRCA(curFreq).stats.errUB;
            NRset = readyRCA(curFreq).stats.NR_pOpt;
            NRmodel = readyRCA(curFreq).stats.hModel;
            for c=1:length(curConds)
                valH(c)=plot(binVals,valSet(:,rcNum,curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
                hE = ErrorBars(binVals,valSet(:,rcNum,curConds(c)),[errSet1(:,rcNum,curConds(c)),errSet1(:,rcNum,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
                uistack(valH(c),'bottom')
                cellfun(@(x) uistack(x,'bottom'), hE);
                hold on
                if ~isnan(NRset(1,rcNum,curConds(c)))
                    % plot Naka-Rushton
                    nFine = 1e2;
                    nrX = linspace( min(binVals), max(binVals), nFine )';
                    nrVals = NRmodel( nrX, NRset(:,rcNum,curConds(c)));
                    hNR{c} = plot( nrX, nrVals, '-k', 'LineWidth',lWidth);
                else
                end
            end
            cellfun(@(x) uistack(x,'bottom'), hNR);

            if f == 2
                yUnit = 1;
                yMax = 5.0;
            elseif f == 1
                yUnit = 0.5;
                yMax = 3.5;
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
            meanNoise = max(readyRCA(curFreq).stats.noiseVals(:,rcNum,curConds),[],3)';
            yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
            xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
            pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
            uistack(pH,'bottom')

            if s==1
                text(0.04,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,sprintf('%s: Exp. %0.0f',figLabel{e},adultExp(e)),'fontsize',fSize*1.5,'fontname','Helvetica');

                %if r > nComp
                %    text(0.0008,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,'OZ','fontsize',fSize,'fontname','Helvetica');
                %else
                %    text(0.0008,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,['RC ',num2str(r)],'fontsize',fSize,'fontname','Helvetica');
                %end
            else
                text(16,max(get(gca,'ylim')),sprintf('n = %0d',size(readyRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Helvetica');
            end
            if e==1
                title(titleStr,'fontsize',fSize,'fontname','Helvetica');
            elseif e==numExp
                if s==1
                    ylabel('amplitude (\muV)')
                    xlabel('distance (arcmins)');
                else
                    lH = legend(valH,condLabels(curConds),'location','northeast','orientation','horizontal');
                    legend boxoff
                    set(lH,'PlotBoxAspectRatio',[12 1 1]);
                    lPos = get(lH,'position');
                    lPos(1) = lPos(1) + .28;
                    lPos(2) = lPos(2) - .13;
                    set(lH,'position',lPos);
                    ch = get(lH, 'Children');
                    textCh = ch(strcmp(get(ch, 'Type'), 'text'));
                    for iText = 1:numel(textCh)
                        set(textCh(iText), 'Position', get(textCh(iText), 'Position') + [-0.03 0 0])
                    end
                end
            end
            hold off;
        end
    end
    drawnow;
end
for f = 1:nFreq
    figure(f)
    for e = 1:numExp
        addX = 1.4;
        addY = 1.4;
        newPos = get(egiH(e,f),'position');
        newPos(1) = newPos(1)-(newPos(3)*addX*.7);
        newPos(2) = newPos(2)-(newPos(4)*addY*.6);
        newPos(3) = newPos(3)*(1+addX);
        newPos(4) = newPos(4)*(1+addY);
        set(egiH(e,f),'position',newPos);
    end
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = figHeight;
    figPos(3) = figWidth;
    set(gcf,'pos',figPos);
    
    export_fig(sprintf('%s/paper_figures/adultExp_harm%0.0f.pdf',figFolder,f),'-pdf','-transparent',gcf);
    
%     if plotSNR        
%         
%     else
%         export_fig(sprintf('%s/exp%d_rc%d_%s.pdf',saveLocation,doExp,f,plotType),'-pdf','-transparent',gcf);
%     end
end
       
        
    
        


    
    
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

    

    

    