%prep_motion2D3D;
%analyze_motion2D3D(1);
%analyze_motion2D3D(2);
%analyze_motion2D3D(4);
%analyze_motion2D3D(5);
%analyze_motion2D3D(8);
clear all;
close all;

adultExp = [1,2,3,4,5];
plotType = 'freq';
figFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures';
%figFolder = '/Users/kohler/Desktop/figures';
%% PLOT RCs
% set figure size in the beginning
figHeight = 25;
figWidth = 17.8;
condsToUse = 1:8;
% figure params
lWidth = 1.5;
fSize = 12;
cBrewer = load('colorBrewer.mat');
color1 = [cBrewer.rgb20(3,:); cBrewer.rgb20(4,:)];
color2 = [cBrewer.rgb20(5,:); cBrewer.rgb20(6,:)];
mainColors = repmat([color1; color2],2,1);
mainColors = mainColors(condsToUse,:);
alternaColors = mainColors(1:4,:);
alternaColors(2,:) = cBrewer.rgb20(9,:);
alternaColors(4,:) = cBrewer.rgb20(19,:);
alternaColors = repmat(alternaColors,2,1);
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};
condLabels = repmat({'rel-Mot','abs-Mot','rel-Disp','abs-Disp'},1,2);
condLabels = condLabels(condsToUse);
xMin = .4; xMax = 18;

nFreq = 4;
numExp = length(adultExp);
rcNum = 1;
figLabel = {'A','B','C','D','E','F','G','H','I','J'};
flipIdx = [1,-1,1,1,1;
           -1,-1,1,1,1;
           -1,-1,1,1,-1;
           -1,-1,-1,-1,-1];
for e = 1:numExp
    load(sprintf('%s/exp%0.0f/plottingData.mat',figFolder,adultExp(e)));
    if any(e == [3,4])
        subColors = alternaColors;
    else
        subColors = mainColors;
    end
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
        if e == 1
            if f == 1
                clear rcaColorBar;
            else
            end
            rcaColorBar(:,f) = [min(readyRCA(curFreq).A(:,rcNum)),max(readyRCA(curFreq).A(:,rcNum))];
            newExtreme = round(max(abs(rcaColorBar(:,f)))*10)./10;
            rcaColorBar(:,f) = [-newExtreme,newExtreme*1.001];
        else
        end
        [figH(e,f),cH(e,f)] = mrC.plotOnEgi(readyRCA(curFreq).A(:,rcNum)*flipIdx(f,e),rcaColorBar(:,f),true);
        
        hold off
            
        %title(['Full RCA ' num2str(rcNum)],'fontsize',fSize,'fontname','Helvetica');
        for s=1:2 % vertical or horizontal motion
            if s==1
                sH(f,e,1) = subplot(numExp,3,1+(e-1)*3);
                curConds = find(ismember(condsToUse,1:4));
                titleStr = ['horizontal: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            else
                sH(f,e,2) = subplot(numExp,3,2+(e-1)*3);
                curConds = find(ismember(condsToUse,5:8));
                titleStr = ['vertical: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            end
            hold on

            valSet = readyRCA(curFreq).stats.Amp;
            errSet1 = readyRCA(curFreq).stats.ErrLB;
            errSet2 = readyRCA(curFreq).stats.ErrUB;
            NRset = readyRCA(curFreq).stats.NR_Params;
            NRmodel = readyRCA(curFreq).stats.hModel;
            
            % store NR for later
            NakaRushton(e,f).Params = readyRCA(curFreq).stats.NR_Params;
            NakaRushton(e,f).R2 = readyRCA(curFreq).stats.NR_R2;
            NakaRushton(e,f).JKSE = readyRCA(curFreq).stats.NR_JKSE;
            NakaRushton(e,f).JKParams = readyRCA(curFreq).stats.NR_JKParams;
            
            % store t-values for later
            rc_tSqrdP(e,f,s) = {permute(squeeze(readyRCA(curFreq).stats.tSqrdP(:,rcNum,curConds)),[2,1])};
            rc_tSqrdVal(e,f,s) = {permute(squeeze(readyRCA(curFreq).stats.tSqrdVal(:,rcNum,curConds)),[2,1])};
            % make new p-values
            [rcaDataReal,rcaDataImag] = getRealImag(readyRCA(curFreq).data(curConds,:));
            rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataReal,'uni',false);
            rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
            rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataImag,'uni',false);
            rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
            
            tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp
            
            for pT=1:length(tIdx) % do four paired tests
                for b = 1:length(binVals)
                    xyData = permute(cat(1,rcaDataReal(b,:,tIdx(pT,:)),rcaDataImag(b,:,tIdx(pT,:))),[2,1,3]);
                    tempStrct = tSquaredFourierCoefs(xyData);
                    tempP(:,b) = tempStrct.pVal;
                    tempT(:,b) = tempStrct.tSqrd;
                end
                rc_tSqrdP(e,f,s) = {cat(1,rc_tSqrdP{e,f,s},tempP)};
                rc_tSqrdVal(e,f,s) = {cat(1,rc_tSqrdVal{e,f,s},tempT)};
                clear temp*
            end
                
                    
            
            
            %tSquaredFourierCoefs(xyData(~nanVals,:));
            
            for c=1:length(curConds)
                valH(c)=plot(binVals,valSet(:,rcNum,curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
                hE = ErrorBars(binVals,valSet(:,rcNum,curConds(c)),[errSet1(:,rcNum,curConds(c)),errSet2(:,rcNum,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
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
                yUnit = 1;
                yMax = 4;
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
            meanNoise = max(readyRCA(curFreq).stats.NoiseAmp(:,rcNum,curConds),[],3)';
            yNoiseVals = [0,meanNoise(1),meanNoise,meanNoise(end),0]; % start and end points just repeats of first and last
            xNoiseVals = [xMin,xMin,binVals',xMax,xMax];
            pH = patch(xNoiseVals,yNoiseVals,[.75 .75 .75],'edgecolor','none');
            uistack(pH,'bottom')
            
            text(0.2,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,figLabel{e+(s-1)*5},'fontsize',fSize*1.5,'fontname','Helvetica');
            text(0.45,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.05,sprintf('Exp. %0.0f',adultExp(e)),'fontsize',fSize,'fontname','Helvetica');

            if s==2
                text(16,max(get(gca,'ylim')),sprintf('n = %0d',size(readyRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Helvetica');
            else
            end
            if e==1
                title(titleStr,'fontsize',fSize,'fontname','Helvetica','interpreter','tex');
            elseif e==numExp
                if s==1
                    ylabel('amplitude (\muV)')
                    xlabel('displacement (arcmins)');
                else
                    %lH = legend(valH,condLabels(curConds),'location','northeast','orientation','horizontal');
                    %legend boxoff
                    %set(lH,'PlotBoxAspectRatio',[12 1 1]);
                    %lPos = get(lH,'position');
                    %lPos(1) = lPos(1) + .28;
                    %lPos(2) = lPos(2) - .13;
                    %set(lH,'position',lPos);
                    %ch = get(lH, 'Children');
                    %textCh = ch(strcmp(get(ch, 'Type'), 'text'));
                    %for iText = 1:numel(textCh)
                    %    set(textCh(iText), 'Position', get(textCh(iText), 'Position') + [-0.03 0 0])
                    %end
                end
            end
            hold off;
        end
    end
    drawnow;
end

%% ADJUST TOPO PLOTS
for f = 1:nFreq
    figure(f)
    addX = 1.8;
    addY = 1.8;
    for e = 1:numExp
         if e == numExp
            set(cH(e,f),'fontsize',fSize,'fontname','Helvetica','YTick',linspace(min(rcaColorBar(:,f)),min(rcaColorBar(:,f))*-1,5));
            xlabel(cH(e,f),'weights','fontsize',fSize,'fontname','Helvetica')
            set(cH(e,f),'location','southoutside');
            set(cH(e,f),'units','centimeters');
            %cBarPos = get(cH(e,f),'position');
            cBarPos = [10,1.4,4,.4];
            set(cH(e,f),'position',cBarPos);
            oldPos = newPos;
            newPos = get(egiH(e,f),'position');
            newPos(1) = oldPos(1);
            newPos(2) = newPos(2)-(newPos(4)*addY*.6);
            newPos(3) = oldPos(3);
            newPos(4) = oldPos(4);
        else
            set(cH(e,f),'visible','off');
            newPos = get(egiH(e,f),'position');
            newPos(1) = newPos(1)-(newPos(3)*addX*.7);
            newPos(2) = newPos(2)-(newPos(4)*addY*.6);
            newPos(3) = newPos(3)*(1+addX);
            newPos(4) = newPos(4)*(1+addY);
        end
        set(egiH(e,f),'position',newPos); 
    end
    drawnow;
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

%% MAKE STATS FIGURES
close all;
fSize = 10;
lWidth = 1.5;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

figHeight = 25;
figWidth = 17.8*(2/3);
numExp = 5;
numTests = 8; % 4 conds, in-phase vs anti-phase, absolute vs relative
xVals = 1:10;
testSymbols = {'o','o','o','o','o','o','o','o'};
markerOpts = {'MarkerSize',6,'LineWidth',lWidth,'markerfacecolor',[1,1,1]};
freqLabels = {'1F','2F','3F','4F'};
% significance colors
tHotColMap = jmaColors('pval');
tHotColMap(end,:) = [1 1 1];
pCutOff = 0.05;
for e = 1:numExp
    if any(e == [3,4])
        subColors = alternaColors;
    else
        subColors = mainColors;
    end
    markerColors{1} = [subColors(1:4,:);subColors([2,4,3,4],:)];
    markerColors{2} = subColors([1,3,1,2],:);
    for f = 1:nFreq
        figure(f);
        freqName = freqLabels{f};
        for s=1:2 % vertical or horizontal motion
            if s==1
                sH(f,e,1) =  subplot(numExp,2,s+(e-1)*2);
                curConds = find(ismember(condsToUse,1:4));
                titleStr = ['horizontal: \it\fontname{Arial}',sprintf('%s',freqName)];
            else
                sH(f,e,2) =  subplot(numExp,2,s+(e-1)*2);
                curConds = find(ismember(condsToUse,5:8));
                titleStr = ['vertical: \it\fontname{Arial}',sprintf('%s',freqName)];
            end
            hold on
           
            curP = rc_tSqrdP{e,f,s};
            
            %hP = plot([0,20],ones(1,2)*.5,'k-','linewidth',lWidth);
            rectangle('position',[-1.5,4.5,2,4],'FaceColor',[.75 .75 .75],'EdgeColor','none');
            rectangle('position',[-1.5,0.5,2,4],'FaceColor',[.5 .5 .5],'EdgeColor','none');
            hP = arrayfun(@(x) plot([-2,20],ones(1,2)*x,'k-','linewidth',lWidth),0.5:1:numTests);
            hP = [hP,arrayfun(@(x) plot(ones(1,2)*x,[0,10],'k-','linewidth',lWidth),0.5)];

            mP = arrayfun(@(x) plot(0,x,sprintf('k%s',testSymbols{x}),'color',markerColors{1}(x,:),markerOpts{:}),1:length(testSymbols));
            mP = [mP,arrayfun(@(x) plot(-1,x+numTests/2,sprintf('k%s',testSymbols{x}),'color',markerColors{2}(x,:),markerOpts{:}),1:(numTests/2))];
            hImg = image([1,10],[1,numTests],curP, 'CDataMapping', 'scaled','Parent',gca);
            colormap(gca,tHotColMap );   
            cMapMax = pCutOff+2*pCutOff/(size(tHotColMap,1));
            set(gca,'xtick',2:2:10,'xticklabel',arrayfun(@(x) sprintf('%0.1f',binVals(x)),2:2:10,'uni',false),...
                'ytick',[],'ycolor',[1,1,1],'yaxislocation','left','yticklabel','','CLim', [ 0 cMapMax ], gcaOpts{:} ); % set range for color scale
            xlim([-1.5,10.5]);
            ylim([.5,numTests+.5]);
            uistack(hImg,'bottom');
            uistack(hP,'top');
            %uistack(mP,'top');
            text(-2,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,figLabel{e+(s-1)*5},'fontsize',fSize*1.5,'fontname','Helvetica');
            if s == 1
                ylabel(sprintf('Exp. %0.0f',adultExp(e)),'color','k','fontsize',fSize,'fontname','Helvetica');
                if e==numExp
                    xlabel('displacement (arcmins)');
                else
                end
            else
            end
            if e ==1                 
                title(titleStr,'fontsize',fSize,'fontname','Helvetica');
            else
            end
            
            hold off;
        end
    end
end

for f = 1:nFreq
    
    figure(f)
    drawnow;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = figHeight;
    figPos(3) = figWidth;
    set(gcf,'pos',figPos);
%     for e = 1:numExp
%         subPos = get(sH(f,e,1), 'Position'); 
%         subPosNew = get(sH(f,e,2), 'Position');
%         %subPosNew(1) = subPosNew(1) - 0.04;
%         %subPosNew(2) = subPosNew(2) - 0.02;
%         %subPosNew(4) = subPosNew(4) + 0.04;
%         %subPosNew(2) = subPos(2)-0.02;
%         set(sH(f,e,2), 'Position', subPosNew);
%     end
        export_fig(sprintf('%s/paper_figures/adultExp_stats%0.0f.pdf',figFolder,f),'-pdf','-transparent',gcf);

end


%% MAKE NAKA-RUSHTON FIGURES
close all;
for e = 1:numExp
    if any(e == [4,5])
        subColors = alternaolors;
    else
        subColors = mainColors;
    end
    for f = 1:nFreq
        figure(f);
        freqName = freqLabels{f};
        for s=1:2 % vertical or horizontal motion
            if s==1
                curConds = find(ismember(condsToUse,1:4));
                titleStr = ['horizontal: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            else
                curConds = find(ismember(condsToUse,5:8));
                titleStr = ['vertical: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            end
           
            NRvals = cat(1,NakaRushton(e,f).Params,permute(NakaRushton(e,f).R2,[3,1,2]));
            NRvals = squeeze(NRvals(:,rcNum,curConds));
            NRerrors = squeeze(NakaRushton(e,f).JKSE(:,rcNum,curConds));
            
            % do paired tests
            testVal = squeeze(NakaRushton(e,f).JKParams(:,:,rcNum,curConds));
            testVal = permute(testVal,[2,1,3]); % move subjects to first dim
            
            tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp
            jkDf = size(testVal,1)-1;
            for pT=1:length(tIdx) % do four paired tests
                 % only look at c50 and rMax
                diffErr = jackKnifeErr(testVal(:,[1,3,4],tIdx(pT,1))-testVal(:,[1,3,4],tIdx(pT,2)));
                grandDiff = NRvals([1,3,4],tIdx(pT,1)) - NRvals([1,3,4],tIdx(pT,2));
                paramPairedT(:,pT+(s-1)*4,f,e) = grandDiff'./diffErr;
                paramPairedP(:,pT+(s-1)*4,f,e) = 2*tcdf( -abs(paramPairedT(:,pT+(s-1)*4,f,e)) , jkDf);
            end
            
            subPlotOrder = [1,2,5,6];
            paramList = [1,3,4,6];
            for z = 1:4 % five parameters (including R2)
                sH(f,e,s) = subplot(4,1,z);
                hold on
                xVals = (1:4)+5*(s-1)+10*(e-1);
                arrayfun(@(x) bar(xVals(x),NRvals(paramList(z),x),'facecolor',subColors(x,:),'edgecolor','none'),1:4,'uni',false);
                if paramList(z) ~= 6
                    arrayfun(@(x) ...
                        ErrorBars(xVals(x),NRvals(paramList(z),x),NRerrors(paramList(z),x),'color',[0,0,0],'type','bar','cap',false,'barwidth',lWidth),...
                        1:4,'uni',false);
                else
                end
                if e == numExp
                    switch paramList(z)
                        case 1
                            ylabel('c50');
                            yMin = 0; yMax = 15; yUnit = 5;
                        case 2
                            ylabel('exponent');
                            yMin = 0; yMax = 1; yUnit = 1;
                        case 3
                            ylabel('rMax');
                            yMin = 0; yMax = 4; yUnit = 1;
                        case 4
                            ylabel('baseline');
                            yMin = 0; yMax = 2; yUnit = 1;
                        case 6
                            ylabel('R2');
                            yMin = 0; yMax = 1; yUnit = 0.5;
                        otherwise
                    end
                    gcaOpts{4} = [0.01,0.01];
                    ylim([yMin,yMax]);
                    set(gca, gcaOpts{:});               
                    if z == 1
                        set(gca,'YTickLabel',sprintf('%0.0f|',yMin:yUnit:yMax))
                        set(gca,'xtick',[],'ytick',yMin:yUnit:yMax);
                        arrayfun(@(x) text(5+(x-1)*10,max(get(gca,'ylim')), ...
                                [sprintf('Exp. %0.0f: ',adultExp(x)),...
                                '\it\fontname{Arial}',sprintf('%s',freqName(1:2))],...
                                'horizontalalignment','center'),...
                                1:numExp,'uni',false);
                    elseif z == 4
                        set(gca,'YTickLabel',sprintf('%0.1f|',yMin:yUnit:yMax))
                        set(gca,'xtick',[2.5:5:50],'ytick',yMin:yUnit:yMax,'xticklabel',{'Hori','Vert'},'clipping','off');
                    else
                        set(gca,'YTickLabel',sprintf('%0.1f|',yMin:yUnit:yMax))
                        set(gca,'xtick',[],'ytick',yMin:yUnit:yMax);
                    end

                    hold on
                    plot(get(gca,'xlim'),zeros(2,1),'k','color',[0 0 0],'linewidth',lWidth)
                    if z == 4
                        a2 = axes;
                        hold on
                        arrayfun(@(x) plot(ones(1,2)*x,[-1,1],'color',[0 0 0],'linewidth',lWidth),1:4);
                        xlim([0,5]);
                        ylim([-1,1]);
                        set(a2, 'Color', 'none', 'Visible', 'off');
                    else
                    end
                    hold off
                else
                end
            end
        end
    end
end
for f = 1:nFreq
    figure(f)
    drawnow;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = figWidth;
    figPos(3) = figHeight;
    set(gcf,'pos',figPos);
    export_fig(sprintf('/Users/kohler/Desktop/adultExp_Naka%0.0f.pdf',f),'-pdf','-transparent',gcf);        
end

paramNames = {'c50','rMax','b'};
orientNames = {'H','V'};
testNames =  {'InRefVsNo', 'AntiRefVsNo', 'RefInVsAnti', 'NoRefInVsAnti'};

% make table of paired results
freqNum = 2;
for z = 1:4
    for s = 1:2
        curIdx = z+(s-1)*4;
        motion2D3D_stats{curIdx} = table(...
            [squeeze(paramPairedT(1,curIdx,freqNum,:)), squeeze(paramPairedP(1,curIdx,freqNum,:)),squeeze(paramPairedP(1,curIdx,freqNum,:))<0.05], ...
            [squeeze(paramPairedT(2,curIdx,freqNum,:)), squeeze(paramPairedP(2,curIdx,freqNum,:)),squeeze(paramPairedP(2,curIdx,freqNum,:))<0.05], ...
            [squeeze(paramPairedT(3,curIdx,freqNum,:)), squeeze(paramPairedP(3,curIdx,freqNum,:)),squeeze(paramPairedP(3,curIdx,freqNum,:))<0.05] );
        motion2D3D_stats{curIdx}.Properties.VariableNames = paramNames;
        motion2D3D_stats{curIdx}.Properties.Description =sprintf('%s:%s',testNames{z},orientNames{s});
    end
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

    

    

    