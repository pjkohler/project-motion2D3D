% main_prep_m2D3D

clear all;
close all;

adultExp = [1,2,3,4,5];
plotType = 'freq';
figFolder = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures';
%figFolder = '/Users/kohler/Desktop/figures';

% PROJECTED DATA OR NOT?
projectedData = true;

%% PLOT RCs
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

nFreq = 4;
numExp = length(adultExp);
rcNum = 1;
figLabel = {'A','B','C','D','E','F','G','H','I','J'};
flipIdx = [1,-1,1,1,1;
           -1,-1,-1,-1,1;
           -1,-1,1,1,-1;
           -1,-1,-1,-1,-1];

% topo stuff
freqLabels = {'1F','2F','3F','4F'};
multX = 2.2;
multY = 2.2;
       
paperFolder = sprintf('%s/paper_figures/exp1-5',figFolder);       
if ~exist(paperFolder,'dir')
    mkdir(paperFolder);
else
end    

%% MAKE FIGURES AND SET FIGURE SIZE IN THE BEGINNING
sweepH = 24;
sweepW = 16;
topoH = 25;
topoW = 5;

for f = 1:nFreq
    sweepFig(f) = figure;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = sweepH;
    figPos(3) = sweepW;
    set(gcf,'pos',figPos);
    topoFig(f) = figure;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(4) = topoH;
    figPos(3) = topoW;
    set(gcf,'pos',figPos);
end

%% FILL FIGURES

for e = 1:numExp
    load(sprintf('%s/exp%0.0f/rcaData.mat',figFolder,adultExp(e)),'readyRCA');
    if any(e == [4,5])
        subColors = alternaColors;
    else
        subColors = mainColors;
    end
    for f = 1:nFreq
        figure(sweepFig(f));
        if strcmp(plotType,'all'); 
            curFreq = f+8;
            % grab frequency name from condition 1, 
            % will be the same for all conditions
            freqName = readyRCA(curFreq).settings.freqLabels{f};
        elseif strcmp(plotType,'full'); 
            curFreq = f+4;
            freqName = readyRCA(curFreq).settings.freqLabels{f};
        elseif strcmp(plotType,'freq'); 
            curFreq = f;
            freqName = readyRCA(curFreq).settings.freqLabels{1};
        else
            msg = sprintf('\n unknown plot type: %s\n',plotType);
            error(msg);
        end
        
        % define shared params
        if e == 1 && f == 1
            binVals = cellfun(@(x) str2num(x), readyRCA(curFreq).settings.binLabels);
            logStep = diff(reallog(binVals(1:2))); % step size
            xMin = exp(reallog(binVals(1))-logStep*.5);
            xMax = exp(reallog(binVals(end))+logStep*2.5); % add 2.5 steps
            extraBins = arrayfun(@(x) exp(reallog(binVals(end))+x), [logStep,logStep*2]);
        else
        end
        
        % compute PCA variance explained
        [ rcaRelExpl(f,e), rcaVarExpl(f,e) ,pcaVarExpl(f,e), pcaEigs(:,f,e) ] = rcaExplained(readyRCA(curFreq),1);
            
        for s=1:2 % vertical or horizontal motion
            if s==1
                sH(f,e,1) = subplot(numExp,2,1+(e-1)*2);
                curConds = find(ismember(condsToUse,1:4));
                titleStr = ['horizontal: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            else
                sH(f,e,2) = subplot(numExp,2,2+(e-1)*2);
                curConds = find(ismember(condsToUse,5:8));
                titleStr = ['vertical: \it\fontname{Arial}',sprintf('%s',freqName(1:2))];
            end
            hold on

            % store NR for later
            NakaRushton(e,f).Params = readyRCA(curFreq).stats.NR_Params;
            NakaRushton(e,f).R2 = readyRCA(curFreq).stats.NR_R2;
            NakaRushton(e,f).JKSE = readyRCA(curFreq).stats.NR_JKSE;
            NakaRushton(e,f).JKParams = readyRCA(curFreq).stats.NR_JKParams;
            
            NRset = squeeze(readyRCA(curFreq).stats.NR_Params(:,rcNum,curConds));
            NRmodel = readyRCA(curFreq).stats.hModel;
            
            % compute new signal values, averaged over bins
            [rcaDataReal,rcaDataImag] = getRealImag(readyRCA(curFreq).data(curConds,:));
            rcaDataReal = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataReal,'uni',false);
            rcaDataReal = cell2mat(permute(rcaDataReal,[3,2,1]));
            rcaDataImag = cellfun(@(x) squeeze(nanmean(x(:,rcNum,:),3)),rcaDataImag,'uni',false);
            rcaDataImag = cell2mat(permute(rcaDataImag,[3,2,1]));
            realBinMean = squeeze(nanmean(rcaDataReal));
            imagBinMean = squeeze(nanmean(rcaDataImag));
            % compute new noise values, averaged over bins
            [noiseLoReal,noiseLoImag] = getRealImag(readyRCA(curFreq).noiseData.lowerSideBand(curConds,:));
            [noiseHiReal,noiseHiImag] = getRealImag(readyRCA(curFreq).noiseData.lowerSideBand(curConds,:));
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
                valSet = squeeze(readyRCA(curFreq).stats.Amp(:,rcNum,curConds));
                errSet1 = squeeze(readyRCA(curFreq).stats.ErrLB(:,rcNum,curConds));
                errSet2 = squeeze(readyRCA(curFreq).stats.ErrUB(:,rcNum,curConds));
                % compute vector means of mean bins, and errors
                valSet(length(binVals)+1,:) = sqrt(nanmean(realBinMean,1).^2+nanmean(imagBinMean,1).^2);
                
                % store t-values for later
                rc_tSqrdP(e,f,s) = {permute(squeeze(readyRCA(curFreq).stats.tSqrdP(:,rcNum,curConds)),[2,1])};
                rc_tSqrdVal(e,f,s) = {permute(squeeze(readyRCA(curFreq).stats.tSqrdVal(:,rcNum,curConds)),[2,1])};
                
                % compute elliptical error and do Hotelling's T2 against zero
                for c=1:length(curConds)
                    xyBinMean = cat(2,realBinMean(:,c),imagBinMean(:,c));
                    nanVals = sum(isnan(xyBinMean),2)>0;      
                    [binErrs,~,~,errorEllipse] = fitErrorEllipse(xyBinMean(~nanVals,:),'SEM');
                    ellipseData(e,f,curConds(c)).ellipse = errorEllipse;
                    ellipseData(e,f,curConds(c)).means = nanmean(xyBinMean);
                    errSet1(length(binVals)+1,c) = binErrs(1);
                    errSet2(length(binVals)+1,c) = binErrs(2);
                    % compute t-values
                    tStruct = tSquaredFourierCoefs(xyBinMean(~nanVals,:));
                    rc_tSqrdP{e,f,s}(c,length(binVals)+1) = tStruct.pVal;
                    rc_tSqrdVal{e,f,s}(c,length(binVals)+1) = tStruct.tSqrd;
                    % note, just using mean df for all values, would need
                    % to be fixed if multi data was ever used seriously
                    rc_tSqrdDF{e,f,s}(c,1:length(binVals)+1) = tStruct.df2;
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
                rc_tSqrdP{e,f,s} = permute(temp_p,[3,2,1]);
                rc_tSqrdVal{e,f,s} = permute(temp_stats.tstat,[3,2,1]);
                rc_tSqrdDF{e,f,s} = permute(temp_stats.df,[3,2,1]);
                clear temp_*;
            end
            % maximum noise across all four conditions being plotted
            % since this is just means, we can compute it the same way for
            % projected and not-projected
            noiseSet = max(readyRCA(curFreq).stats.NoiseAmp(:,rcNum,curConds),[],3);
            % max of mean over bins
            noiseSet(length(binVals)+1) = max(sqrt(nanmean(realBinMeanNoise,1).^2+nanmean(imagBinMeanNoise,1).^2));

            tIdx = [1,2;3,4;1,3;2,4]; % mot ref vs no ref, disp ref vs no ref, relMot vs no relDisp, absMot vs no absDisp
            
            for pT=1:length(tIdx) % do four paired tests      
                for b = 1:(length(binVals)+1)
                    if ~projectedData
                        if b <= (length(binVals))
                            xyData = permute(cat(1,rcaDataReal(b,:,tIdx(pT,:)),rcaDataImag(b,:,tIdx(pT,:))),[2,1,3]);
                        else
                            xyData = cat(3,realBinMean(:,tIdx(pT,:)), imagBinMean(:,tIdx(pT,:)));
                        end
                        tempStrct = tSquaredFourierCoefs(xyData);
                        tempP(:,b) = tempStrct.pVal;
                        tempT(:,b) = tempStrct.tSqrd;
                        tempD(:,b) = tempStrct.mahalanobisD;
                        tempU(:,b) = tempStrct.cohenNonOverlap;
                        tempMu1(:,b) = valSet(b,tIdx(pT,1));
                        tempMu2(:,b) = valSet(b,tIdx(pT,2));
                        tempDF(:,b) = tempStrct.df2;
                    else
                        curData = squeeze(project_amps(:,b,tIdx(pT,:)));
                        not_nan = ~any(isnan(curData),2);
                        [~,tempP(:,b),~,tempStrct] = ttest(curData(not_nan,1),curData(not_nan,2),'alpha',0.05,'dim',1,'tail','both');
                        tempT(:,b) = tempStrct.tstat;
                        tempD(:,b) = tempT(:,b)./sqrt(tempStrct.df+1); % Cohen's D
                        OVL = 2*normcdf(-abs(tempD(:,b))/2);
                        tempU(:,b) = 1-OVL/(2-OVL); 
                        tempMu1(:,b) = mean(curData(not_nan,1));
                        tempMu2(:,b) = mean(curData(not_nan,2));
                        tempDF(:,b) = tempStrct.df;
                        %tempMu1(:,b) = project_valSet(b,tIdx(pT,1));
                        %tempMu2(:,b) = project_valSet(b,tIdx(pT,2));
                    end
                end
                rc_tSqrdP(e,f,s) = {cat(1,rc_tSqrdP{e,f,s},tempP)};
                rc_tSqrdVal(e,f,s) = {cat(1,rc_tSqrdVal{e,f,s},tempT)};
                rc_tSqrdDF(e,f,s) = {cat(1,rc_tSqrdDF{e,f,s},tempDF)};
                if pT == 1
                    rc_tSqrdD(e,f,s) = {tempD};
                    rc_tSqrdU(e,f,s) = {tempU};
                    rc_tSqrdMu1(e,f,s) = {tempMu1};
                    rc_tSqrdMu2(e,f,s) = {tempMu2};
                else
                    rc_tSqrdD(e,f,s) = {cat(1,rc_tSqrdD{e,f,s},tempD)};
                    rc_tSqrdU(e,f,s) = {cat(1,rc_tSqrdU{e,f,s},tempU)};
                    rc_tSqrdMu1(e,f,s) = {cat(1,rc_tSqrdMu1{e,f,s},tempMu1)};
                    rc_tSqrdMu2(e,f,s) = {cat(1,rc_tSqrdMu2{e,f,s},tempMu2)};
                end
                clear temp*
            end
                    
            for c=1:length(curConds)
                valH(c)=plot([binVals;extraBins((c>2)+1)],valSet(:,c),'o','MarkerSize',5,'LineWidth',lWidth,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
                hE = ErrorBars([binVals;extraBins((c>2)+1)],valSet(:,c),[errSet1(:,c),errSet2(:,c)],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth);
                uistack(valH(c),'bottom')
                cellfun(@(x) uistack(x,'bottom'), hE);
                hold on
                if ~isnan(NRset(1,c))
                    % plot Naka-Rushton
                    nFine = 1e2;
                    nrX = linspace( min(binVals), max(binVals), nFine )';
                    nrVals = NRmodel( nrX, NRset(:,c));
                    hNR{c} = plot( nrX, nrVals, '-','color',subColors(curConds(c),:),'LineWidth',lWidth);
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
                yUnit = 1;
                yMax = 2;
            end
            set(gca,gcaOpts{:},'XScale','log','XMinorTick','off','xtick',[0.1,0.5,1,2,4,8,16],'ytick',0:yUnit:yMax,'Layer','top','clipping','off');
            xlim([xMin,xMax]);
            ylim([0,yMax])
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
            
            text(0.2,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,figLabel{e+(s-1)*5},'fontsize',fSize*1.5,'fontname','Helvetica');
            text(0.45,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.05,sprintf('Exp. %0.0f',adultExp(e)),'fontsize',fSize,'fontname','Helvetica');
            if projectedData
                numSubs = rc_tSqrdDF{e,f,s}(1,end) + 1;
            else
                numSubs = rc_tSqrdDF{e,f,s}(1,end) + 2;
            end
            text(0.45,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.18,sprintf('n = %0d',numSubs),'fontsize',fSize,'fontname','Helvetica');
            
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
            cur_axis = gca;
            hold off
            % plot icons
            hold on
            if s == 1
                icon_location = sprintf('%s/icons/exp%02d_h.png',paperFolder,e);
            else
                icon_location = sprintf('%s/icons/exp%02d_v.png',paperFolder,e);
            end
            % load image
            [im_data,im_cmap,im_alpha] = imread(icon_location);
            im_data(repmat(im_alpha,1,1,3)==0) = 255;
            % find image dimensions, given current plot
            ax_pos = get(cur_axis,'position');
            im_aspect = size(im_data,2)/size(im_data,1);
            im_height = 1;
            plot_aspect = ax_pos(4)/ax_pos(3);
            im_width = im_height*im_aspect;%*plot_aspect;            
            % create new axis
            im_axis = axes('position', ax_pos);
            % put axis scale on plot scale
            
            imX = [1-im_width,1];
            imY = [0,1];
            set(im_axis,'xdir','normal','ydir','normal');
            imH = image(imX,imY,im_data,'parent',im_axis);
            axis(im_axis,'square')
            set(im_axis,'ylim',[0,1],'xlim',[0,1],'visible','off')
            set(imH, 'AlphaData', im_alpha,'AlphaDataMapping','none');
            
            im_scale = 3/4; % image scaling, relative to axis
            orig_pos = ax_pos;
            ax_pos(2) = orig_pos(2) + orig_pos(4)*.5;
            ax_pos(1) = orig_pos(1) + orig_pos(3)*.56;
            ax_pos(3:4) = orig_pos(3:4) .* im_scale;
            set(im_axis,'position',ax_pos);
            uistack(im_axis,'bottom');
            set(cur_axis,'color','none');
            set(im_axis,'color','none');
            hold off;
            if e == numExp
                % not very pretty code for putting labels on graph
                if s == 1
                    ax_pos = get(cur_axis,'position');
                    ax_pos(1) = ax_pos(1)-ax_pos(3)*.075;
                    ax_pos(2) = ax_pos(2)-ax_pos(4)*1.3;
                    ax_pos(3) = ax_pos(3)*1.2;
                    text_pos = ax_pos;
                else
                    text_pos(1) = text_pos(1)+text_pos(3);
                end
                text_axis = axes('position', text_pos);
                set(text_axis,'xdir','normal','ydir','normal','ylim',[0,1],'xlim',[0,1],'visible','off');
                textOpts = {'horizontalalignment','center','verticalalignment','middle','parent',text_axis,'margin',1,'fontangle','italic','fontname','Arial'};

                textColors = {[1,1,1],[0,0,0],[1,1,1],[0,0,0],[1,1,1],[1,1,1]};
                patchColors = {cBrewer.rgb20(7,:),cBrewer.rgb20(13,:),cBrewer.rgb20(1,:),...
                    [204,204,204]./255,[128,128,128]./255,[0,0,0]};
                patchH = .3;
                patchX = [0,1/3,1/3,0];
                patchY = [0 0 patchH patchH];
                textX = [1/6,1/2,5/6];
                textStr = {['100% Corr', char(10), '100% Coh'],...
                           ['0% Corr', char(10), '100% Coh'],...
                           ['-100% Corr', char(10), ' 100% Coh'],...
                           ['100% Corr', char(10), '100% Coh'],...
                           ['0% Corr', char(10), ' 0% Coh'],...
                           'blank'};
                for p = 1:3
                    % plot patch
                    patch(patchX+(p-1)/3,patchY+patchH*1.5,patchColors{p+(s-1)*3},'edgecolor','none','parent',text_axis);
                    % plot text
                    text(textX(p),patchH*2,textStr{p+(s-1)*3},'color',textColors{p+(s-1)*3},textOpts{:});
                end
            else
            end
        end
        
        set(0,'CurrentFigure',topoFig(f));
        egiH(e,f) = subplot(5,1,e);
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
        text(-1.3,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.05,sprintf('Exp. %0.0f',adultExp(e)),'fontsize',fSize,'fontname','Helvetica');
        if e == 1
            text(0.9,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.05,freqLabels{f},'fontsize',fSize,'fontname','Helvetica');
        else
        end
        hold off
        if e == numExp
            % adjust positions
            for e = 1:numExp
                if e == numExp
                    set(cH(e,f),'fontsize',fSize,'fontname','Helvetica','YTick',linspace(min(rcaColorBar(:,f)),min(rcaColorBar(:,f))*-1,5));
                    xlabel(cH(e,f),'weights','fontsize',fSize,'fontname','Helvetica')
                    set(cH(e,f),'location','southoutside');
                    set(cH(e,f),'units','centimeters');
                    cBarPos = get(cH(e,f),'position');
                    cBarPos(1) = cBarPos(1)-cBarPos(1)*multX*.2;
                    cBarPos(2) = cBarPos(2)-cBarPos(2)*multY*.1;
                    cBarPos(3) = cBarPos(3) * multX*.75;
                    cBarPos(4) = cBarPos(4) * multY*.75;
                    oldPos = newPos;
                    newPos = get(egiH(e,f),'position');
                    newPos(1) = oldPos(1);
                    newPos(2) = newPos(2)-(oldPos(4)*.3);
                    newPos(3) = oldPos(3);
                    newPos(4) = oldPos(4);
                else
                    set(cH(e,f),'visible','off');
                    newPos = get(egiH(e,f),'position');
                    newPos(1) = newPos(1)-(newPos(3)*multX*.1);
                    newPos(2) = newPos(2)-(newPos(4)*multY*.1);
                    newPos(3) = newPos(3)*multX;
                    newPos(4) = newPos(4)*multY;
                end
                set(egiH(e,f),'position',newPos);
            end
            set(cH(e,f),'position',cBarPos);
        else
        end
    end
end

% now save figures
tic
for f = 1:nFreq
    if projectedData
        if strcmp(plotType,'freq') || f == 2
            export_fig(sprintf('%s/adultExp_harm%0.0f_proj_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',sweepFig(f));
            export_fig(sprintf('%s/adultExp_topo%0.0f_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',topoFig(f));
        else
        end
    else
        if strcmp(plotType,'freq') || f == 2
            export_fig(sprintf('%s/adultExp_harm%0.0f_multi_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',sweepFig(f));
            export_fig(sprintf('%s/adultExp_topo%0.0f_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',topoFig(f));
        else
        end
    end 
end
toc

if ~strcmp(plotType,'freq')
    % stop here, don't do stats
    return
else
end

%% MAKE TALK FIGURES
for f = 1:nFreq
    figure(sweepFig(f));
    for e = 1:numExp
        if e < 4
            subColors = mainColors(1:4,:);
            deleteOrder = [4,2,3,1];
        else
            subColors = alternaColors(1:4,:);
            deleteOrder = [3,1,4,2];
        end
        for s = 1:2
            subplot(numExp,2,s+(e-1)*2);
            subAx = gca;
            newFig = figure;
            newAx = gca;
            subH = get(subAx,'children'); %get handle to all the children in the figure
            copyobj(subH,gca); %copy children to new parent axes i.e. the subplot axes
            textH = findall(gca,'type','text');
            for t= 1:length(textH)
                if length(get(textH(t),'string')) == 1 
                    delete(textH(t));
                else
                end
            end
            copySettings = {'xlim','ylim','xscale','yscale','xtick','ytick','xminortick','yticklabel','fontsize','fontname','linewidth','box','ticklength','tickdir'};
            cellfun(@(x) set(gca,x,get(subAx,x)), copySettings,'uni',false)
            ylabel('amplitude (\muV)')
            xlabel('displacement (arcmins)');
            set(newFig,'units','centimeter');
            figPos = get(newFig,'pos');
            figPos(3) = sweepW/2;
            figPos(4) = sweepH/5;
            set(newFig,'position',figPos);
            for c = 1:5
                if c == 1
                    deleteColors = [];
                else
                    deleteColors = subColors(deleteOrder(c-1),:);
                end
                if ~isempty(deleteColors)
                    lineObjs = findobj(gca, 'type', 'line');
                    deleteIdx = arrayfun(@(x) ismember(get(x,'color'),deleteColors,'rows'),lineObjs);
                    delete(lineObjs(deleteIdx))
                end
                saveLocation = sprintf('/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/exp%d',e);
                if ~exist(sprintf('%s/exp%d_talkFigures',saveLocation,e),'dir')
                    mkdir(sprintf('%s/exp%d_talkFigures',saveLocation,e));
                else
                end
                if s == 1
                    export_fig(sprintf('%s/exp%d_talkFigures/exp%d_harm%0.0f_%s_hori%0.0f.pdf',saveLocation,e,e,f,plotType,(5-c)+1),'-pdf','-transparent',newFig);
                else
                    export_fig(sprintf('%s/exp%d_talkFigures/exp%d_harm%0.0f_%s_vert%0.0f.pdf',saveLocation,e,e,f,plotType,(5-c)+1),'-pdf','-transparent',newFig);
                end
            end
            close(newFig);
        end
    end
end


close all


%% MAKE ELLIPSE FIGURES
if ~projectedData
    for f = 2%1:nFreq
        ellipseFig(f) = figure;
        set(gcf, 'units', 'centimeters');
        figPos = get(gcf,'pos');
        figPos(4) = 25;
        figPos(3) = 8;
        set(gcf,'pos',figPos);
        for e = 1:numExp
            if any(e == [4,5])
                subColors = alternaColors;
            else
                subColors = mainColors;
            end
            for s=1:2 % vertical or horizontal motion
                if s==1
                    sH(f,e,1) = subplot(numExp,2,1+(e-1)*2);
                    curConds = 1:4;
                else
                    sH(f,e,2) = subplot(numExp,2,2+(e-1)*2);
                    curConds = 5:8;
                end
                hold on
                axis equal;
                for c = 1:length(curConds)
                    curMeans = ellipseData(e,f,curConds(c)).means;
                    curEllipse = ellipseData(e,f,curConds(c)).ellipse;
                    if e > 2
                        tempAngle = angle(complex(curMeans(1),curMeans(2)))+pi;
                        tempComplex = abs(complex(curMeans(1),curMeans(2))).*exp(1i*tempAngle);
                        curMeans = [real(tempComplex),imag(tempComplex)];
                        tempAngle = angle(complex(curEllipse(:,1),curEllipse(:,2)))+pi;
                        tempComplex = abs(complex(curEllipse(:,1),curEllipse(:,2))).*exp(1i*tempAngle);
                        curEllipse = [real(tempComplex),imag(tempComplex)];
                    end
                    plot([0 curMeans(1)],[0 curMeans(2)],'-','Color',subColors(curConds(c),:),'LineWidth',1);
                    %plot(ellipseData(e,f,c).means(1),ellipseData(e,f,c).means(2),'x','Color',subColors(curConds(c),:),'LineWidth',1);  
                    plot(curEllipse(:,1), curEllipse(:,2),'-','Color',subColors(curConds(c),:),'LineWidth',1); 
                end
                %drawnow;
                set(gca,gcaOpts{:},'clipping','off','visible','on');
                xlim([-4,1]);
                ylim([-1,4]);
                text(-4.25,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,figLabel{e+(s-1)*5},'fontsize',fSize*1.5,'fontname','Helvetica');
                text(-3.75,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.05,sprintf('Exp. %0.0f',adultExp(e)),'fontsize',fSize,'fontname','Helvetica');
                text(-3.75,max(get(gca,'ylim'))-diff(get(gca,'ylim'))*.18,sprintf('n = %0d',degFree(e,f)),'fontsize',fSize,'fontname','Helvetica');
                if s == 1 && e == 5
                    xlabel('real','fontsize',fSize,'fontname','Helvetica');
                    ylabel('imag','fontsize',fSize,'fontname','Helvetica');
                else
                end
                hold off
            end
        end
        export_fig(sprintf('%s/adultExp_ellipse%0.0f.pdf',paperFolder,f),'-pdf','-transparent',gcf);
    end
end

%% MAKE PCA FIGURES
close all
figure;
for f = 1:2
    for e = 1:5
        subplot(5,2,f+(e-1)*2);
        plot(pcaEigs(1:50,f,e),'-o');
        title(sprintf('exp %0.0f: %s',e,freqLabels{f}),'fontsize',fSize,'fontname','Helvetica');
        set(gca,gcaOpts{:});
        xlim([0.5,50.5]);
    end
    if f == 1
        expIdx = [1,2,4,5]; % exclude no signal experiment 3
    else
        expIdx = 1:5;
    end
    meanVarExpl(f) = mean(rcaVarExpl(f,expIdx),2);
    stdVarExpl(f) = std(rcaVarExpl(f,expIdx),0,2);
    meanRelExpl(f) = mean(rcaRelExpl(f,expIdx),2);
    stdRelExpl(f) = std(rcaRelExpl(f,expIdx),0,2);
end
set(gcf, 'units', 'centimeters');
figPos = get(gcf,'pos');
figPos(4) = sweepH;
figPos(3) = sweepW;
set(gcf,'pos',figPos);
export_fig(sprintf('%s/paper_figures/exp1-5/adultExp_pca.pdf',figFolder),'-pdf','-transparent',gcf);
close gcf;

%% MAKE STATS FIGURES
close all;
fSize = 10;
lWidth = 1.5;
gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'box','off','fontsize',fSize,'fontname','Helvetica','linewidth',lWidth};

figHeight = 25;
figWidth = 17.8*(2/3);
numExp = 5;
numTests = 8; % 4 conds, in-phase vs anti-phase, absolute vs relative
xVals = 1:11;
testSymbols = {'o','o','o','o','o','o','o','o'};
markerOpts = {'MarkerSize',6,'LineWidth',lWidth,'markerfacecolor',[1,1,1]};
freqLabels = {'1F','2F','3F','4F'};
% significance colors
tHotColMap = jmaColors('pval');
tHotColMap(end,:) = [1 1 1];
pCutOff = 0.05;
for e = 1:numExp
    if any(e == [4,5])
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
            hImg = image([1,11],[1,numTests],curP, 'CDataMapping', 'scaled','Parent',gca);
            colormap(gca,tHotColMap );   
            cMapMax = pCutOff+2*pCutOff/(size(tHotColMap,1));
            set(gca,'xtick',2:2:10,'xticklabel',arrayfun(@(x) sprintf('%0.1f',binVals(x)),2:2:10,'uni',false),...
                'ytick',[],'ycolor',[1,1,1],'yaxislocation','left','yticklabel','','CLim', [ 0 cMapMax ], gcaOpts{:} ); % set range for color scale
            xlim([-1.5,11.5]);
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
    if projectedData
        export_fig(sprintf('%s/adultExp_stats%0.0f_proj_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',gcf);
    else
        export_fig(sprintf('%s/adultExp_stats%0.0f_multi_%s.pdf',paperFolder,f,plotType),'-pdf','-transparent',gcf);
    end
end

%% MAKE NAKA-RUSHTON FIGURES
close all;
for e = 1:numExp
    if any(e == [4,5])
        subColors = alternaColors;
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
            paramList = [1,2,3,4];
            for pT=1:length(tIdx) % do four paired tests
                 % only look at c50 and rMax
                diffErr = jackKnifeErr(testVal(:,paramList,tIdx(pT,1))-testVal(:,paramList,tIdx(pT,2)));
                grandDiff = NRvals(paramList,tIdx(pT,1)) - NRvals(paramList,tIdx(pT,2));
                paramPairedT(:,pT+(s-1)*4,f,e) = grandDiff'./diffErr;
                paramPairedP(:,pT+(s-1)*4,f,e) = 2*tcdf( -abs(paramPairedT(:,pT+(s-1)*4,f,e)) , jkDf);
                paramPairedDf(:,pT+(s-1)*4,f,e) = jkDf;
            end
            
            subPlotOrder = [1,2,5,6];
            paramList = [paramList,6]; % add R2
            for z = 1:length(paramList) % five parameters (including R2)
                sH(f,e,s) = subplot(5,1,z);
                hold on
                xVals = (1:4)+5*(s-1)+10*(e-1);
                arrayfun(@(x) bar(xVals(x),NRvals(paramList(z),x),'facecolor',subColors(x,:),'edgecolor','none'),1:4,'uni',false);
                if paramList(z) ~= 6
                    eH = arrayfun(@(x) ...
                        ErrorBars(xVals(x),NRvals(paramList(z),x),NRerrors(paramList(z),x),'color',[0,0,0],'type','bar','cap',false,'barwidth',lWidth),...
                        1:4,'uni',false);
                   cellfun(@(x) set(x{:},'clipping','on'),eH);
                else
                end
                if e == numExp
                    switch paramList(z)
                        case 1
                            ylabel('d_5_0');
                            yMin = 0; yMax = 8; yUnit = 2;
                        case 2
                            ylabel('n');
                            yMin = 0; yMax = 8; yUnit = 2;
                        case 3
                            ylabel('R_m_a_x');
                            yMin = 0; yMax = 4; yUnit = 1;
                        case 4
                            ylabel('b');
                            yMin = 0; yMax = 2; yUnit = 1;
                        case 6
                            ylabel('R^2');
                            yMin = 0; yMax = 1; yUnit = 0.5;
                        otherwise
                    end
                    gcaOpts{4} = [0.01,0.01];
                    ylim([yMin,yMax]);
                    set(gca, gcaOpts{:});               
                    if z == 1
                        set(gca,'YTickLabel',sprintf('%0.0f|',yMin:yUnit:yMax))
                        set(gca,'xtick',[],'ytick',yMin:yUnit:yMax);
                        arrayfun(@(x) text(5+(x-1)*10,max(get(gca,'ylim'))*1.1, ...
                                [sprintf('Exp. %0.0f: ',adultExp(x)),...
                                '\it\fontname{Arial}',sprintf('%s',freqName(1:2))],...
                                'horizontalalignment','center'),...
                                1:numExp,'uni',false);
                    elseif z == length(paramList)
                        set(gca,'YTickLabel',sprintf('%0.0f|',yMin:yUnit:yMax))
                        set(gca,'xtick',[2.5:5:50],'ytick',yMin:yUnit:yMax,'xticklabel',{'Hori','Vert'},'clipping','off');
                    else
                        set(gca,'YTickLabel',sprintf('%0.0f|',yMin:yUnit:yMax))
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
    figPos(4) = figWidth*(5/4);
    figPos(3) = figHeight;
    set(gcf,'pos',figPos);
    export_fig(sprintf('%s/adultExp_Naka%0.0f.pdf',paperFolder,f),'-pdf','-transparent',gcf);
end

close all;
%% MAKE TABLES

% Naka-Rushtin table

paramNames = {'D50','n','Rmax','b','df'};
binNames = {'T2','p','D','NonO','mean1','mean2','df'};

orientNames = {'Horizontal','Vertical'};
testList =  {'In_RefVsNo', 'Anti_RefVsNo', 'Ref_InVsAnti', 'NoRef_InVsAnti','IOVD_InVsAnti'};
expList =  {'Exp. 1','Exp. 2','Exp. 3','Exp. 4','Exp. 5'};

% add the t2 on mean bin-vals, to params
tempT = cellfun(@(x) x(5:end,end), rc_tSqrdVal,'uni',false);
tempP = cellfun(@(x) x(5:end,end), rc_tSqrdP,'uni',false);
tempD = cellfun(@(x) x(:,end), rc_tSqrdD,'uni',false);
tempU = cellfun(@(x) x(:,end), rc_tSqrdU,'uni',false);
tempMu1 = cellfun(@(x) x(:,end), rc_tSqrdMu1,'uni',false);
tempMu2 = cellfun(@(x) x(:,end), rc_tSqrdMu2,'uni',false);
for e = 1:numExp
    for f = 1:nFreq
        paramPairedP(5,:,f,e) = cat(1,tempP{e,f,:});
        paramPairedT(5,:,f,e) = cat(1,tempT{e,f,:});
        paramPairedD(1,:,f,e) = cat(1,tempD{e,f,:});
        paramPairedU(1,:,f,e) = cat(1,tempU{e,f,:});
        paramPairedMu1(1,:,f,e) = cat(1,tempMu1{e,f,:});
        paramPairedMu2(1,:,f,e) = cat(1,tempMu2{e,f,:});
    end
end

% convert p-vals to strings
stringPairedP = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedP,'uni',false);
stringPairedT = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedT,'uni',false);
stringPairedD = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedD,'uni',false);
stringPairedU = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedU,'uni',false);
stringPairedMu1 = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedMu1,'uni',false);
stringPairedMu2 = arrayfun(@(x) num2str(x,'%0.4f'),paramPairedMu2,'uni',false);

sigIdx = cell2mat(arrayfun(@(x) x < 0.0001,paramPairedP,'uni',false));
stringPairedP(sigIdx) = {'<0.0001'};   

% make table of paired results
freqNum = 2;
for z = 1:length(testList)
    switch testList{z}
        case 'Ref_InVsAnti'
            expIdx = 1:5;
            testIdx = 3;
        case 'IOVD_InVsAnti'
            expIdx = [4,5];
            testIdx = 4;
        case 'NoRef_InVsAnti'
            expIdx = 2;
            testIdx = 4;
        case 'In_RefVsNo'
            expIdx = 1:3;
            testIdx = 1;
        case 'Anti_RefVsNo'
            expIdx = 1:3;
            testIdx = 2;
        otherwise
    end        
    for s = 1:2
        testIdx = testIdx+(s-1)*4; % horizontal and vertical
        curIdx = z+(s-1)*length(testList);
        motion2D3D_stats{curIdx} = table(...
            expList(expIdx)',...
            [squeeze(stringPairedT(1,testIdx,freqNum,expIdx)), squeeze(stringPairedP(1,testIdx,freqNum,expIdx))], ...
            [squeeze(stringPairedT(2,testIdx,freqNum,expIdx)), squeeze(stringPairedP(2,testIdx,freqNum,expIdx))], ...
            [squeeze(stringPairedT(3,testIdx,freqNum,expIdx)), squeeze(stringPairedP(3,testIdx,freqNum,expIdx))], ...
            [squeeze(stringPairedT(4,testIdx,freqNum,expIdx)), squeeze(stringPairedP(4,testIdx,freqNum,expIdx))],...
            squeeze(paramPairedDf(1,testIdx,freqNum,expIdx)));
        motion2D3D_stats{curIdx}.Properties.VariableNames = [orientNames{s},paramNames];
        motion2D3D_stats{curIdx}.Properties.Description =sprintf('%s:%s',testList{z},orientNames{s});
        csv_nr{s} = sprintf('%s/paper_figures/exp1-5/NakaRush_tables/nrTest2D3D_2F_%s_%02d.csv',figFolder,testList{z},s);
        writetable(motion2D3D_stats{curIdx},csv_nr{s},'WriteRowNames',true);
        % now do bin stats
        motion2D3D_bins{curIdx} = table(...
            expList(expIdx)',...
            squeeze(stringPairedT(5,testIdx,freqNum,expIdx)), ...
            squeeze(stringPairedP(5,testIdx,freqNum,expIdx)), ...
            squeeze(stringPairedD(1,testIdx,freqNum,expIdx)), ...
            squeeze(stringPairedU(1,testIdx,freqNum,expIdx)), ...
            squeeze(stringPairedMu1(1,testIdx,freqNum,expIdx)), ...
            squeeze(stringPairedMu2(1,testIdx,freqNum,expIdx)), ...
            squeeze(paramPairedDf(1,testIdx,freqNum,expIdx)));
        motion2D3D_bins{curIdx}.Properties.VariableNames = [orientNames{s},binNames];
        motion2D3D_bins{curIdx}.Properties.Description =sprintf('%s:%s',testList{z},orientNames{s});
        csv_ttest{s} = sprintf('%s/paper_figures/exp1-5/ttest_tables/t2Test2D3D_2F_%s_%02d.csv',figFolder,testList{z},s);
        writetable(motion2D3D_bins{curIdx},csv_ttest{s},'WriteRowNames',true);
    end
    combined_nr = sprintf('%s/paper_figures/exp1-5/NakaRush_tables/nrTest2D3D_2F_%s_combined.csv',...
        figFolder,testList{z});
    if exist(combined_nr,'file')
        delete(combined_nr);
    else
    end
    combined_ttest = sprintf('%s/paper_figures/exp1-5/ttest_tables/tTest2D3D_2F_%s_combined.csv',...
        figFolder,testList{z});
    if exist(combined_ttest,'file')
        delete(combined_ttest);
    else
    end
    catCmd{1} = sprintf('cat %s %s > %s',csv_nr{1},csv_nr{2},combined_nr);
    catCmd{2} = sprintf('cat %s %s > %s',csv_ttest{1},csv_ttest{2},combined_ttest);
    catCmd{3} = sprintf('rm %s; rm %s' ,csv_nr{1},csv_nr{2});
    catCmd{4} = sprintf('rm %s; rm %s' ,csv_ttest{1},csv_ttest{2});  
    system(strjoin(catCmd, '; '));
end

%% MAKE BIN-BY-BIN TABLE

dfTest = cellfun(@(x) all(x(:,end) == x(1,1)), rc_tSqrdDF);

if ~all(dfTest(:))
    error('different dfs for different tests');
else
end

% make table of paired results
freqNum = 2;
for z = 1:length(testList)
    switch testList{z}
        case 'Ref_InVsAnti'
            expIdx = 1:5;
            testIdx = 3;
        case 'IOVD_InVsAnti'
            expIdx = [4,5];
            testIdx = 4;
        case 'NoRef_InVsAnti'
            expIdx = 2;
            testIdx = 4;
        case 'In_RefVsNo'
            expIdx = 1:3;
            testIdx = 1;
        case 'Anti_RefVsNo'
            expIdx = 1:3;
            testIdx = 2;
        otherwise
    end 
    finishedArray = [];
    stringBinVals = cat(1,arrayfun(@(x) num2str(x,'%0.2f'),binVals,'uni',false),{'n/a'});
    for s = 1:2
        numBinP = cell2mat(cellfun(@(x) x(4+testIdx,:),rc_tSqrdP(expIdx,freqNum,s),'uni',false));
        stringBinP = arrayfun(@(x) num2str(x,'%0.4f'), numBinP,'uni',false);
        sigIdx = cell2mat(arrayfun(@(x) x < 0.0001,numBinP,'uni',false));
        stringBinP(sigIdx) = {'<0.0001'};
        numBinT = cell2mat(cellfun(@(x) x(4+testIdx,:),rc_tSqrdVal(expIdx,freqNum,s),'uni',false));
        stringBinT = arrayfun(@(x) num2str(x,'%0.4f'), numBinT,'uni',false);
        numBinD = cell2mat(cellfun(@(x) x(testIdx,:),rc_tSqrdD(expIdx,freqNum,s),'uni',false));
        stringBinD = arrayfun(@(x) num2str(x,'%0.4f'), numBinD,'uni',false);
        numBinMu1 = cell2mat(cellfun(@(x) x(testIdx,:),rc_tSqrdMu1(expIdx,freqNum,s),'uni',false));
        stringBinMu1 = arrayfun(@(x) num2str(x,'%0.4f'), numBinMu1,'uni',false);
        numBinMu2 = cell2mat(cellfun(@(x) x(testIdx,:),rc_tSqrdMu2(expIdx,freqNum,s),'uni',false));
        stringBinMu2 = arrayfun(@(x) num2str(x,'%0.4f'), numBinMu2,'uni',false);
        readyArray = stringBinVals';
        rowLabels = {'disp (arcmins)'};
        for q = 1:length(expIdx)
            readyArray = [readyArray; stringBinP(q,:); stringBinT(q,:); stringBinD(q,:)];
            rowLabels = {rowLabels{:},num2str(expIdx(q),'Exp. %0.0f p'),num2str(expIdx(q),'Exp %0.0f t-statistic'),num2str(expIdx(q),'Exp %0.0f Cohen''s D')};
        end
        readyArray = cat(2,rowLabels',readyArray);
        varLabels = cat(2,{sprintf('%s (df=%0.0f)',orientNames{s},rc_tSqrdDF{1}(1,end))},...
              arrayfun(@(x) num2str(x,'bin%0.0f'),1:length(binVals),'uni',false),...
              {'ave'});
        readyArray = cat(1,varLabels,readyArray);
        finishedArray = cat(1,finishedArray,readyArray);
    end
    combined_bin_file= sprintf('%s/paper_figures/exp1-5/bin_tables/bin_tTest2D3D_2F_%s_combined.csv',...
        figFolder,testList{z});
    if exist(combined_bin_file,'file')
        delete(combined_bin_file);
    else
    end
    binTable = array2table(finishedArray);
    writetable(binTable,combined_bin_file,'WriteRowNames',false,'WriteVariableNames',false);
end
   

    

    

    