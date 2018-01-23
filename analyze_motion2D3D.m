function analyze_motion2D3D(doExp,trialError,plotType,talkFigs)
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
        doExp = [1,2,3,4,5];
    else
    end
    if nargin < 2
        trialError = false;
    else
    end
    if nargin < 3
        plotType = 'freq';
    else
    end
    if nargin < 4
        talkFigs = false;
    else
    end
    
    if numel(doExp) > 1
        for z = 1:length(doExp)
            analyze_motion2D3D(doExp(z),trialError,plotType,talkFigs);
        end
        return;
    else
    end
    
    fprintf('\n ... running exp %0.0f ...\n',doExp);
    
    saveLocation = sprintf('/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/exp%d',doExp);
    saveFileName = sprintf('%s/rcaData.mat',saveLocation);
    load(saveFileName,'freqRCA');

    %% COMPUTE VALUES FOR PLOTTING
    readyRCA = freqRCA;
    readyRCA(5:8) = freqRCA(5); % repeat fullRCA, this is only done to preserve topographies in the code below
    readyRCA(9:12) = freqRCA(6); % repeat allRCA, this is only done to preserve topographies in the code below
    clear freqRCA;
    keepConditions = true;
    errorType = 'SEM';
    doNR = false(4,8,8); % 4 freqs, 8 RCs (with comparison), 8 conditions
    doNR([1,2,4],1,:) = true; % do fitting for first RC, first, second and fourth harmonic, all conditions
    condsToUse = 1:8;
    for f=1:6
        fprintf('%d\n',f);
        if f==6
            idxList = 9:12;
        elseif f==5
            idxList = 5:8;
        else
            idxList = f;
        end
        rcStruct = aggregateData(readyRCA(f),keepConditions,errorType,trialError,doNR);
        for i = 1:length(idxList)
            curIdx = idxList(i);
            % RC
            readyRCA(curIdx).stats.Amp = squeeze(rcStruct.ampBins(:,i,:,:));
            readyRCA(curIdx).stats.SubjectAmp = squeeze(rcStruct.subjectAmp(:,i,:,:,:));
            readyRCA(curIdx).stats.ErrLB = squeeze(rcStruct.ampErrBins(:,i,:,:,1));
            readyRCA(curIdx).stats.ErrUB = squeeze(rcStruct.ampErrBins(:,i,:,:,2));
            readyRCA(curIdx).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins(:,i,:,:));
            readyRCA(curIdx).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise(:,i,:,:,:));
            % Naka-Rushton
            readyRCA(curIdx).stats.NR_Params = squeeze(rcStruct.NakaRushton.Params(:,i,:,:));
            readyRCA(curIdx).stats.NR_R2 = squeeze(rcStruct.NakaRushton.R2(:,i,:,:));
            readyRCA(curIdx).stats.NR_JKSE = squeeze(rcStruct.NakaRushton.JackKnife.SE(:,i,:,:));
            readyRCA(curIdx).stats.NR_JKParams = squeeze(rcStruct.NakaRushton.JackKnife.Params(:,:,i,:,:));
            readyRCA(curIdx).stats.hModel = rcStruct.NakaRushton.hModel;
            % t-values
            readyRCA(curIdx).stats.tSqrdP = squeeze(rcStruct.tSqrdP(:,i,:,:));
            readyRCA(curIdx).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig(:,i,:,:));
            readyRCA(curIdx).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal(:,i,:,:)); 
        end
    end
    % shut down parallel pool, which was used for fitting Naka-Rushton
    delete(gcp('nocreate'));

    %% PLOT RCs
    close all;
    nFreq = length(readyRCA(end).settings.freqsToUse);
    nComp = readyRCA(end).settings.nComp;
    plotSNR = false;
    % plot settings
    % set figure size in the beginning
    figHeight = 30;
    figWidth = 20;

    binVals = cellfun(@(x) str2num(x), readyRCA(end).settings.binLevels{1});

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

    

    for f=1:nFreq
        if strcmp(plotType,'all'); 
            curFreq = f+8;
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
        figure;
        for r = 1:(nComp+1)

            if r<6
                egiH(r) = subplot(nComp+1,3,3+(r-1)*3);
                hold on
                rcaColorBar = [min(readyRCA(curFreq).A(:,r)),max(readyRCA(curFreq).A(:,r))];
                newExtreme = max(abs(rcaColorBar));
                rcaColorBar = [-newExtreme,newExtreme];
                mrC.plotOnEgi(readyRCA(curFreq).A(:,r),rcaColorBar);
                hold off
            else
            end
            %title(['Full RCA ' num2str(r)],'fontsize',fSize,'fontname','Arial');
            for s=1:2 % vertical or horizontal motion
                if s==1
                    sH(f,r,1) = subplot(nComp+1,3,1+(r-1)*3);
                    curConds = find(ismember(condsToUse,1:4));
                    % grab frequency labels from condition 1, 
                    % will be the same for all conditions
                    titleStr = sprintf('horizontal: %s',freqName);
                else
                    sH(f,r,2) = subplot(nComp+1,3,2+(r-1)*3);
                    curConds = find(ismember(condsToUse,5:8));
                    % grab frequency labels from condition 1, 
                    % will be the same for all conditions
                    titleStr = sprintf('vertical: %s',freqName);
                end
                hold on
  
                valSet = readyRCA(curFreq).stats.Amp;
                errSet1 = readyRCA(curFreq).stats.ErrLB;
                errSet2 = readyRCA(curFreq).stats.ErrUB;
                NRset = readyRCA(curFreq).stats.NR_Params;
                NRmodel = readyRCA(curFreq).stats.hModel;
                for c=1:length(curConds)
                    valH(c)=plot(binVals,valSet(:,r,curConds(c)),'o','MarkerSize',5,'LineWidth',lWidth*1.5,'Color',subColors(curConds(c),:),'markerfacecolor',[1,1,1]);
                    hE = ErrorBars(binVals,valSet(:,r,curConds(c)),[errSet1(:,r,curConds(c)),errSet1(:,r,curConds(c))],'color',subColors(curConds(c),:),'type','bar','cap',false,'barwidth',lWidth*1.5);
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
                meanNoise = max(readyRCA(curFreq).stats.NoiseAmp(:,r,curConds),[],3)';
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
                        ylabel('Amplitude (\muV)')
                        xlabel('Distance (arcmins)');
                    else
                        lH = legend(valH,condLabels(curConds),'location','northeast');
                        legend boxoff
                        lPos = get(lH,'position');
                        lPos(1) = lPos(1) + .2;
                        lPos(2) = lPos(2) + .05;
                        set(lH,'position',lPos);            
                        tH = text(60,-0.1,sprintf('n = %0d',size(readyRCA(curFreq).data,2)),'fontsize',fSize,'fontname','Arial');
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
            export_fig(sprintf('%s/exp%d_rc%d_%s_snr.pdf',saveLocation,doExp,f,plotType),'-pdf','-transparent',gcf);
        else
            export_fig(sprintf('%s/exp%d_rc%d_%s.pdf',saveLocation,doExp,f,plotType),'-pdf','-transparent',gcf);
        end
        
        if f < 3 && talkFigs
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
                        export_fig(sprintf('%s/exp%d_talkFigures/exp%d_rc%d_%s_hori%d.pdf',saveLocation,doExp,doExp,f,plotType,z),'-pdf','-transparent',fH);
                    else
                        export_fig(sprintf('%s/exp%d_talkFigures/exp%d_rc%d_%s_vert%d.pdf',saveLocation,doExp,doExp,f,plotType,z),'-pdf','-transparent',fH);
                    end
                    close(fH);
                end
            end
        else
        end
    end
    save(sprintf('%s/plottingData.mat',saveLocation),'readyRCA');
    close all;
end
