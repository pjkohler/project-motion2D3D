%% ADD PATHS
% codeFolder = '/Users/kohler/code';
% rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
% addpath(genpath(rcaCodePath));
% addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
% addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
% setenv('DYLD_LIBRARY_PATH','')

codeFolder = '/Users/labmanager/Desktop/LabManager/MatAnal';
rcaCodePath = sprintf('%s/rcaBase',codeFolder);
addpath(genpath(rcaCodePath));
addpath(genpath(sprintf('%s/mrC',codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
addpath(genpath(sprintf('%s/export_fig',codeFolder)));
setenv('DYLD_LIBRARY_PATH','')

%mainPath = '/Users/labmanager/Desktop/LabManager/WM_Data/2017_2DBaby/';
mainPath = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D';
figureFolder = sprintf('%s/figures/exp3',mainPath);

%% IDENTIFY DATA LOCATION
folderNames=[];
expIdx=[];
dataLocation = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/exp3_infants';
%dataLocation = '/Users/labmanager/Desktop/LabManager/WM_Data/2017_2DBaby/2DBaby_Data/TESTDATA';
folderNames = [folderNames; subfolders(sprintf('%s/*20*',dataLocation),1)]; %change back to 128 between wild cards
expIdx = [expIdx; ones(size(subfolders(sprintf('%s/*20*',dataLocation))))];

% dataLocation = subfolders(sprintf('%s/exp3_infants',dataLocation),1); %setup for working from Denali 
% folderNames=[folderNames;subfolders(sprintf('%s/*20*',dataLocation),1)];
% expIdx=[expIdx; ones(size(subfolders(sprintf('%s/*20*',dataLocation))))];
    
     for f = 1:length(folderNames)
         tempFolders = subfolders(folderNames{f},1);
         tempFolders = tempFolders(cellfun(@(x) ~isempty(strfind(x,'Exp')),tempFolders));
         folderNames{f} = sprintf('%s',tempFolders{end});
         %folderNames{f} = tempFolders;
     end

%%  SET UP INPUTS
binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= 2; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
condsToUse = 1:4; % if you want to include all conditions, create a vector here listing all condition numbers
condSep = [3,7];
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
nFreq = length(freqsToUse);
nCond = length(condsToUse);
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = true; % generate source data for first instance of rca?
keepConditions = true;
errorType = 'SEM';
doNR = true(1,5,4); % 5 freqs, 5 RCs + OZ, 4 conditions
trialError = true;
doExp = 1;

%%  RUN RCA
for doExp = min(expIdx):max(expIdx)
    if ~exist(figureFolder,'dir')
        mkdir(figureFolder);
    else
    end
    
    subsToUse = ismember(expIdx,doExp);
    
    %do RCA on all DATA
    allRCA = rcaSweep(folderNames, binsToUse, freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,true);
    export_fig(sprintf('%s/BabyCombined%dF%dF_cov.pdf',figureFolder,min(freqsToUse),max(freqsToUse)),'-pdf','-transparent',gcf);
    
    for f=1:length(freqsToUse)
        freqRCA(f) = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,true);
        export_fig(sprintf('%s/BabyCombinedRCA%dF1%dF1_cov.pdf',figureFolder,min(freqsToUse),max(freqsToUse)),'-pdf','-transparent',gcf);
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/BabyCombined%dRCA_cov.pdf',figureFolder,freqsToUse(f)),'-pdf','-transparent',rcaH);
        close all;
    end
    
    allRCA.inputData = allRCA.inputData(:,subsToUse);
    nanDims = [1,2]; % if all time points are zero, or all channels are zero
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};
    
    for z=1:length(structVars)
        if strfind(lower(structVars{z}),'noise')
            for n = 1:length(noiseVars)
                allRCA.(structVars{z}).(noiseVars{n}) = allRCA.(structVars{z}).(noiseVars{n})(:,subsToUse);
                allRCA.(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}).(noiseVars{n}),'uni',false);
                for f=1:length(freqsToUse)
                    freqRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                end
            end
        else
            % populate allRCA
            allRCA.(structVars{z}) = allRCA.(structVars{z})(:,subsToUse);
            allRCA.(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}),'uni',false);
            for f=1:length(freqsToUse)
                freqRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}),'uni',false);
            end
        end
    end
    
    % now save the data, seperately for each experiment
    saveFileName = sprintf('%s/rcaDataBaby%dF',figureFolder, freqsToUse);
    if ~exist([saveFileName,'.mat'], 'file')
        save(saveFileName,'*ToUse') % create file if it does not exist
    else
        save(saveFileName,'*ToUse','-append') % otherwise append
    end
    save(saveFileName,'allRCA','-append');
    save(saveFileName,'freqRCA','-append')
end

%% Compute Values for Plotting 
for f=1:length(freqsToUse)+1
    fprintf('%d\n',f);
    if f==length(freqsToUse)+1
        curRCA = allRCA;
        %curIdx = (length(freqsToUse)+1):(length(freqsToUse)+2);
        curIdx = (length(freqsToUse)+1);
    else
        curRCA = freqRCA(f);
        curIdx = f;
    end
    
    tempDataStrct = aggregateData(curRCA.data,curRCA.settings,keepConditions,errorType,trialError,doNR);
    ampVals(:,curIdx,1:nComp,:) = tempDataStrct.ampBins;
    errLB(:,curIdx,1:nComp,:)   =tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,curIdx,1:nComp,:)   =tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(curRCA.noiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    tempNoiseStrct2 = aggregateData(curRCA.noiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    [snrVals(:,curIdx,1:nComp,:),noiseVals(:,curIdx,1:nComp,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
    NR_pOpt(:,curIdx,1:nComp,:) = tempDataStrct.NakaRushton.pOpt;
    NR_JKSE(:,curIdx,1:nComp,:) = tempDataStrct.NakaRushton.JKSE;
    NR_R2(curIdx,1:nComp,:) = tempDataStrct.NakaRushton.R2;
    if ~isempty(tempDataStrct.NakaRushton.hModel) && ~exist('hModel','var')
        hModel = tempDataStrct.NakaRushton.hModel;
    else
    end
    delete(gcp('nocreate'));
    
    % COMPARISON
    tempDataStrct = aggregateData(curRCA.comparisonData,curRCA.settings,keepConditions,errorType,trialError,doNR);
    ampVals(:,curIdx,(nComp+1),:) = tempDataStrct.ampBins;
    errLB(:,curIdx,(nComp+1),:)=tempDataStrct.ampBins-tempDataStrct.ampErrBins(:,:,:,:,1);
    errUB(:,curIdx,(nComp+1),:)=tempDataStrct.ampErrBins(:,:,:,:,2)-tempDataStrct.ampBins;
    tempNoiseStrct1 = aggregateData(curRCA.comparisonNoiseData.lowerSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    tempNoiseStrct2 = aggregateData(curRCA.comparisonNoiseData.higherSideBand,curRCA.settings,keepConditions,errorType,trialError,doNR);
    [snrVals(:,curIdx,(nComp+1),:),noiseVals(:,curIdx,6,:)] = computeSnr(tempDataStrct,tempNoiseStrct1,tempNoiseStrct2,false);
    NR_pOpt(:,curIdx,(nComp+1),:) = tempDataStrct.NakaRushton.pOpt;
    NR_JKSE(:,curIdx,(nComp+1),:) = tempDataStrct.NakaRushton.JKSE;
    NR_R2(curIdx,(nComp+1),:) = tempDataStrct.NakaRushton.R2;
    delete(gcp('nocreate'));
end

for f= length(allRCA.data(:,1)):-1:1
    if isnan(allRCA.data{f,1}(1,1,1))
        allRCA.data(f,:)=[];
        allRCA.comparisonData(f,:)=[];
    end
end

freqRCA((length(freqsToUse)+1):(length(freqsToUse)*2)) = allRCA;

save(sprintf('%s/BabyDataOutput_%dF1', figureFolder, freqsToUse));