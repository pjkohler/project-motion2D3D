%% ADD PATHS
clear all except codeFolder; 
close all;
if ~exist('codeFolder','var')
    codeFolder = '/Users/kohler/code';
    rcaCodePath = sprintf('%s/git/rcaBase',codeFolder);
    addpath(genpath(rcaCodePath));
    addpath(genpath(sprintf('%s/git/mrC',codeFolder)));
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib',codeFolder)));
else
end
setenv('DYLD_LIBRARY_PATH','')

%% IDENTIFY DATA LOCATION
folderNames =[];
expIdx = [];
adultExp = [1,2,3,4,5];
topFolder = '/Volumes/svndl/FinishedExperiments/2018_Kohler_NatureCommunications';
for z=1:length(adultExp) % load all adult experiments
    dataLocation(z) = subfolders(sprintf('%s/exp%0.0f*',topFolder,adultExp(z)),1);
    folderNames=[folderNames;subfolders(sprintf('%s/*20*',dataLocation{z}),1)];
    expIdx=[expIdx; adultExp(z) * ones(size(subfolders(sprintf('%s/*20*',dataLocation{z}))))];
end
% go deeper into the folder names
for f = 1:length(folderNames)
    tempFolders = subfolders(folderNames{f},1);
    tempFolders = tempFolders(cellfun(@(x) isempty(strfind(x,'mff')),tempFolders));
    folderNames{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',tempFolders{end});
end

%% SETUP INPUTS
binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= 1:4; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
condsToUse = 1:8; % if you want to include all conditions, create a vector here listing all condition numbers
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=7; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
nFreq = length(freqsToUse);
nCond = length(condsToUse);
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false; % generate source data for first instance of rca?

%% RUN RCA
for e = 1:length(adultExp)
    doExp = adultExp(e);
    fprintf('\n ... running exp %0.0f ...\n',doExp);
    saveLocation = sprintf('%s/figures/exp%d',topFolder,doExp);
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation);
    else
    end
    subsToUse = ismember(expIdx,doExp);
    
    % do RCA on all DATA
    if doExp == 1
        % only run allRCA if doing Experiment 1
        allRCA = rcaSweep(folderNames,binsToUse,freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
        allLocation = saveLocation;
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/allRCA_cov.pdf',saveLocation),'-pdf','-transparent',rcaH);
        close all;
    else
    end
    
    % do full RCA
    fullRCA = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,false);
    rcaH = grabCovFig(gcf);
    export_fig(sprintf('%s/fullRCA_cov.pdf',saveLocation),'-pdf','-transparent',rcaH);
    close all;
    
    for f=1:length(freqsToUse)
        readyRCA(f) = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle);
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/freq%dRCA_cov.pdf',saveLocation,f),'-pdf','-transparent',rcaH);
        close all;
    end
    % absord fullRCA and allRCA in readyRCA
    readyRCA = [readyRCA,fullRCA];
    clear fullRCA;
    
    % select subjects from current experiment
    for f = 1:length(allRCA);
        fIdx = length(freqsToUse)*2+f;
        readyRCA(fIdx) = allRCA(f);
        readyRCA(fIdx).data = allRCA(f).data(:,subsToUse);
        readyRCA(fIdx).comparisonData = allRCA(f).comparisonData(:,subsToUse);
        readyRCA(fIdx).inputData = allRCA(f).inputData(:,subsToUse);
        readyRCA(fIdx).noiseData.lowerSideBand = allRCA(f).noiseData.lowerSideBand(:,subsToUse);
        readyRCA(fIdx).noiseData.higherSideBand = allRCA(f).noiseData.higherSideBand(:,subsToUse);
        readyRCA(fIdx).comparisonNoiseData.lowerSideBand = allRCA(f).comparisonNoiseData.lowerSideBand(:,subsToUse);
        readyRCA(fIdx).comparisonNoiseData.higherSideBand = allRCA(f).comparisonNoiseData.higherSideBand(:,subsToUse);
    end
    
    zeroData = cellfun(@(x) any(x(:)==0), cat(1,readyRCA(:).data));
    
    if any(zeroData(:))
        error('data values exactly zero, this should not happen');
    else
    end
    
    %% COMPUTE VALUES FOR PLOTTING
    
    keepConditions = true;
    errorType = 'SEM';
    trialError = false;
    doNR = false(4,8,8); % 4 freqs, 8 RCs (with comparison), 8 conditions
    doNR([1,2,4],1,:) = true; % do fitting for first RC, first, second and fourth harmonic, all conditions
    
    for f = 1:length(readyRCA)
        fprintf('\n ... rc no. %0.0f ...\n',f);
        rcStruct = aggregateData(readyRCA(f),keepConditions,errorType,trialError,doNR);
        % RC
        readyRCA(f).stats.Amp = squeeze(rcStruct.ampBins);
        readyRCA(f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp);
        readyRCA(f).stats.ErrLB = squeeze(rcStruct.ampErrBins(:,:,:,:,1));
        readyRCA(f).stats.ErrUB = squeeze(rcStruct.ampErrBins(:,:,:,:,2));
        readyRCA(f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins);
        readyRCA(f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise);
        % Naka-Rushton
        readyRCA(f).stats.NR_Params = squeeze(rcStruct.NakaRushton.Params);
        readyRCA(f).stats.NR_R2 = squeeze(rcStruct.NakaRushton.R2);
        readyRCA(f).stats.NR_JKSE = squeeze(rcStruct.NakaRushton.JackKnife.SE);
        readyRCA(f).stats.NR_JKParams = squeeze(rcStruct.NakaRushton.JackKnife.Params);
        readyRCA(f).stats.hModel = rcStruct.NakaRushton.hModel;
        % t-values
        readyRCA(f).stats.tSqrdP = squeeze(rcStruct.tSqrdP);
        readyRCA(f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig);
        readyRCA(f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal);
    end
    % shut down parallel pool, which was used for fitting Naka-Rushton
    delete(gcp('nocreate'));
    clc;
    
    %% SAVE THE DATA, SEPARATELY FOR EACH EXPERIMENT
    saveFileName = sprintf('%s/rcaData',saveLocation);
    save(saveFileName,'readyRCA')
    clear readyRCA;
end
