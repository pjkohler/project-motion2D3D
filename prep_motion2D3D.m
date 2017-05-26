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
for z=1:3 % load all experiments
    switch z
        case 1
            % experiment 1: original
            dataLocation{z} = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/2D3DHV_SQWAVE_20bar1meter/';
        case 2
            % experiment 2: blank ref
            dataLocation{z} = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/BLANKREF_2D3DHV_SQWAVE_20bar1meter';
        case 3
            % experiment 3: near far
            dataLocation{z} = '/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/NEARFAR_2D3DHV20bar1meter';
    end
    folderNames=[folderNames;subfolders(sprintf('%s/*20*',dataLocation{z}),1)];
    expIdx=[expIdx; z * ones(size(subfolders(sprintf('%s/*20*',dataLocation{z}))))];
end
% go deeper into the folder names
for f = 1:length(folderNames)
    tempFolders = subfolders(folderNames{f},1);
    folderNames{f} = sprintf('%s/Exp_TEXT_HCN_128_Avg',tempFolders{end});
end

%% SETUP INPUTS
binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= 1:4; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
condsToUse = 1:8; % if you want to include all conditions, create a vector here listing all condition numbers
condSep = [3,7];
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
nFreq = length(freqsToUse);
nCond = length(condsToUse);
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false; % generate source data for first instance of rca?

%% RUN RCA
for doExp = 1:3
    saveLocation = sprintf('/Volumes/Denali_4D2/kohler/EEG_EXP/DATA/motion2D3D/figures/exp%d',doExp);
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
    else
        % otherwise, just load the data
        load(sprintf('%s/rcaData.mat',allLocation));
    end
    
    % do full RCA
    fullRCA = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,false);
    
    % do full RCA, on horizontal and vertical relative disparity only
    for c = 1:length(condSep)
        condRCA(c) = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse(1:2:end),condSep(c),trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,false);
    end
    rcaH = grabCovFig(gcf);
    export_fig(sprintf('%s/fullRCA_cov.pdf',saveLocation),'-pdf','-transparent',rcaH);
    close all;
    for f=1:length(freqsToUse)
        freqRCA(f) = rcaSweep(folderNames(subsToUse),binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle);
        rcaH = grabCovFig(gcf);
        export_fig(sprintf('%s/freq%dRCA_cov.pdf',saveLocation,f),'-pdf','-transparent',rcaH);
        close all;
    end

    % rca replaces NaNs with zeroes, correct this
    allRCA.inputData = allRCA.inputData(:,subsToUse);
    nanDims = [1,2]; % if all time points are zero, or all channels are zero
    structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
    noiseVars = {'lowerSideBand','higherSideBand'};
    for z=1:length(structVars)
        if strfind(lower(structVars{z}),'noise')
            for n = 1:length(noiseVars)
                allRCA.(structVars{z}).(noiseVars{n}) = allRCA.(structVars{z}).(noiseVars{n})(:,subsToUse);
                allRCA.(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}).(noiseVars{n}),'uni',false);
                fullRCA.(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),fullRCA.(structVars{z}).(noiseVars{n}),'uni',false);
                for f=1:length(freqsToUse)
                    freqRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
                end
                for c = 1:length(condSep)
                    condRCA(c).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),condRCA(c).(structVars{z}).(noiseVars{n}),'uni',false);
                end
            end
        else
            % populate allRCA
            allRCA.(structVars{z}) = allRCA.(structVars{z})(:,subsToUse);
            allRCA.(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}),'uni',false);
            fullRCA.(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),fullRCA.(structVars{z}),'uni',false);
            for f=1:length(freqsToUse)
                freqRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}),'uni',false);
            end
            for c = 1:length(condSep)
                condRCA(c).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),condRCA(c).(structVars{z}),'uni',false);
            end
        end
    end
    % now save the data, seperately for each experiment
    saveFileName = sprintf('%s/rcaData',saveLocation);
    if ~exist([saveFileName,'.mat'], 'file')
        save(saveFileName,'*ToUse') % create file if it does not exist
    else
        save(saveFileName,'*ToUse','-append') % otherwise append
    end
    save(saveFileName,'allRCA','-append');
    save(saveFileName,'fullRCA','condRCA','-append')
    save(saveFileName,'freqRCA','-append')
end
