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

clear all
close all

mainPath = '/Volumes/svndl/FinishedExperiments/2018_Kohler_NatureCommunications';
figureFolder = sprintf('%s/figures/infant_exp',mainPath);

%% IDENTIFY DATA LOCATION
folderNames=[];
expIdx=[];
dataLocation = '/Volumes/svndl/FinishedExperiments/2018_Kohler_NatureCommunications/infant_exp';
folderNames = [folderNames; subfolders(sprintf('%s/*20*',dataLocation),1)]; %change back to 128 between wild cards
expIdx = [expIdx; ones(size(subfolders(sprintf('%s/*20*',dataLocation))))];

for f = 1:length(folderNames)
    tempFolders = subfolders(folderNames{f},1);
    tempFolders = tempFolders(cellfun(@(x) ~isempty(strfind(x,'Exp')),tempFolders));
    folderNames{f} = sprintf('%s',tempFolders{end});
end

%%  SET UP INPUTS
binsToUse=1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse= [1,2]; % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)
trialsToUse = []; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
condsToUse = 1:4; % if you want to include all conditions, create a vector here listing all condition numbers
condSep = [3,7];
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=7; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
nFreq = length(freqsToUse);
nCond = length(condsToUse);
chanToCompare = 75; % channel to use for a performance evaluation, can be []
dataType = 'RLS'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false; % generate source data for first instance of rca?
keepConditions = true;
errorType = 'SEM';
doNR = false(4,8,8); % 4 freqs, 8 RCs (with comparison), 8 conditions
doNR([1,2,4],[1,5],:) = true; % do fitting for first and fifth RC, first, second and fourth harmonic, all conditions
trialError = false;

%%  RUN RCA
if ~exist(figureFolder,'dir')
    mkdir(figureFolder);
else
end

%do RCA on all DATA
babyRCA((1:length(freqsToUse))+length(freqsToUse)) = rcaSweep(folderNames, binsToUse, freqsToUse,condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
rcaH = grabCovFig(gcf);
export_fig(sprintf('%s/BabyCombined%dF%dF_cov.pdf',figureFolder,min(freqsToUse),max(freqsToUse)),'-pdf','-transparent',rcaH);
close all;

for f=1:length(freqsToUse)
    babyRCA(f) = rcaSweep(folderNames,binsToUse,freqsToUse(f),condsToUse,trialsToUse,nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,false);
    rcaH = grabCovFig(gcf);
    export_fig(sprintf('%s/BabyCombined%dF_cov.pdf',figureFolder,freqsToUse(f)),'-pdf','-transparent',rcaH);
    close all;
end

zeroData = cellfun(@(x) any(x(:)==0), cat(1,babyRCA(:).data));
    
if any(zeroData(:))
    error('data values exactly zero, this should not happen');
else
end
    
%     NB REMOVING NANS NO LONGER NECESSARY!
%     nanDims = [1,2]; % if all time points are zero, or all channels are zero
%     structVars = {'data','noiseData','comparisonData','comparisonNoiseData'};
%     noiseVars = {'lowerSideBand','higherSideBand'};
%     
%     for z=1:length(structVars)
%         if strfind(lower(structVars{z}),'noise')
%             for n = 1:length(noiseVars)
%                 allRCA.(structVars{z}).(noiseVars{n}) = allRCA.(structVars{z}).(noiseVars{n})(:,subsToUse);
%                 allRCA.(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}).(noiseVars{n}),'uni',false);
%                 for f=1:length(freqsToUse)
%                     freqRCA(f).(structVars{z}).(noiseVars{n}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}).(noiseVars{n}),'uni',false);
%                 end
%             end
%         else
%             % populate allRCA
%             allRCA.(structVars{z}) = allRCA.(structVars{z})(:,subsToUse);
%             allRCA.(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),allRCA.(structVars{z}),'uni',false);
%             for f=1:length(freqsToUse)
%                 freqRCA(f).(structVars{z}) = cellfun(@(x) Zero2NaN(x,nanDims),freqRCA(f).(structVars{z}),'uni',false);
%             end
%         end
%     end

%% COMPUTE VALUES FOR PLOTTING
for f = 1:length(babyRCA)
    rcStruct = aggregateData(babyRCA(f),keepConditions,errorType,trialError,doNR);
    % RC
    babyRCA(f).stats.Amp = squeeze(rcStruct.ampBins);
    babyRCA(f).stats.SubjectAmp = squeeze(rcStruct.subjectAmp);
    babyRCA(f).stats.ErrLB = squeeze(rcStruct.ampErrBins(:,:,:,:,1));
    babyRCA(f).stats.ErrUB = squeeze(rcStruct.ampErrBins(:,:,:,:,2));
    babyRCA(f).stats.NoiseAmp = squeeze(rcStruct.ampNoiseBins);
    babyRCA(f).stats.SubjectNoiseAmp = squeeze(rcStruct.subjectAmpNoise);
    % Naka-Rushton
    babyRCA(f).stats.NR_Params = squeeze(rcStruct.NakaRushton.Params);
    babyRCA(f).stats.NR_R2 = squeeze(rcStruct.NakaRushton.R2);
    babyRCA(f).stats.NR_JKSE = squeeze(rcStruct.NakaRushton.JackKnife.SE);
    babyRCA(f).stats.NR_JKParams = squeeze(rcStruct.NakaRushton.JackKnife.Params);
    babyRCA(f).stats.hModel = rcStruct.NakaRushton.hModel;
    % t-values
    babyRCA(f).stats.tSqrdP = squeeze(rcStruct.tSqrdP);
    babyRCA(f).stats.tSqrdSig = squeeze(rcStruct.tSqrdSig);
    babyRCA(f).stats.tSqrdVal = squeeze(rcStruct.tSqrdVal);
end
clear rcStruct;

% shut down parallel pool, which was used for fitting Naka-Rushton
delete(gcp('nocreate'));

save(sprintf('%s/BabyDataOutput.mat', figureFolder),'babyRCA');