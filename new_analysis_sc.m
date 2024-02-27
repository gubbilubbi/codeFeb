%% TODO

% ---Done---
% Add code for surrogate data √

% ---TODO---
% Remove the unnecessary old variables such as capoutliers and boolFC
% adapt the folder path accordingly

%% 0. Settings
clear all

set(0,'DefaultFigureVisible','on')

rng("default")

fs=200; % Sampling freq.

% The row indicies of the unordered data
order=[19,119,11,111,3,103,22,122,10,110,2,102,21,121,17,117,12,112,...
    14,114,13,113,5,105,20,120,7,107,18,118,9,109,8,108,1,101,6,106,...
    16,116,15,115,4,104];
%Order so that first 22 rows is left hemisphere, last 22 right hemisphere
[orderVals, orderInds]=sort(order); 

% Visual: blue, Sensorimotor: red, Auditory: yellow, Temporal: green, 
% Posterior: cyan, Anterior: magenta
areaColors=["b";"r";"y";"g";"c";"m";"b";"r";"y";"g";"c";"m"];

% Each ROI belongs to a region
% Taken from https://bitbucket.org/dpat/tools/raw/master/REF/ATLASES/HCP-MMP1_cortices.txt
% HCP_MMP1_cortices.txt
hemisphereColors=[repmat("b",5,1);repmat("r",4,1);repmat("y",3,1);...
               repmat("g",2,1);repmat("c",4,1);repmat("m",4,1)];
cortexColors=[hemisphereColors;hemisphereColors];

% TODO: I can probably comment out all of these later if not used again
boolPrec=0; % Precision matrix for FC
boolGroupFC=1; %Whether each group or scan has one FC
boolParCorr=0; % Partial correlation for FC
FCtoAdj="abs"; % Value to determine FC to Adjacency conversion method
methodFC="corr";
boolFDR=0; %Apply FDR on the FC matrix and the correlation values
capOutliers=0; % To cap outliers in the data to a value (4) if zscored
normMethod="zscore";
% Start and end of data timepoints to use
startInd=1;
endInd=94000;

%%% Important %%%
methodSC="normLaplacian"; % Laplacian or normLaplacian

%%
boolSavePlots=1;

%IDs for PD and HC scans
idsPD=[1,2,3,5,7,8,9,10,12,13,14,20,27,30,32,36];
idsHC=[4,6,11,15,16,17,18,19,21,22,23,24,25,26,28,29,31,33,34,35];
%idsPD=[1,3];
%idsHC=[2,4];
nPD=length(idsPD);
nHC=length(idsHC);

%Strings for plotting
plotLegends=["PD ses1","PD ses2","HC ses1","HC ses2"];
% PD first column (1,3), HC second column (2,4) in subplots
subplotsInds={1,3,2,4};

% Loading required data
L=endInd-startInd+1;

startSes1HC=nPD*2+1; % First HC after all PD
endSes1HC=nPD*2+nHC; % HC ses1 ending in data
startSes2HC=endSes1HC+1;
endSes2HC=endSes1HC+nHC;
% Start and end indicies for PD and HC in the dataset for ses1 and ses2
indsPDHC={1:nPD,nPD+1:nPD*2,startSes1HC:endSes1HC,startSes2HC:endSes2HC};
nGroups=length(indsPDHC); % 4 different groups in total

folderName=normMethod+"_SC="+methodSC+...
    "_capOutliers="+capOutliers+"_start="+startInd+"_end="+endInd;
dataFolder="./gsp_data/dataFeb/"+folderName;
%% Load data
FC=load(dataFolder+"/FC.mat").FC;
SC=load(dataFolder+"/SC.mat").SC;

nROI=size(SC,1);

% Number of scans (i.e. 2*nPatients)
% Each patient is scanned twice
nScans=size(SC,3); 
%% Calculate Laplacian

A=SC;

degreeMat=zeros(nROI,nROI,nScans);
Anorm=zeros(nROI,nROI,nScans); % Cleared afterwards
LaplacianMat=zeros(nROI,nROI,nScans); % Cleared afterwards
U=zeros(nROI,nROI,nScans);
LambdaL=zeros(nROI,nScans);


for s=1:nScans
    degreeMat(:,:,s)=diag(sum(A(:,:,s),2));
    switch methodSC
        case "Laplacian"
            LaplacianMat(:,:,s)=degreeMat(:,:,s)-A(:,:,s);
        case "normLaplacian"
            % Normalized Adjacency Matrix 
            Anorm(:,:,s)=degreeMat(:,:,s)^(-1/2)*A(:,:,s)*degreeMat(:,:,s)^(-1/2); 
            % STRUCTURAL CONNECTOME DECOMPOSITION INTO STRUCTURAL HARMONICS
            % Compute symmatrically normalized Laplacian
            LaplacianMat(:,:,s)=eye(nROI)-Anorm(:,:,s);
    end
    
    % Laplacian eigendecomposition (eigenvectors and values)
    [tempU,tempLambdaL] = eig(LaplacianMat(:,:,s));
    % eigenvalues sorted ascendingly = spatial frequencies
    [LambdaL(:,s), Ind]=sort(diag(tempLambdaL));
    % eigenvectors = structural harmonics
    U(:,:,s)=tempU(:,Ind);

    clear tempU; clear tempLambdaL; 
    clear Ind; clear Anorm; clear LaplacianMat;
end

%% Plot degree distribution
% Plot the degree disitribution 
% If it is heterogenous, then perhaps it is better to use normLaplacian
for i=1:nGroups
    figure;
    plot(mean(sum(A(:,:,indsPDHC{i}),2),3))
end

%% Spider plots of exponential function fitted to Autocorrelation Function 
fit_to_acf

%% Cutoff Selection
%%% Select one of the below options
%%% Manual, 50 or 33
%cutoffChoice='Manual';
cutoffChoice='50';
%cutoffChoice='33';

switch cutoffChoice
    case 'Manual'
        nFreqGroups=3;
        freqGroups=[1,2,3]; % low, medium and high

        lower=12;
        upper=33;
        
        % Manually selected cutoffs for each group
        K=[lower,lower,lower,lower;upper,upper,upper,upper];
        

        saveFolder=dataFolder+"/cutoff_"+lower+"_"+upper;

    case '50' %Power-spectrum density split 50/50 between low and high
        nFreqGroups=2;
        freqGroups=[1,3]; % low and high

        psdFolder=dir(dataFolder+"/psd50_*").name;
        psdFolderSplit=split(psdFolder,'_');
        
        % Cutoffs for 50/50 split
        lowCutoffs=[str2double(psdFolderSplit{2}),str2double(psdFolderSplit{3}),...
            str2double(psdFolderSplit{4}),str2double(psdFolderSplit{5})];
        K=[lowCutoffs;lowCutoffs+1];
        
        saveFolder=dataFolder+"/psd50_"+K(1,1)+"_"+K(1,2)+"_"+K(1,3)+"_"+K(1,4);

        

    case '33' %Power-spectrum density split 33/33/33 between l/m/h
        nFreqGroups=3;
        freqGroups=[1,2,3]; % low, medium and high

        psdFolder= dir(dataFolder+"/psd33_*").name;
        psdFolderSplit=split(psdFolder,'_');

        lowCutoffs=[str2double(psdFolderSplit{2}),str2double(psdFolderSplit{4}),...
            str2double(psdFolderSplit{6}),str2double(psdFolderSplit{8})];
        highCutoffs=[str2double(psdFolderSplit{3}),str2double(psdFolderSplit{5}),...
            str2double(psdFolderSplit{7}),str2double(psdFolderSplit{9})];

        K=[lowCutoffs;highCutoffs];
        
        saveFolder=dataFolder+"/psd33_"+K(1,1)+"_"+K(2,1)+"_"+K(1,2)+"_"+K(2,2)+"_"+K(1,3)+"_"+K(2,3)+"_"+K(1,4)+"_"+K(2,4);
        

    otherwise
        error('Incorrect cutoff choice')
end

%% 
% Cutoffs for freqs and groups (fgK frequency group cutoff value)
% First row low-, second row for medium- and last row for high- freq.s
fgK={1:K(1,1),1:K(1,2),1:K(1,3),1:K(1,4);...
    K(1,1)+1:K(2,1)-1,...
    K(1,2)+1:K(2,2)-1,...
    K(1,3)+1:K(2,3)-1,...
    K(1,4)+1:K(2,4)-1;...
    K(2,1):nROI,K(2,2):nROI,K(2,3):nROI,K(2,4):nROI};

% Split U into low, medium and high freqs (Laplacian eigenvectors)
% and reconstruct the MEG signal into low/medium/high signals
Ufreqs=zeros([size(U),nFreqGroups]); %Fourth dimension for low,medium,high
for i=1:nFreqGroups
    for g=1:nGroups
        Ufreqs(:,fgK{freqGroups(i),g},indsPDHC{g},i)=U(:,fgK{freqGroups(i),g},indsPDHC{g});
    end
end

% Group mean and std of the abs. value of the graph signals (filtered)
% nROI x nGroups x nFreqs (can be used for plotting)
meanX=load(saveFolder+"/meanX.mat").meanX; 
stdX=load(saveFolder+"/stdX.mat").stdX;
% nROI x nScans x nFreqs with the l2-norm
energyX=load(saveFolder+"/energyX.mat").energyX;

%% Histograms for single scans PSD cutoff choice

psdNN_lh=load(dataFolder+"/psd_single/"+"/psdNN_lh.mat",'psdNN_lh').psdNN_lh;
figure;
histogram(psdNN_lh(indsPDHC{1}))
figure;
histogram(psdNN_lh(indsPDHC{2}))
figure;
histogram(psdNN_lh(indsPDHC{3}))
figure;
histogram(psdNN_lh(indsPDHC{4}))
psdNN_lmh=load(dataFolder+"/psd_single/"+"/psdNN_lmh.mat",'psdNN_lmh').psdNN_lmh;
figure;
histogram(psdNN_lmh(indsPDHC{1}))
figure;
histogram(psdNN_lmh(indsPDHC{2}))
figure;
histogram(psdNN_lmh(indsPDHC{3}))
figure;
histogram(psdNN_lmh(indsPDHC{4}))



%% Display figures
    % The selected connectivity matrix (Adjacency or FC)
    % Distribution of eigenvalues
    % Mean eigenvectors' values for each group
    % Individual and group mean 
        % eigvals
        % zero-crossings
        % weighted zero-crossings

plots_sc
%% Codebook for the brain (atlas)
% Run if the codebook is not created yet
create_codebook 
%% Create brain plots of the eigenvectors
brain_plots

%% Plot brain signals
% meanX and stdX

%% Plot energy concentrations per regions
energy_conc_plots

%% Statistical analysis for eigenvalues between different groups
testType="mwu"; % t, mwu or ks for t-test, Mann-Whitney U-test or KS-test
[eigTestRejectNull,eigTestPvals]=stat_tests(LambdaL,nROI,1,indsPDHC,testType,saveFolder,boolSavePlots,"eigVals");

%% Statistical analysis for the energy between different groups
testType="mwu"; % t, mwu or ks for t-test, Mann-Whitney U-test or KS-test
[energyTestRejectNull,energyTestPvals]=stat_tests(energyX,nROI,nFreqGroups,indsPDHC,testType,saveFolder,boolSavePlots,"energy");

%% Brain plots for statistical significance
%% TODO GO OVER THESE 
%% Plot statistical tests on the brain

switch nROI
    case 44
        disp("44 regions")
        CodeBookpath=which('myCodeBook.mat');
    case 360
        disp("360 regions")
        CodeBookpath=which('Glasser360_2mm_codebook.mat');
    otherwise
        error(['No available atlas for the number of regions being: ' nROI])
end

CodeBook=load(CodeBookpath);
CodeBook=CodeBook.myCodeBook;

rejectMat=double(cell2mat(testRejectNull));
rejectPMat=double(cell2mat(testPvals));

%for t=1:6
    f=figure;
    t=3;
    plot_on_brain([rejectMat(:,1,t),rejectMat(:,2,t),rejectMat(:,3,t)],0,CodeBook,[0,1], ...
            ["Low","Medium","High"],0,0.001,2,0.2);
    hold on
    disp("Done 1")
    plot_on_brain([rejectMat(:,1,t),rejectMat(:,2,t),rejectMat(:,3,t)],0,CodeBook,[0,1], ...
            ["Low","Medium","High"],0,0.001,2,0);
    hold on
    disp("Done 2")
    plot_on_brain([rejectPMat(:,1,t),rejectPMat(:,2,t),rejectPMat(:,3,t)],0,CodeBook,[0,1], ...
            ["Low","Medium","High"],0,0.001,1.5,0);
%clf(f);
%end
disp("Done 3")

%% Calculate SDI for the empirical data
showSDIPlots=0;
calc_SDI

%% Perform the NeuroSynth meta-analysis
% Determine path
findpath=what('gsp_data');
%findpath=findpath(1:end-18);
mypath=findpath.path; %set code folder path
saveFolderChars=convertStringsToChars(saveFolder);
mypath=strcat(mypath,saveFolderChars(11:end));

% Do it once per group..?
for g=1:nGroups
    disp(g)
    neurosynth_inputs
end
