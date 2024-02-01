%% TODO

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
boolGroupFC=1; %Whether each group or subject has one FC
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

nFreqGroups=3;

%IDs for PD and HC subjects
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
nSubjs=size(SC,3);
%% Calculate Laplacian

A=SC;

DegreeMat=zeros(nROI,nROI,nSubjs); % Cleared afterwards
Anorm=zeros(nROI,nROI,nSubjs); % Cleared afterwards
LaplacianMat=zeros(nROI,nROI,nSubjs); % Cleared afterwards
U=zeros(nROI,nROI,nSubjs);
LambdaL=zeros(nROI,nSubjs);

for s=1:nSubjs
    DegreeMat(:,:,s)=diag(sum(A(:,:,s),2));
    switch methodSC
        case "Laplacian"
            LaplacianMat(:,:,s)=DegreeMat(:,:,s)-A(:,:,s);
        case "normLaplacian"
            % Normalized Adjacency Matrix 
            Anorm(:,:,s)=DegreeMat(:,:,s)^(-1/2)*A(:,:,s)*DegreeMat(:,:,s)^(-1/2); 
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
    clear Ind; clear DegreeMat; clear Anorm; clear LaplacianMat;
end

%% Spider plots of exponential function fitted to Autocorrelation Function 
fit_to_acf

%% Find best cutoff for L/M/H
% As there is no easy way to determine where low freqs end and medium
% freqs start etc. we here try to find cutoff values for low and high freqs
% such that the results are similiar to the results if the cutoff values
% were chosen slightly different (+- bound). As a change in eig.vec.s
% change the GFT of the signal and thus the results, a robust parameter
% choice at least indicates that the results are somewhat stable for our
% choice of parameters and that the choices of cutoffs therefore lead to
% some sort of "low/medium/high" freqs instead of something more arbitrary.

% Manually set value for in which range the correlation should be
% calculated from a given cutoff value
bound=2;  
robustParamsAll

%% Plot figures for mean correlation values
plotRobustParamsAll

%% Values for selected parameters

%%% Settings %%%
% Change i and j to get a different cutoff freq for low and high
i=10; % Low cutoff idnex
j=1;  % High cutoff index
% Lowest acceptable correlation (boundary correlation value)
bCorr=0.90;
% Groups to display
dGroups=[1,2,3,4]; %all groups: [1,2,3,4] or [1,3] for different SCs
%%% Settings %%%

k=nHighs*(i-1)+j;
lower=lows(i);
upper=highs(j);
disp("-------")
disp("Low cutoff: "+lower);
disp("High cutoff: "+upper);
disp("Medium index: "+k);
disp("-------")

% Get the correlation values
corrChosenParams

%% Load data 
%% TODO: Change this so that we can have different cutoffs for each group?
% Load selected data based on picked cutoff values from robustParams
    % Need empty array for meanX, stdX and energyX
    % Need a new savefolder 

% Needed if I want to load in different cutoff data for each group
%meanX=zeros(nROI,nSubjs,nFreqGroups);
%stdX=zeros(nROI,nSubjs,nFreqGroups);
%energyX=zeros(nROI,nSubjs,nFreqGroups);

%%% Settings %%%
lower=10;
upper=35;

% Manually selected cutoffs for each group
NN=[lower,lower,lower,lower;upper,upper,upper,upper];

% Cutoffs for freqs and groups
% First row low-, second row for medium- and last row for high- freq.s
fgNN={1:NN(1,1),1:NN(1,2),1:NN(1,3),1:NN(1,4);...
    NN(1,1)+1:NN(2,1)-1,...
    NN(1,2)+1:NN(2,2)-1,...
    NN(1,3)+1:NN(2,3)-1,...
    NN(1,4)+1:NN(2,4)-1;...
    NN(2,1):nROI,NN(2,2):nROI,NN(2,3):nROI,NN(2,4):nROI};

% Split U into low, medium and high freqs (Laplacian eigenvectors)
% and reconstruct the MEG signal into low/medium/high signals
Ufreqs=zeros([size(U),nFreqGroups]); %Third dimension for low,medium,high
for i=1:nFreqGroups
    for g=1:nGroups
        Ufreqs(:,fgNN{i,g},indsPDHC{g},i)=U(:,fgNN{i,g},indsPDHC{g});
    end
end

saveFolder=dataFolder+"/cutoff_"+lower+"_"+upper;
meanX=load(saveFolder+"/meanX.mat").meanX;
stdX=load(saveFolder+"/stdX.mat").stdX;
energyX=load(saveFolder+"/energyX.mat").energyX;
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
%% Create brain plots
brain_plots

%% Plot brain signals


%% Investigate energy and statistical analysis between different groups


%% Null Models


%% Calculate SDI


