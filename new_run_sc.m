
%% MAIN 2022

% ---Done---
% Fix FC parameters for connectivity matrix? Or can I skip? √
% Remove Adjacency calc for FC √
% Don't calculate create_corr_timeseries √
% Check calc_acf √
% Save the folder with methodSC and the data in it √
% Calculate lmh (which boundaries?) √
% Calculate energy 50/50 split and then GFT √
% Save the mPSD √
% Find the plotting of the acf √
% Add code for surrogate data (null models) √
    % informed SC √
    % ignorant SC √
    % Save necessary things √
% Implement 33/33/33 energy split √

% ---TODO---
% Fixed plotting for PSD??
% Load SC
% Remove all the code that is not run at KI


% ---Post data retrieval---
% Calculate SDI and save it √


%%  Graph Signal Processing (GSP) for Neuroimaging - Giulia Preti 17.12.2022

% Author: Valter Lundegårdh

% Inspiration from Giulia Preti's Graph Signal Processing
% https://github.com/gpreti/GSP_StructuralDecouplingIndex/tree/master
% DOI: 10.1109/ICASSP49357.2023.10095285

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
[~, orderInds]=sort(order); 

% Visual: blue, Sensorimotor: red, Auditory: yellow, Temporal: green, 
% Posterior: cyan, Anterior: magenta
areaColors=["b";"r";"y";"g";"c";"m";"b";"r";"y";"g";"c";"m"];

% Each ROI belongs to a region
hemisphereColors=[repmat("b",5,1);repmat("r",4,1);repmat("y",3,1);...
               repmat("g",2,1);repmat("c",4,1);repmat("m",4,1)];
cortexColors=[hemisphereColors;hemisphereColors];

plotTitles=["PD1","PD2","HC1","HC2"];

%% Settings
boolPrec=0; % Precision matrix for FC
boolGroupFC=1; %Whether each group or scan has one FC
boolParCorr=0; % Partial correlation for FC
FCtoAdj="abs"; % Value to determine FC to Adjacency conversion method

methodFC="corr";
boolFDR=0; %Apply FDR on the FC matrix and the correlation values

capOutliers=0; % To cap outliers in the data to a value (4) if zscored
normMethod="zscore";

%%% Change these %%%
dataFolderName="dataFeb";
methodSC="normLaplacian"; % Laplacian or normLaplacian

%%
boolSavePlots=0;

% Start and end of data timepoints to use
startInd=1;
endInd=94000;

%IDs for PD and HC scans
idsPD=[1,2,3,5,7,8,9,10,12,13,14,20,27,30,32,36];
idsHC=[4,6,11,15,16,17,18,19,21,22,23,24,25,26,28,29,31,33,34,35];
%idsPD=[1,3];
%idsHC=[2,4];
nPD=length(idsPD);
nHC=length(idsHC);

%Strings for plotting
plotLegends=["PD ses1","PD ses2","HC ses1","HC ses2"];
subplotsInds={1,3,2,4};

% 1. Loading required data
dataLength=endInd-startInd+1;

data=load_data('/archive/21098_pd_dimensionality/MEG/',44,startInd:endInd,idsHC);% time x regions x scans
%data=load_data('scilife/MEG/',44,startInd:endInd,idsHC);% time x regions x scans

% Reorder data so in ascending order based on ID
data=data(:,orderInds,:);

nScans=size(data,3); % Number of scans = nPatients*2 (each scanned twice)
nROI=size(data,2); %44

% Start and end indicies for PD and HC in the dataset for ses1 and ses2
startSes1HC=nPD*2+1;
endSes1HC=nPD*2+nHC;
startSes2HC=endSes1HC+1;
endSes2HC=endSes1HC+nHC;
indsPDHC={1:nPD,nPD+1:nPD*2,startSes1HC:endSes1HC,startSes2HC:endSes2HC};
nGroups=length(indsPDHC); % 4 different groups in total

% z-score the data: Normalize each column (region) to mean 0 & std 1
data=zscore(data); 
disp("zscore normalisation")

if(capOutliers)
    %thr=3; %0.00135 => ~1% of data for every region
    thr=4; %0.00003 => ~0.3% of data for every region
    data(data>thr)=thr;
    data(data<-thr)=-thr;
end

% Create folder for psd figures, autocorrelation function values
% and the codebook used for plotting
saveMiscFolder="./gsp_data/"+dataFolderName+"/misc/";
if ~exist(saveMiscFolder, 'dir')
     mkdir(saveMiscFolder)
end
%% %% ONLY RUN THIS ONCE PER
%saveDataFolder="./scilife/testData/"+folderName;
if ~exist(saveMiscFolder, 'dir')
       mkdir(saveMiscFolder)
end
%[corrMeans,corrStds]=create_corr_timeseries(data,methodFC,indsPDHC,saveMiscFolder);

% params for exponential fitted to acf
acfVals=calc_acf(data);
save(saveMiscFolder+"/acfVals.mat","acfVals")

%% Calculate the FC matrix

% Calculate FC matrix for each scan and mean FC for an entire group
[FC,meanFC,pvalsMat]=create_connectivity_matrix(data,methodFC,indsPDHC,boolGroupFC,FCtoAdj,boolFDR,boolParCorr,boolPrec);

%% TODO: Load the Strucutral connectivity matrices
SC=abs(meanFC);
A=SC;

for s=1:nScans
    % Check that the adjacency matrix is symmetric
    if(~issymmetric(A(:,:,s)))
        error("Adjacency not symmetric")
    end
    
    % Check that the graph is connected
    dirGraph = digraph(A(:,:,s));
    bins = conncomp(dirGraph, 'Type', 'weak');
    isConnected = all(bins == 1);
    if(~isConnected)
        error('Adjacency not connected')
    end
end

%% Calculate Laplacian

degreeMat=zeros(nROI,nROI,nScans);
Anorm=zeros(nROI,nROI,nScans);
LaplacianMat=zeros(nROI,nROI,nScans); % Cleared afterwards
U=zeros(nROI,nROI,nScans);
LambdaL=zeros(nROI,nScans);

for s=1:nScans
    degreeMat(:,:,s)=diag(sum(A(:,:,s),2));
    switch methodSC
        case "Laplacian"
            LaplacianMat(:,:,s)=degreeMat(:,:,s)-A(:,:,s);
        case "normLaplacian"
            % normalized Adjacency Matrix 
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
    clear Ind; clear LaplacianMat;
end 

%%% 2. GRAPH FOURIER TRANSFORM
X_hat=zeros(nROI,dataLength,nScans);
for s=1:nScans %graph fourier transform
    %nROIxnROI*nROIxtime
    % U Complex Transpose
    X_hat(:,:,s)=U(:,:,s)'*data(:,:,s).';
end

% Plot eigenvalues (i.e. variation of eigenvectors) and select cutoff
plot_scans_and_mean(LambdaL,"Laplacian eigenvalues/Variation of eigenvectors" ...
    ,"Spectral index","Variation",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    NaN,0)

% Calculate (weigthed) zero-crossings 
ZC=zeros(nROI,nScans);
wZC=zeros(nROI,nScans);

for s=1:nScans
    for u=1:nROI %for each eigenvector
        UU=U(:,u,s);
        sumZC=0; %ZC summ
        sumWZC=0; %wZC summ
        %for each eigvec, check all sign change between its regions.
        for i=1:nROI-1 
            for j=i+1:nROI
                if (UU(i)*UU(j))<0 %if signals are of opposite signs
                    sumWZC=sumWZC+A(i,j,s);
                    sumZC=sumZC+1;
                end
                ZC(u,s)=sumZC; % Num of zero-crossings
                wZC(u,s)=sumWZC; %wZC
            end
        end
    end
end

folderName=normMethod+"_SC="+methodSC+...
    "_capOutliers="+capOutliers+"_start="+startInd+"_end="+endInd;
saveFolder="./gsp_data/"+dataFolderName+"/"+folderName;
if ~exist(saveFolder, 'dir')
       mkdir(saveFolder)
end

% Save matrices
save(saveFolder+"/SC.mat","SC")
save(saveFolder+"/FC.mat","FC")
save(saveFolder+"/meanFC.mat","meanFC")
save(saveFolder+"/pvalsMat.mat","pvalsMat")

%% 3. SPECTRAL FILTERING OF FUNCTIONAL SIGNALS INTO LOW- AND HIGH-PASS
% Average energy spectral density of functional data projected on the 
% structural harmonics

% nROI x length x nScans
pow=abs(X_hat).^2; % Modulus of the complex number squared = Energy

% Mean over timeseries => one value for each scan in each region
PSD=squeeze(mean(pow,2)); %% nROI x nScans

% Values for plotting and calculating AUC for each group
avg=zeros(nROI,nGroups);
stdPSD=zeros(nROI,nGroups);
upper1=zeros(nROI,nGroups);
lower1=zeros(nROI,nGroups);
idx=zeros(nROI,nGroups);
mPSD=zeros(nROI,nGroups);
AUCTOT=zeros(nGroups,1);
n=zeros(nGroups,1);
psdNN=zeros(nGroups,1);

for g=1:nGroups
    % These are for optional plotting
    avg(:,g)=mean(PSD(:,indsPDHC{g}),2);
    stdPSD(:,g)=std(PSD(:,indsPDHC{g}),0,2); 
    upper1(:,g)=avg(:,g)+stdPSD(:,g);
    lower1(:,g)=avg(:,g)-stdPSD(:,g);
    idx(:,g)=max(PSD(:,indsPDHC{g}),[],2)>0 & min(PSD(:,indsPDHC{g}),[],2)>0 & mean(PSD(:,indsPDHC{g}),2)>0;
    % Mean power spectrum density over each group
    mPSD(:,g)=mean(PSD(:,indsPDHC{g}),2);

    AUCTOT(g)=trapz(mPSD(:,g)); %total area under the curve 
    AUC=0;
    n(g)=0;

    minDiff=2; %theoretical total difference can at most be 2
    for lowC=1:nROI
        lAUC=trapz(mPSD(1:lowC,g))/AUCTOT(g); %low AUC
        hAUC=trapz(mPSD(lowC+1:end,g))/AUCTOT(g); % high AUC
        totDiff=abs(lAUC-hAUC);

        % If the absolute difference is smaller then update the 
        % best cutoff values for the most even AUC split (33/33/33)
        if(totDiff < minDiff)
            psdNN(g,1)=lowC;
            minDiff=totDiff;
        end
    end

    % Plot the PSD
    %patch takes [x-coords],[y-coords],[colour]
    psdFig=figure;
    patch([LambdaL(:,indsPDHC{g}(1))', fliplr(LambdaL(:,indsPDHC{g}(1))')], ...
       [lower1(:,g)' fliplr(upper1(:,g)')], [0.8 0.8 0.8]);
    hold on;
    % One lambda per group
    plot(LambdaL(:,indsPDHC{g}(1)),mPSD(:,g));
    %xlim([0.05 2]);ylim([1 50]);
    title("Spectral Energy Density "+plotTitles(g));xlabel('Harmonic Frequency');ylabel('Energy')
    %xline(LambdaL(psdNN(g))) % Plot cutoff
    set(gca, 'XScale', 'log', 'YScale','log')
    % Save the PSD plots
    %if(boolSavePlot)
    saveas(psdFig,saveMiscFolder+"/psd_group_"+plotTitles(g)+".png")
    %end
end

% vector for possibility to perform custom cutoff for each group
strCutoffs=strjoin(reshape((string([psdNN(1),psdNN(2),psdNN(3),psdNN(4)])).',1,[]),"_");
psdFolder=saveFolder+"/psd50_"+strCutoffs+"/";
if ~exist(psdFolder, 'dir')
       mkdir(psdFolder);
end

%%% Create 2 versions of the U matrix, containing only low- or high- frequency Laplacian eigenvectors
Ulow=zeros(size(U));
Uhigh=zeros(size(U));
for g=1:nGroups
    Ulow(:,1:psdNN(g),indsPDHC{g})=U(:,1:psdNN(g),indsPDHC{g}); %low spatial frequencies  
    Uhigh(:,psdNN(g)+1:end,indsPDHC{g})=U(:,psdNN(g)+1:end,indsPDHC{g}); %high spatial frequencies
end

%%% reconstruct functional signals containing only low / high frequencies
% This both filters and inverts the data back to time domain 
% because V^-1 * H^ (H binary diag matrix) takes only the specified values
% from V^-1, i.e. low or high components.

Xc=zeros(dataLength,nROI,nScans);
Xd=zeros(dataLength,nROI,nScans);

for s=1:nScans
    % reconstruction of coupled signal, low, medium and high freqs.
    Xc(:,:,s)=(Ulow(:,:,s)*X_hat(:,:,s)).';
    Xd(:,:,s)=(Uhigh(:,:,s)*X_hat(:,:,s)).'; 
end

% Energy (l2-norm) over time for each region and scan
% For SDI and comparing regions
energyX=zeros(nROI,nScans,2);
energyX(:,:,1)=squeeze(vecnorm(Xc,2,1));
energyX(:,:,2)=squeeze(vecnorm(Xd,2,1));

% For mean and std energy plots on the brain
meanX=zeros(nROI,nGroups,2); %low,high
stdX=zeros(nROI,nGroups,2); %low,high

for g=1:nGroups
    % Mean and std of reconstructed signals over time and session group 
    meanX(:,g,1)=mean(abs(Xc(:,:,indsPDHC{g})),[1,3]);
    meanX(:,g,2)=mean(abs(Xd(:,:,indsPDHC{g})),[1,3]);
    stdX(:,g,1)=std(abs(Xc(:,:,indsPDHC{g})),0,[1,3]);
    stdX(:,g,2)=std(abs(Xd(:,:,indsPDHC{g})),0,[1,3]);
end

save(psdFolder+"/mPSD.mat",'mPSD')
save(psdFolder+"/meanX.mat",'meanX')
save(psdFolder+"/stdX.mat",'stdX')
save(psdFolder+"/energyX.mat",'energyX')

clear meanX
clear stdX
clear energyX





%% 3. SPECTRAL FILTERING OF FUNCTIONAL SIGNALS INTO LOW-, MEDIUM- and HIGH-PASS
% Average energy spectral density of functional data projected on the 
% structural harmonics

% nROI x length x nScans
pow=abs(X_hat).^2; % Modulus of the complex number squared = Energy

% Mean over timeseries => one value for each scan in each region
PSD=squeeze(mean(pow,2)); %% nROI x nScans

% Values for plotting and calculating AUC for each group
avg=zeros(nROI,nGroups);
stdPSD=zeros(nROI,nGroups);
upper1=zeros(nROI,nGroups);
lower1=zeros(nROI,nGroups);
idx=zeros(nROI,nGroups);
mPSD=zeros(nROI,nGroups);
AUCTOT=zeros(nGroups,1);
AUC=zeros(nGroups,1);
n=zeros(nGroups,1);
psdNN=zeros(nGroups,2);

for g=1:nGroups
    % These are for optional plotting
    avg(:,g)=mean(PSD(:,indsPDHC{g}),2);
    stdPSD(:,g)=std(PSD(:,indsPDHC{g}),0,2); 
    upper1(:,g)=avg(:,g)+stdPSD(:,g);
    lower1(:,g)=avg(:,g)-stdPSD(:,g);
    idx(:,g)=max(PSD(:,indsPDHC{g}),[],2)>0 & min(PSD(:,indsPDHC{g}),[],2)>0 & mean(PSD(:,indsPDHC{g}),2)>0;
    
    % Mean power spectrum density over each group
    mPSD(:,g)=mean(PSD(:,indsPDHC{g}),2);

    AUCTOT(g)=trapz(mPSD(:,g)); %total area under the curve 
    AUC(g)=0;
    n(g)=0;

    minDiff=2; %theoretical total difference can at most be 2

    for lowC=1:nROI
        for highC=lowC+1:nROI
            lAUC=trapz(mPSD(1:lowC,g))/AUCTOT(g); %low AUC
            mAUC=trapz(mPSD(lowC+1:highC-1,g))/AUCTOT(g); % medium AUC
            hAUC=trapz(mPSD(highC:end,g))/AUCTOT(g); % high AUC
            totDiff=abs(lAUC-mAUC)+abs(mAUC-hAUC)+abs(lAUC-hAUC);

            % If the absolute difference is smaller then update the 
            % best cutoff values for the most even AUC split (33/33/33)
            if(totDiff < minDiff)
                psdNN(g,1)=lowC;
                psdNN(g,2)=highC;
                minDiff=totDiff;
            end
        end
    end
    
    % splitLow=trapz(mPSD(1:psdNN(g,1),g))/AUCTOT(g);
    % splitMedium=trapz(mPSD(psdNN(g,1)+1:psdNN(g,2)-1,g))/AUCTOT(g);
    % splitHigh=trapz(mPSD(psdNN(g,2):end,g))/AUCTOT(g);
   
    % disp("Energy split "+plotTitles(g))
    % disp(splitLow)
    % disp(splitMedium)
    % disp(splitHigh)
    % disp("Total difference")
    % disp(abs(splitLow-splitMedium)+abs(splitMedium-splitHigh)+abs(splitLow-splitHigh));
    
    % Check that the TOTAUC is 1
    %disp((trapz(mPSD(1:psdNN(g,1),g))+trapz(mPSD(psdNN(g,1):psdNN(g,2),g))+trapz(mPSD(psdNN(g,2):end,g)))/AUCTOT(g))
    
end

% vector for possibility to perform custom cutoff for each group
strCutoffs=strjoin(reshape((string([psdNN(1,1),psdNN(1,2),psdNN(2,1),psdNN(2,2),psdNN(3,1),psdNN(3,2),psdNN(4,1),psdNN(4,2)])).',1,[]),"_");
psdFolder=saveFolder+"/psd33_"+strCutoffs+"/";
if ~exist(psdFolder, 'dir')
       mkdir(psdFolder);
end

%%% Create 3 versions of the U matrix, containing only low-, medium- or high- frequency Laplacian eigenvectors
Ulow=zeros(size(U));
Umedium=zeros(size(U));
Uhigh=zeros(size(U));
for g=1:nGroups
    Ulow(:,1:psdNN(g,1),indsPDHC{g})=U(:,1:psdNN(g,1),indsPDHC{g}); %low spatial frequencies  
    Umedium(:,psdNN(g,1)+1:psdNN(g,2)-1,indsPDHC{g})=U(:,psdNN(g,1)+1:psdNN(g,2)-1,indsPDHC{g}); %medium spatial frequencies
    Uhigh(:,psdNN(g,2):end,indsPDHC{g})=U(:,psdNN(g,2):end,indsPDHC{g}); %high spatial frequencies
end

%%% reconstruct functional signals containing only low / high frequencies
% This both filters and inverts the data back to time domain 
% because V^-1 * H^ (H binary diag matrix) takes only the specified values
% from V^-1, i.e. low or high components.

Xl=zeros(dataLength,nROI,nScans);
Xm=zeros(dataLength,nROI,nScans);
Xh=zeros(dataLength,nROI,nScans);

for s=1:nScans
    % reconstruction of coupled signal, low, medium and high freqs.
    Xl(:,:,s)=(Ulow(:,:,s)*X_hat(:,:,s)).';
    Xm(:,:,s)=(Umedium(:,:,s)*X_hat(:,:,s)).';
    Xh(:,:,s)=(Uhigh(:,:,s)*X_hat(:,:,s)).'; 
end

% Energy (l2-norm) over time for each region and scan
% For SDI and comparing regions
energyX=zeros(nROI,nScans,3);
energyX(:,:,1)=squeeze(vecnorm(Xl,2,1));
energyX(:,:,2)=squeeze(vecnorm(Xm,2,1));
energyX(:,:,3)=squeeze(vecnorm(Xh,2,1));

% For mean and std energy plots on the brain
meanX=zeros(nROI,nGroups,3); %low,high
stdX=zeros(nROI,nGroups,3); %low,high

for g=1:nGroups
    % Mean and std of reconstructed signals over time and session group 
    meanX(:,g,1)=mean(abs(Xl(:,:,indsPDHC{g})),[1,3]);
    meanX(:,g,2)=mean(abs(Xm(:,:,indsPDHC{g})),[1,3]);
    meanX(:,g,3)=mean(abs(Xh(:,:,indsPDHC{g})),[1,3]);
    stdX(:,g,1)=std(abs(Xl(:,:,indsPDHC{g})),0,[1,3]);
    stdX(:,g,2)=std(abs(Xm(:,:,indsPDHC{g})),0,[1,3]);
    stdX(:,g,3)=std(abs(Xh(:,:,indsPDHC{g})),0,[1,3]);
end

save(psdFolder+"/mPSD.mat",'mPSD')
save(psdFolder+"/meanX.mat",'meanX')
save(psdFolder+"/stdX.mat",'stdX')
save(psdFolder+"/energyX.mat",'energyX')

clear meanX
clear stdX
clear energyX

%% 3. SPECTRAL FILTERING OF FUNCTIONAL SIGNALS FOR SINGLE SCANS
% Low/High and Low/Medium/High

% nROI x length x nScans
pow=abs(X_hat).^2; % Modulus of the complex number squared = Energy

% Mean over timeseries => one value for each scan in each region
PSD=squeeze(mean(pow,2)); %% nROI x nScans

% Values for plotting and calculating AUC for each scan
AUCTOT=zeros(nScans,1);
psdNN_lh=zeros(nScans,1);
psdNN_lmh=zeros(nScans,2);

for s=1:nScans
    % power spectrum density for the scan
    AUCTOT(s)=trapz(PSD(:,s)); %total area under the curve 

    minDiff=1; %theoretical total difference can at most be 1 for 50/50
    for lowC=1:nROI
        lAUC=trapz(PSD(1:lowC,s))/AUCTOT(s); %low AUC
        hAUC=trapz(PSD(lowC+1:end,s))/AUCTOT(s); % high AUC
        totDiff=abs(lAUC-hAUC);

        % If the absolute difference is smaller then update the 
        % best cutoff values for the most even AUC split (33/33/33)
        if(totDiff < minDiff)
            psdNN_lh(s)=lowC;
            minDiff=totDiff;
        end
    end
  
    minDiff=2; %theoretical total difference can at most be 2 for 33/33/33

    for lowC=1:nROI
        for highC=lowC+1:nROI
            lAUC=trapz(PSD(1:lowC,s))/AUCTOT(s); %low AUC
            mAUC=trapz(PSD(lowC+1:highC-1,s))/AUCTOT(s); % medium AUC
            hAUC=trapz(PSD(highC:end,s))/AUCTOT(s); % high AUC
            totDiff=abs(lAUC-mAUC)+abs(mAUC-hAUC)+abs(lAUC-hAUC);

            % If the absolute difference is smaller then update the 
            % best cutoff values for the most even AUC split (33/33/33)
            if(totDiff < minDiff)
                psdNN_lmh(s,1)=lowC;
                psdNN_lmh(s,2)=highC;
                minDiff=totDiff;
            end
        end
    end

end

psdSingleFolder=saveFolder+"/psd_single/";
if ~exist(psdSingleFolder, 'dir')
       mkdir(psdSingleFolder);
end

save(psdSingleFolder+"/psdNN_lh.mat",'psdNN_lh');
save(psdSingleFolder+"/psdNN_lmh.mat",'psdNN_lmh');

%% Create SC-informed graph signal surrogates by randomization of real harmonics Fourier coefficients
% Inspiration from Petri and Van De Ville
% https://doi.org/10.1038/s41467-019-12765-7
% https://github.com/gpreti/GSP_StructuralDecouplingIndex

nSurr=19;

Xsurr=zeros(nROI,dataLength,nSurr);
XsurrHat=zeros(nROI,dataLength,nSurr);
XsurrCoupled=zeros(nROI,dataLength,nSurr);
XsurrDecoupled=zeros(nROI,dataLength,nSurr);
energySurrCoupled=zeros(nROI,nSurr,nScans);
energySurrDecoupled=zeros(nROI,nSurr,nScans);

% SPATIAL RANDOMIZATION
for s=1:nScans
    disp(s)
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(U,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);
        
        % Create surrogate data
        Xsurr(:,:,n)=U(:,:,s)*PHI*U(:,:,s)'*data(:,:,s).';
        % GFT of surrogate data with real harmonics
        % to compare against empirical SDI
        XsurrHat(:,:,n)=U(:,:,s)'*Xsurr(:,:,n);
        XsurrCoupled(:,:,n)=Ulow(:,:,s)*XsurrHat(:,:,n);
        XsurrDecoupled(:,:,n)=Uhigh(:,:,s)*XsurrHat(:,:,n);

        % norms  of the weights
        for r=1:nROI
            energySurrCoupled(r,n,s)=norm(XsurrCoupled(r,:,n));
            energySurrDecoupled(r,n,s)=norm(XsurrDecoupled(r,:,n));
        end
    end
end

save(psdFolder+"/energySurrCoupled.mat",'energySurrCoupled')
save(psdFolder+"/energySurrDecoupled.mat",'energySurrDecoupled')

%% Only for normLaplacian (VDV uses normLaplacian)
%% Create SC-ignorant graph signal surrogates by randomization of configuration model derived-SC harmonics Fourier coefficients

% Graph structure "randomization": configuration model (cm) used to generate degree-preserving graph 
Wnew=Anorm;

for s=1:nScans
    Wcm(:,:,s)=sum(Wnew(:,:,s),2)*sum(Wnew(:,:,s),2).'./sum(sum(Wnew(:,:,s),2));
    
    % Create surrogate signals based on configuration model graph : XrandSran
    Lcm(:,:,s)=diag(sum(Wnew(:,:,s),2))-Wcm(:,:,s);
    %%Laplacian Decomposition
    [tempUcm,tempLambdaLcm] = eig(Lcm(:,:,s));
    [Lambdacm(:,s), IndLcm]=sort(diag(tempLambdaLcm));
    Ucm(:,:,s)=tempUcm(:,IndLcm);
end

%%% TODO I SHOULD ACTUALLY PLOT THE wZC %%%

%plot_scans_and_mean(Lambdacm,"Laplacian eigvals/Variation of eigvecs" ...
%    ,"Spectral index","eigVal",nGroups,indsPDHC,subplotsInds,plotLegends, ...
%    saveFolder,boolSavePlots)

%%
nSurr=19;
XsurrCM=zeros(nROI,dataLength,nSurr);
XsurrHatCM=zeros(nROI,dataLength,nSurr);
XsurrCoupledCM=zeros(nROI,dataLength,nSurr);
XsurrDecoupledCM=zeros(nROI,dataLength,nSurr);
energySurrCoupledCM=zeros(nROI,nSurr,nScans);
energySurrDecoupledCM=zeros(nROI,nSurr,nScans);

% SPATIAL RANDOMIZATION
for s=1:nScans
    disp(s)
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(Ucm,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);
        
        % Create surrogate data
        XsurrCM(:,:,n)=Ucm(:,:,s)*PHI*Ucm(:,:,s)'*data(:,:,s).';
        % GFT of surrogate data with real harmonics
        % to compare against empirical SDI
        XsurrHatCM(:,:,n)=U(:,:,s)'*XsurrCM(:,:,n);
        XsurrCoupledCM(:,:,n)=Ulow(:,:,s)*XsurrHatCM(:,:,n);
        XsurrDecoupledCM(:,:,n)=Uhigh(:,:,s)*XsurrHatCM(:,:,n);

        % norms  of the weights
        for r=1:nROI
            energySurrCoupledCM(r,n,s)=norm(XsurrCoupledCM(r,:,n));
            energySurrDecoupledCM(r,n,s)=norm(XsurrDecoupledCM(r,:,n));
        end
    end
end

save(psdFolder+"/energySurrCoupledCM.mat",'energySurrCoupledCM')
save(psdFolder+"/energySurrDecoupledCM.mat",'energySurrDecoupledCM')


%% Plot timeseries

% Indicies of the data used
times=startInd:endInd;

% Only need the indicies for the first session for PD and HC
plot_timeseries(data,Xc,Xd,times,indsPDHC{1},nPD,"PD-patient-", ...
                saveFolder,boolSavePlots)
plot_timeseries(data,Xc,Xd,times,indsPDHC{3},nHC,"HC-patient-", ...
                saveFolder,boolSavePlots)

%% Visualization of the Connectivity Matrix's Mean and Variance %%
%%% PLOTTING %%%
stringMat='Adj'; % Which connectivity matrix to plot
% Only for Adj as log(-FC) does not work
boolLog=0; % Log-plots (The log scale allows you to better see contrast)

% Connectivity matrix to be plotted
switch stringMat
    case 'Adj'
        cMat=A; 
        cStr=" Adj ";
    case 'FC'
        cMat=FC;
        cStr=" FC ";
    otherwise
end

maxMatVar=max(var(cMat,0,3),[],"all");
plotType=["Mean","Var"];
% Plot Mean connectivity 
for i=0:1
    f=figure;
    for g=1:nGroups
        if(~i)
            toPlot=mean(cMat(:,:,indsPDHC{g}),3);
        else
            toPlot=var(cMat(:,:,indsPDHC{g}),0,3);
        end
        logStr="";
        if(boolLog)
            logStr="Log ";
            meanPlot=log(meanPlot);
            varPlot=log(varPlot);
        end
        subplot(3,2,subplotsInds{g});imagesc(toPlot);
        title(logStr+plotType(i+1)+cStr+plotLegends(g));
        xlabel('regions');ylabel('regions');colorbar; 
        if(~i) %meanFC limited to [0,1]
            clim([0,1]);
        else %Var limited from 0 to maxVar + 5%
            clim([0,maxMatVar+maxMatVar*0.05]);
        end
    end
    if(~i)
        subplot(3,2,5);
        imagesc(abs(mean(cMat(:,:,indsPDHC{1}),3)-mean(cMat(:,:,indsPDHC{2}),3)));
        colorbar;
        clim([0,1])
        subplot(3,2,6);
        imagesc(abs(mean(cMat(:,:,indsPDHC{3}),3)-mean(cMat(:,:,indsPDHC{4}),3)));
        colorbar;
        clim([0,1])
    end
    if(boolSavePlots)
        saveas(f,saveFolder+"/"+plotType(i+1)+"FC.png")
    end
end

%% Plot eigenvalues and structural harmonics and wZC

% Plot distribution of eigenvalues
eigFig=figure;
sgtitle("Eigenvalue distributions")
for g=1:nGroups
    subplot(2,2,subplotsInds{g})
    histfit(reshape(LambdaL(:,indsPDHC{g}),1,[]),44);
    title(plotLegends{g})
    xlabel("Eigenvalue");ylabel("Count")
    %ylim([0,])
end

saveas(eigFig,saveFolder+"/eigDistr.png")

% Plot eigenvalues (i.e. variation of eigenvectors)
plot_scans_and_mean(LambdaL,"Laplacian eigvals/Variation of eigvecs" ...
    ,"Spectral index","eigVal",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    saveFolder,boolSavePlots)

% Plot all of the eigenvectors values, mean over each session
meanEigVec=figure;
sgtitle('Mean Laplacian Eigvecs');
minU=min(U,[],"all");
maxU=max(U,[],"all");
for g=1:nGroups
    subplot(2,2,subplotsInds{g});imagesc(mean(U(:,:,indsPDHC{g}),3));
    title(plotLegends(g))
    xlabel('Spectral index');ylabel('regions')  
    colorbar;
    clim([minU,maxU]);
end

saveas(meanEigVec,saveFolder+"/meanEigVec.png")

% Plot num zero-corssings
plot_scans_and_mean(ZC,"Zero Crossings (ZC)" ...
    ,"Connectome harmonics","ZC",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    saveFolder,boolSavePlots)
% Plot weighted zero-crossings
plot_scans_and_mean(wZC,"Weigthed Zero Crossings (wZC)" ...
    ,"Connectome harmonics","wZC",nGroups,indsPDHC,subplotsInds,plotLegends, ...
    saveFolder,boolSavePlots)
%% Spectral X
%for each timepoint, you have one coefficient for each eigenvector, i.e., the contribution of each eigenvector to the functional signal

%figure;imagesc(X_hat(:,:,4));title('Spectral Coefficients - scan 1');colorbar;xlabel('time');ylabel('spectral index')

% These are simply the values of the signal in the freq domain. It makes
% sense that we get +- large values for low spectral indices because those
% are the indices belonging to aligned eigenvectors => we would expect the energy to be
% the largest among these if the signal is overall fairly aligned with the
% network. But for the PD RS data we don't see this...

%% PLOTTING
for g=1:nGroups
    % If we want the same colorbar we can add ,[-0.2,0.2] after 
    % Uhigh in imagesc
    figure;
    sgtitle(["Mean high and low eigvecs " plotLegends{g}])
    subplot(3,1,1);
    imagesc(mean(Uhigh(:,:,indsPDHC{g}),3));colorbar;title('High Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
    subplot(3,1,2);
    imagesc(mean(Umid(:,:,indsPDHC{g}),3));colorbar;title('Mid Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
    subplot(3,1,3);
    imagesc(mean(Ulow(:,:,indsPDHC{g}),3));colorbar;title('Low Frequency Structural Harmonics');
    xlabel('spectral index');ylabel('regions')
end
%% 
function plot_timeseries(meg,coupled,decoupled,times,ids,nPats,titStr,folder,boolSave)
    % meg: 3d Matrix to plot  time x regions x scans
    % coupled: coupled timeseries regions x time x scans
    % decoupled: decoupled timeseries regions x time x scans
    % times: plotStart:plotEnd range of ints Int
    % ids: ids for patients
    % nPats: number of PD or HC patients in total
    % titStr: string for the title ("PD id " or "HC id ")

    % Plots the MEG,coupled and decoupled timeseries for the given patients
    % left column ses1, right column ses 2.

    sess=["(ses 1)","(ses 2)"];
    SPinds={[1,3,5],[2,4,6]}; %indices subplot
    
    for i=1:nPats % Patients
        f = figure;
        f.Position = [100 100 1400 700];
        
        %Find min and max over sessions for the scan
        subj=ids(i);
        subjs=[subj,subj+nPats];
        megYmin=min(meg(times,:,subjs),[],"all");
        megYmax=max(meg(times,:,subjs),[],"all");
        cYmin=min(coupled(times,:,subjs),[],"all");
        cYmax=max(coupled(times,:,subjs),[],"all");
        dYmin=min(decoupled(times,:,subjs),[],"all");
        dYmax=max(decoupled(times,:,subjs),[],"all");
        
        % Column
        for j=1:2
            subplot(3,2,SPinds{j}(1));
            plot(meg(times,:,subj+(j-1)*nPats));
            title("MEG timecourses "+sess(j));xlabel('time');ylabel('amplitude')
            ylim([megYmin,megYmax]);

            subplot(3,2,SPinds{j}(2));
            plot(coupled(times,:,subj+(j-1)*nPats))
            title("Coupled signal timecourses "+sess(j));xlabel('time');ylabel('amplitude');
            ylim([megYmin,megYmax]);
            %ylim([cYmin,cYmax]);

            subplot(3,2,SPinds{j}(3));
            plot(decoupled(times,:,subj+(j-1)*nPats))
            title("Decoupled signal timecourses "+sess(j));xlabel('time');ylabel('amplitude');
            ylim([megYmin,megYmax]);
            %ylim([dYmin,dYmax]);
        end
        sgtitle(titStr+string(i))
        timeFolder=folder+"/time/";
        if ~exist(timeFolder, 'dir')
           mkdir(timeFolder)
        end
        if(boolSave)
            saveas(f,timeFolder+titStr+string(i)+"_timeseries.png")
        end
        %clf(f)
    end
end

function plot_scans_and_mean(toBePlotted,sgtit,xlab,ylab,nGroups,indsPDHC, ...
    subplotsInds,plotLegends,folder,savePlot)

indvFig=figure;
sgtitle("Individual "+sgtit);
% Plot individual scans and group mean
meanG=zeros(size(toBePlotted,1),nGroups);
stdG=zeros(size(toBePlotted,1),nGroups);

yMax=max(toBePlotted,[],"all");
for g=1:nGroups 
    subplot(2,2,subplotsInds{g});plot(toBePlotted(:,indsPDHC{g}));
    hold on;
    meanG(:,g)=mean(toBePlotted(:,indsPDHC{g}),2);
    stdG(:,g)=std(toBePlotted(:,indsPDHC{g}),0,2); 
    subplot(2,2,subplotsInds{g});plot(meanG(:,g),'k');
    title(plotLegends(g));xlabel(xlab);ylabel(ylab);
    ylim([0 yMax+0.05*yMax]);
end
if(savePlot)
    saveas(indvFig,folder+"/indv_"+ylab+".png")
end

% Compare mean between groups
meanFig=figure;
for g=1:nGroups
    errorbar(1:size(meanG, 1), meanG(:,g), stdG(:,g))
    xlabel(xlab);ylabel(ylab);
    hold on;
end
title("Mean "+sgtit);
legend(plotLegends,'Location','southeast');

if(savePlot)
    saveas(meanFig,folder+"/mean_"+ylab+".png")
end

end

