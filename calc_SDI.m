% Author: Valter Lundegårdh
% Inspiration from Petri and Van De Ville
% https://doi.org/10.1038/s41467-019-12765-7

% Calculates the structural decoupling index for the empirical data
% the strucutural connectivity informed/ignorant surrogates and plots
% the results in brain graphs if showPlots is 1.

meanCoupling=zeros(nROI,nGroups);
meanDecoupling=zeros(nROI,nGroups);

% energyX (nROI x nScans*2 x nFreqGroups) contains the l2-norm
% for all scans' regions for every frequency
for i=1:nGroups
    meanCoupling(:,i)=mean(energyX(:,indsPDHC{i},1),2); %mean_c=mean(N_c,2); %average coupling
    meanDecoupling(:,i)=mean(energyX(:,indsPDHC{i},2),2); %mean_d=mean(N_d,2); %average decoupling
end

meanSDI=meanDecoupling./meanCoupling; % mean SDI per group 
indvSDI=energyX(:,:,2)./energyX(:,:,1); % SDI per scan



% Load energy for surrogate data (SC-informed)
% nROI x nSurrs x nScans for surrogate data
energySurrCoupled=load(saveFolder+"/energySurrCoupled.mat").energySurrCoupled; 
energySurrDecoupled=load(saveFolder+"/energySurrDecoupled.mat").energySurrDecoupled;
% SDI for every scan and surrogates
SDI_surr=energySurrDecoupled./energySurrCoupled; 
% mean across surrogates
meanSDIsurr=squeeze(mean(SDI_surr,2)); 
% Mean across surrogates and scans
for g=1:nGroups
    % mean SDI for all surrogates and scans
    meanSDIsurrScans(:,g)=mean(meanSDIsurr(:,indsPDHC{g}),2); 
end


% nROI x nSurrs x nScans for configuration model surrogate data
% only we have normLaplacian because Petri and VDV's CM is only defined
% for normLaplacian.
if(isequal(methodSC,'normLaplacian'))
    energySurrCoupledCM=load(saveFolder+"/energySurrCoupledCM.mat").energySurrCoupledCM; 
    energySurrDecoupledCM=load(saveFolder+"/energySurrDecoupledCM.mat").energySurrDecoupledCM;
    % SDI for every scan and surrogates
    SDI_surr_cm=energySurrDecoupledCM./energySurrCoupledCM; 
    % mean across surrogates
    meanSDIsurrCM=squeeze(mean(SDI_surr_cm,2));
    % Mean across surrogates and scans
    for g=1:nGroups
        % mean SDI for all surrogates and scans
        meanSDIsurrScansCM(:,g)=mean(meanSDIsurrCM(:,indsPDHC{g}),2);
    end
end

% Significant SDI per group vs SC-informed surrogate

% for every scan, find threshold for min and max across surrogates
for s=nScans
    min_SDI_surr(:,s)=min(SDI_surr(:,:,s)')';
    max_SDI_surr(:,s)=max(SDI_surr(:,:,s)')';
    % select significant SDI for each scan, across surrogates 
    % individual thr, first screening
    %for each scan, I threshold the ratio based on individual ratio's surrogate distribution 
    SDI_thr_min(:,s)=indvSDI(:,s)<min_SDI_surr(:,s);
    SDI_thr_max(:,s)=indvSDI(:,s)>max_SDI_surr(:,s);
end

%%threshold empirical mean ratios
meanSDIthr=ones(nROI,nGroups); %1 if not significant

for g=1:nGroups
    % Amounts of detections per region (if region SDI is larger/smaller than
    % the max/min for all surrogates for the group)
    detect_min=sum(SDI_thr_min(:,indsPDHC{g})');
    detect_max=sum(SDI_thr_max(:,indsPDHC{g})'); 
    
    %%for every region, test across subjects 0.05, correcting for the number of
    %%tests (regions), 0.05/nROI (Bonferroni correction)
    x=0:1:100; %100 points
    
    % ???
    y=binocdf(x,100,0.05,'upper'); 
    THRsubjects=x(min(find(y<0.05/nROI))); 
    THRsubjects=floor(length(indsPDHC{g})/100*THRsubjects)+1;

    % ???
    SDI_sig_higher=detect_max>THRsubjects;
    SDI_sig_lower=detect_min>THRsubjects;
    
    % ???
    SDI_sig_tot_positions=[find(SDI_sig_higher==1),find(SDI_sig_lower==1)];
    SDI_sig_tot_positions=sort(unique(SDI_sig_tot_positions));
    
    meanSDIthr(SDI_sig_tot_positions,g)=meanSDI(SDI_sig_tot_positions,g);
end

% Plot SDI figures
% N.B. the "figure" row in PlotBrainGraph called by PlotGraph 
% has been commeneted out to enable the figure to be cleared after
% plotting here. 
if(showSDIPlots)
% Fig. 2A
    % saturate=1;
    % f=figure;
    % CC2=log2(meanSDIsurrScansCM(:,g)); 
    % PlotGraph;title('meanSDIcm '+ plotLegends(g))
    % clf(f);

%for g =1:nGroups
% Fig. 2B
    % saturate=1;
    % f=figure;
    % CC2=log2(meanSDIsurrScans(:,g)); 
    % PlotGraph;title('meanSDIsurr '+ plotLegends(g))
    % clf(f);

% Fig. 2C
    % saturate=1;
    % f=figure;
    % CC2=log2(meanSDIthr(:,g)); 
    % PlotGraph;title('meanSDIthr '+ plotLegends(g))
    % clf(f);
%end

% Fig. S3
    % saturate=0;
    % CC2=log2(meanSDIthr(:,g)./meanSDIsurrScans(:,g)); 
    % PlotGraph;title('meanSDIdiff'+ plotLegends(g))
    % clf(f);
end