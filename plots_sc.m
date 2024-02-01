%% TODO
% Maybe need the change the ylim for plotting SC if it is not bound
% to [0,1]? 


%% Settings
stringMat='Adj'; % Which connectivity matrix to plot
% boolLog is only for Adj as log with negative FC values => error
boolLog=0; % Log-plots (The log scale allows you to better see contrast)

%% Weigthed zero-crossings
% Alternative variation metric
% only considers a change in sign between regions

ZC=zeros(nROI,nSubjs);
wZC=zeros(nROI,nSubjs);

for s=1:nSubjs
    for u=1:nROI %for each eigenvector
        UU=U(:,u,s);
        sumZC=0; %ZC summ
        sumWZC=0; %wZC summ
        %for each eigvec, check all sign changes between all regions.
        for i=1:nROI-1 
            for j=i+1:nROI
                if (UU(i)*UU(j))<0 % if signals are of opposite signs
                    sumWZC=sumWZC+A(i,j,s);
                    sumZC=sumZC+1;
                end
                ZC(u,s)=sumZC; % Num of zero-crossings
                wZC(u,s)=sumWZC; % Total weighted zero-crossings
            end
        end
    end
end

%% Visualization of the Connectivity Matrix's Mean and Variance %%

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

% Find the max variance for each group and use the max over the groups
% as ylim for the variance of the matrices.
maxGroupVar=zeros(4,1);
for g=1:nGroups
    maxGroupVar(g,1)=max(var(cMat(:,:,indsPDHC{g}),0,3),[],"all");
end
maxMatVar=max(maxGroupVar,[],"all");
plotType=["Mean","Var"];

% Plot Mean connectivity 
for i=0:1 % For mean and variance
    f=figure;
    f.Position=[100 100 600 700];
    
    for g=1:nGroups
        if(~i)
            toPlot=mean(cMat(:,:,indsPDHC{g}),3);
        else
            toPlot=var(cMat(:,:,indsPDHC{g}),0,3);
        end
        logStr="";
        if(boolLog)
            logStr="Log";
            toPlot=log(toPlot);
        end
        subplot(3,2,subplotsInds{g});imagesc(toPlot);
        title(logStr+" "+plotType(i+1)+cStr+plotLegends(g));
        xlabel('regions');ylabel('regions');colorbar; 
        if(~i) %meanFC limited to [-1,1]
            switch stringMat
                % TODO maybe change the limit for SC?
                case 'Adj' % If Adj then ylim is [0,1] 
                    clim([0,1]);
                case 'FC' % If FC then ylim is [-1,1] (corr)
                    clim([-1,1]); 
                otherwise
            end
            
        else %Var limited from 0 to maxVar + 5%
            clim([0,maxMatVar+0.001]);
        end
    end
    % Plot differences
    if(~i)
        subplot(3,2,5);
        imagesc(abs(mean(cMat(:,:,indsPDHC{1}),3)-mean(cMat(:,:,indsPDHC{2}),3)));
        colorbar;
        title("Abs Diff in mean "+stringMat+" PD ses1 and ses2");
        xlabel('regions');ylabel('regions');colorbar; 
        clim([0,2])
        subplot(3,2,6);
        imagesc(abs(mean(cMat(:,:,indsPDHC{3}),3)-mean(cMat(:,:,indsPDHC{4}),3)));
        colorbar;
        title("Abs Diff in mean "+stringMat+" HC ses1 and ses2");
        xlabel('regions');ylabel('regions');colorbar; 
        clim([0,2])
    end
    if(boolSavePlots)
        saveas(f,dataFolder+"/"+plotType(i+1)+logStr+stringMat+".png")
    end
end

%% Plot disitrubtion of eigenvalues

nBins=44;

% TODO: Temporary solution to find the max value of all histograms
% to more easily compare the plots
yt=zeros(nGroups,1);
for g=1:nGroups
    histfit(reshape(LambdaL(:,indsPDHC{g}),1,[]),nBins);
    yt(g)=max(get(gca, 'YTick'),[],"all");
end
set(gcf,'Visible','off') %don't show these histograms
yMax=max(yt,[],"all");

% Plot distribution of eigenvalues
eigFig=figure;
sgtitle("Eigenvalue distributions")

% Plot distribution of eigenvalues
for g=1:nGroups
    subplot(2,2,subplotsInds{g})
    histfit(reshape(LambdaL(:,indsPDHC{g}),1,[]),nBins);
    title(plotLegends{g})
    xlabel("Eigenvalue");ylabel("Count")
    ylim([0,yMax*1.05]);
end

saveas(eigFig,dataFolder+"/eigDistr.png")

%% Plot eigenvalues and structural harmonics and (weighted) zero crossings

% Plot eigenvalues (i.e. variation of eigenvectors)
plot_subjects_and_mean(LambdaL,"Laplacian eigvals/Variation of eigvecs" ...
    ,"Spectral index","eigVal",indsPDHC,subplotsInds,plotLegends, ...
    dataFolder,boolSavePlots)

% Plot the eigenvectors' values, mean over each session
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
saveas(meanEigVec,dataFolder+"/meanEigVec.png")

% Plot zero-crossings
plot_subjects_and_mean(ZC,"Zero Crossings (ZC)" ...
    ,"Connectome harmonics","ZC",indsPDHC,subplotsInds,plotLegends, ...
    dataFolder,boolSavePlots)

% Plot weighted zero-crossings
plot_subjects_and_mean(wZC,"Weigthed Zero Crossings (wZC)" ...
    ,"Connectome harmonics","wZC",indsPDHC,subplotsInds,plotLegends, ...
    dataFolder,boolSavePlots)


%% Functions

function plot_subjects_and_mean(toBePlotted,sgtit,xlab,ylab,indsPDHC, ...
    subplotsInds,plotLegends,folder,savePlot)

% Parameters
    % toBePlotted: a matrix (nROIxnSubjs) e.g. eigVals or wZC
    % sgtit: string used for the title of the figure
    % xlab: x-axis label
    % ylab: y-axis label
    % indsPDHC: cell with one entry for each group of which subjects 
    %           that belong to that group.
    % subplotsInds: cell with one entry for each group of which subplot the
    %               group should be plotted in
    % plotLegends: string array, one for each group used in the legend
    % folder: path to save folder
    % savePlot: If to save the figures

    nGroups=length(indsPDHC);

    indvFig=figure;
    sgtitle("Individual "+sgtit);
    % Plot individual subjects and group mean
    meanG=zeros(size(toBePlotted,1),nGroups);
    stdG=zeros(size(toBePlotted,1),nGroups);
    
    yMax=max(toBePlotted,[],"all");
    for g=1:nGroups 
        % Plot all individial values per group
        subplot(2,2,subplotsInds{g});plot(toBePlotted(:,indsPDHC{g}));
        hold on;
        meanG(:,g)=mean(toBePlotted(:,indsPDHC{g}),2);
        stdG(:,g)=std(toBePlotted(:,indsPDHC{g}),0,2);

        % Plot the mean group value
        subplot(2,2,subplotsInds{g});plot(meanG(:,g),'k');
        title(plotLegends(g));xlabel(xlab);ylabel(ylab);
        ylim([0 yMax*1.05]);
    end
    if(savePlot)
        saveas(indvFig,folder+"/indv_"+ylab+".png")
    end
    
    % Compare the mean and std between groups
    meanFig=figure;
    for g=1:nGroups
        errorbar(1:size(meanG,1),meanG(:,g),stdG(:,g))
        xlabel(xlab);ylabel(ylab);
        hold on;
    end
    title("Mean "+sgtit);
    legend(plotLegends,'Location','southeast');
    
    if(savePlot)
        saveas(meanFig,folder+"/mean_"+ylab+".png")
    end
end