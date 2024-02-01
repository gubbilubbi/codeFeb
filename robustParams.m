% Calculates how correlated mean eigvecs and energy are for diff. cutoffs
% Inspired by Fig. 4 in https://ieeexplore.ieee.org/document/7544580/

% For each group calculate the pearson correlation between the mean
% eigvec and mean energy for all frequencies and plot a boxplot for each
% freq with all groups.

% Used to evaluate the best cutoff value for Low/Medium/High across all
% groups. A good value is one with highly correlated eigenvectors 
% for different cutoff values.
%% Settings

% Manually selected cutoffs values
myLow=6;
myHigh=39;
% The lower and upper bound to test within from the chosen cutoff values
bound=1;
lowBound=bound;
highBound=bound;

%%
allLow=myLow-lowBound:myLow+highBound;
allHigh=myHigh-lowBound:myHigh+highBound;

savePlotsFolder="./gsp_data/dataFeb/"+folderName+"/cutoff_"+myLow+"_"+myHigh;

% The energy for the chosen cutoff
myEnergyX=load(savePlotsFolder+"/energyX.mat").energyX;

myLowMedium=myLow+1;
myHighMedium=myHigh-1;

lowRange=allLow(allLow~=myLow);
highRange=allHigh(allHigh~=myHigh);

for g=1:nGroups
    indL=0;
    indH=0;
    indM=0;
    % Correlation for lower mean freqs
    for lower=lowRange
        corrLowU(indL,g)=corr(mean(abs(U(:,1:myLow,indsPDHC{g})),[2,3]),mean(abs(U(:,1:lower,indsPDHC{g})),[2,3]));
        % Energy value to compare with
        newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+lower+"_"+myHigh+"/energyX.mat").energyX;
        corrLowX(indL,g)=corr(mean(myEnergyX(:,indsPDHC{g},1),2),mean(newEnergyX(:,indsPDHC{g},1),2));
        indL=indL+1;
    end
    % Correlation for higher mean freqs
    for upper=highRange
        corrHighU(indH,g)=corr(mean(abs(U(:,myHigh:end,indsPDHC{g})),[2,3]),mean(abs(U(:,upper:end,indsPDHC{g})),[2,3]));
        % Energy value to compare with
        newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+myLow+"_"+upper+"/energyX.mat").energyX;
        corrHighX(indH,g)=corr(mean(myEnergyX(:,indsPDHC{g},3),2),mean(newEnergyX(:,indsPDHC{g},3),2));
        indH=indH+1;
    end
    % Correlation for medium mean freqs
    for lowM=allLow
        for highM=allHigh
            % Energy value to compare to
            newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+lowM+"_"+highM+"/energyX.mat").energyX;
            if(lowM==myLowMedium && highM==myHighMedium)
                % Skip if the medium range is the same as the chosen one
            else
                corrMediumU(indM,g)=corr(mean(abs(U(:,myLowMedium:myHighMedium,indsPDHC{g})),[2,3]),mean(abs(U(:,lowM:highM,indsPDHC{g})),[2,3]));
                corrMediumX(indM,g)=corr(mean(myEnergyX(:,indsPDHC{g},2),2),mean(newEnergyX(:,indsPDHC{g},2),2));
                indM=indM+1;
            end
        end
    end
end
disp("N lower per boxplot: " +indL)
disp("N medium per boxplot: " +indM)
disp("N high per boxplot: " +indH)

% Boxplots for the correlation between eigvecs for different cutoffs
f=figure;
boxplot(corrLowU)
title("Correlation group-wise for low eigenvectors");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrLowU.png")

f=figure;
boxplot(corrMediumU);
title("Correlation group-wise for medium eigenvectors");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrMediumU.png")

f=figure;
boxplot(corrHighU);
title("Correlation group-wise for high eigenvectors");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrHighU.png")

% Boxplots for the correlation between energy for different cutoffs
f=figure;
boxplot(corrLowX)
title("Correlation group-wise for low signals");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrLowX.png")

f=figure;
boxplot(corrMediumX);
title("Correlation group-wise for medium signals");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrMediumX.png")

f=figure;
boxplot(corrHighX);
title("Correlation group-wise for high signals");
ylim([0,1]);
saveas(f,savePlotsFolder+"/corrHighX.png")