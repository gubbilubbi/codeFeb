% Calculates how correlated mean eigvecs and energy are for diff. cutoffs
% Inspired by Fig. 4 in https://ieeexplore.ieee.org/document/7544580/

% For each group calculate the pearson correlation between the mean
% eigvec and mean energy for all frequencies 

% Used to evaluate the best cutoff value for Low/Medium/High across all
% groups. A good value is one with highly correlated eigenvectors 
% for different cutoff values.

% Inspiration from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5400112/
% for correlation of eigenvectors and
% https://ieeexplore.ieee.org/document/9044786 for correlation of energy

% the lower cutoffs tested are the set lowestLow+bound to highestLow-bound
    % and for each there are "bound*2 many" correlation values saved.
    % between the current low and the +-bound low
    % Example with bound 2 =>
    % -----4-5-[6]-7-8------=> 4/6,5/6,6/7,6/8. i.e. all pairs from the
    % currently chosen low cutoff that are within bound steps from it.

% the medium cutoffs tested are lowestLow+bound+1 to highestLow-bound-1 =>
    % there are nLow*nHigh different medium cutoffs because for every low
    % there is a high which gives us a unique medium cutoff. Furthermore,
    % for each medium cutoff there are (bound*2+1)^2-1 corr values saved 
    % as each medium has a low and high boundary that can be moved +-bound
    % many steps. For each low bound there are then bound*2+1 high bounds
    % meaning in total there are (bound*2+1)^2 combinations of medium
    % cutoffs within bound steps, but minus 1 because 1 will be the same as
    % the one we currently have.
    % Example with bound 1 =>
    % -----4-5-6-[7]-8-9----10-11-[12]-13-14-=> 
    % 6/11,6/12,6/13,7/11,7/13,8/11,8/12,8/13
    % i.e. all pairs from the
    % currently chosen medium cutoff that are within bound steps from it
    % in either one or both directions.

% the higher cutoffs tested are the set lowestLow+bound to highestLow-bound
    % and for each there are "bound*2 many" correlation values saved.
    % between the current high and the +-bound high
%% Settings
% The lower and upper bound to test within from the chosen cutoff values.
% Can be chosen to be non-symmetric around the cutoff frequency.
lowBound=bound;
highBound=bound;

lowestLow=3;
highestLow=18;
lowestHigh=27;
highestHigh=41;
%%

% Selected cutoffs
lows=lowestLow+bound:highestLow-bound;
highs=lowestHigh+bound:highestHigh-bound;

nLows=length(lows);
nHighs=length(highs);
nMedium=nLows*nHighs;

% Per medium cutoff we have (bound*2+1)^2 comparisons to make
% as we can go -bound to +bound for low and high => for every
% low +-bound freq cutoff for medium we have (bound*2+1) high ones =>
% (bound*2+1)^2 but -1 as one will be the same
nMediumBounds=((bound*2+1)*(bound*2+1))-1; 

corrLowU=zeros(nLows,bound*2,nGroups);
corrMediumU=zeros(nMedium,nMediumBounds,nGroups);
corrHighU=zeros(nHighs,bound*2,nGroups);

corrLowX=zeros(nLows,bound*2,nGroups);
corrMediumX=zeros(nMedium,nMediumBounds,nGroups);
corrHighX=zeros(nHighs,bound*2,nGroups);

for i=1:length(lows)
    myLow=lows(i);
    for j=1:length(highs)
        myHigh=highs(j);
        
        % The current index for medium matrices
        % Goes from 1 to nLows*nHighs
        k=nHighs*(i-1)+j;

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
            indL=1;
            indH=1;
            indM=1;
            % Correlation for lower mean freqs
            for lower=lowRange
                corrLowU(i,indL,g)=corr(mean(abs(U(:,1:myLow,indsPDHC{g})),[2,3]),mean(abs(U(:,1:lower,indsPDHC{g})),[2,3]));
                % Energy value to compare against
                newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+lower+"_"+myHigh+"/energyX.mat").energyX;
                corrLowX(i,indL,g)=corr(mean(myEnergyX(:,indsPDHC{g},1),2),mean(newEnergyX(:,indsPDHC{g},1),2));
                indL=indL+1;
            end
            % Correlation for higher mean freqs
            for upper=highRange
                corrHighU(j,indH,g)=corr(mean(abs(U(:,myHigh:end,indsPDHC{g})),[2,3]),mean(abs(U(:,upper:end,indsPDHC{g})),[2,3]));
                % Energy value to compare against
                newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+myLow+"_"+upper+"/energyX.mat").energyX;
                corrHighX(j,indH,g)=corr(mean(myEnergyX(:,indsPDHC{g},3),2),mean(newEnergyX(:,indsPDHC{g},3),2));
                indH=indH+1;
            end
            % Correlation for medium mean freqs
            for lowM=allLow+1
                for highM=allHigh-1
                    % Skip if the medium range is the same as the chosen one
                    if(~(lowM==myLowMedium && highM==myHighMedium))
                        % Energy value to compare against
                        newEnergyX=load("./gsp_data/dataFeb/"+folderName+"/cutoff_"+string(lowM-1)+"_"+string(highM+1)+"/energyX.mat").energyX;
                        corrMediumU(k,indM,g)=corr(mean(abs(U(:,myLowMedium:myHighMedium,indsPDHC{g})),[2,3]),mean(abs(U(:,lowM:highM,indsPDHC{g})),[2,3]));
                        corrMediumX(k,indM,g)=corr(mean(myEnergyX(:,indsPDHC{g},2),2),mean(newEnergyX(:,indsPDHC{g},2),2));
                        indM=indM+1;
                    end
                end
            end
        end
    end
end
% Mean and std correlation values for eigenvectors
meanCorrLowU=squeeze(mean(corrLowU,2));
stdCorrLowU=squeeze(std(corrLowU,0,2));
meanCorrMediumU=squeeze(mean(corrMediumU,2));
stdCorrMediumU=squeeze(std(corrMediumU,0,2));
meanCorrHighU=squeeze(mean(corrHighU,2));
stdCorrHighU=squeeze(std(corrHighU,0,2));

% Mean and std correlation values for the energy
meanCorrLowX=squeeze(mean(corrLowX,2));
stdCorrLowX=squeeze(std(corrLowX,0,2));
meanCorrMediumX=squeeze(mean(corrMediumX,2));
stdCorrMediumX=squeeze(std(corrMediumX,0,2));
meanCorrHighX=squeeze(mean(corrHighX,2));
stdCorrHighX=squeeze(std(corrHighX,0,2));

