% Plots one figure for Low/Medium/High correlation of eig.vec.s and energy
% One line plot for each of the 4 different groups with mean and std


% Mean eigenvec correlation
figure;
for g=1:nGroups
    errorbar(1:length(lows),meanCorrLowU(:,g),stdCorrLowU(:,g))
    hold on;
end
xlim([0,length(lows)+1]);
title("Mean Eig.vec. Low corr");
legend(plotLegends,'Location','southeast');

figure;
for g=1:nGroups
    errorbar(1:length(highs),meanCorrHighU(:,g),stdCorrHighU(:,g))
    hold on;
end
xlim([0,length(highs)+1]);
title("Mean Eig.vec. High corr");
legend(plotLegends,'Location','southwest');

figure;
for g=1:nGroups
    errorbar(1:size(corrMediumU,1),meanCorrMediumU(:,g),stdCorrMediumU(:,g))
    hold on;
end
xlim([0,size(corrMediumU,1)+1]);
title("Mean Eig.vec. Medium corr");
legend(plotLegends,'Location','southeast');

% Mean Energy correlation
figure;
for g=1:nGroups
    errorbar(1:length(lows),meanCorrLowX(:,g),stdCorrLowX(:,g))
    hold on;
end
xlim([0,length(lows)+1]);
title("Mean Energy Low Corr");
legend(plotLegends,'Location','southeast');

figure;
for g=1:nGroups
    errorbar(1:length(highs),meanCorrHighX(:,g),stdCorrHighX(:,g))
    hold on;
end
xlim([0,length(highs)+1]);
title("Mean Energy High Corr");
legend(plotLegends,'Location','southwest');

figure;
for g=1:nGroups
    errorbar(1:size(corrMediumX,1),meanCorrMediumX(:,g),stdCorrMediumX(:,g))
    hold on;
end
xlim([0,size(corrMediumX,1)+1]);
title("Mean Energy Medium Corr");
legend(plotLegends,'Location','southeast');