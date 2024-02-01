%% Chosen eig.vec. correlation mean and std for all groups 
chosenMeanCorrLowU=meanCorrLowU(i,dGroups)
chosenStdCorrLowU=stdCorrLowU(i,dGroups)
disp("----------")
chosenMeanCorrMediumU=meanCorrMediumU(k,dGroups)
chosenStdCorrMediumU=stdCorrMediumU(k,dGroups)
disp("----------")
chosenMeanCorrHighU=meanCorrHighU(j,dGroups)
chosenStdCorrHighU=stdCorrHighU(j,dGroups)
disp("----------")
if(~all(chosenMeanCorrLowU-chosenStdCorrLowU >= bCorr))
    disp("Mean-std Low Eigvec Corr < "+bCorr)
end
if(~all(chosenMeanCorrMediumU-chosenStdCorrMediumU >= bCorr))
    disp("Mean-std Medium Eigvec Corr < "+bCorr)
end
if(~all(chosenMeanCorrHighU-chosenStdCorrHighU >= bCorr))
    disp("Mean-std High Eigvec Corr < "+bCorr)
end
disp("----------")
%% Chosen energy correlation mean and std for all groups 
chosenMeanCorrLowX=meanCorrLowX(i,dGroups)
chosenStdCorrLowX=stdCorrLowX(i,dGroups)
disp("----------")
chosenMeanCorrMediumX=meanCorrMediumX(k,dGroups)
chosenStdCorrMediumX=stdCorrMediumX(k,dGroups)
disp("----------")
chosenMeanCorrHighX=meanCorrHighX(j,dGroups)
chosenStdCorrHighX=stdCorrHighX(j,dGroups)
disp("----------")
if(~all(chosenMeanCorrLowX-chosenStdCorrLowX >= bCorr))
    disp("Mean-std Low Energy Corr < "+bCorr)
end
disp("----------")
if(~all(chosenMeanCorrMediumX-chosenStdCorrMediumX >= bCorr))
    disp("Mean-std Medium Energy Corr < "+bCorr)
end
disp("----------")
if(~all(chosenMeanCorrHighX-chosenStdCorrHighX >= bCorr))
    disp("Mean-std High Energy Corr < "+bCorr)
end