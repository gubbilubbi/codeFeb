%% Generate Brain plots
% Atlas + Selection
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

% Last dimension for low,medium,high
meanU=zeros(nROI,nGroups,nFreqGroups);
stdU=zeros(nROI,nGroups,nFreqGroups);

% Split eigenvectors into low, medium and high frequencies 
for i=1:nFreqGroups
    for g=1:nGroups
        meanU(:,g,i)=mean(abs(Ufreqs(:,fgNN{i,g},indsPDHC{g},i)),[2,3]);
        stdU(:,g,i)=std(abs(Ufreqs(:,fgNN{i,g},indsPDHC{g},i)),0,[2,3]);
    end
end

brainFolder=saveFolder+"/brainPlots/";
brainSubFolders={brainFolder+"eigVecs/",... 
                brainFolder+"lEigVec/",...
                brainFolder+"mEigVec/",...
                brainFolder+"hEigVec/",...
                brainFolder+"lSignal/",...
                brainFolder+"mSignal/",...
                brainFolder+"hSignal/"};

if ~exist(brainFolder, 'dir')
   mkdir(brainFolder)
   mkdir(brainSubFolders{1})
   mkdir(brainSubFolders{2})
   mkdir(brainSubFolders{3})
   mkdir(brainSubFolders{4})
   mkdir(brainSubFolders{5})
   mkdir(brainSubFolders{6})
   mkdir(brainSubFolders{7})
end

%%
%%% Settings %%%
thrVal=0.8;
bmeshBool=0.2;
expSphereSize=2;


disp("Plotting Single Group Mean Eigvecs")

vecs=[1,2,3,4,nROI-3,nROI-2,nROI-1,nROI];

% Huang
%CA=[0,max(mean(abs(U(:,vecs,:)),3),[],"all")];
% VDV 
CA=[min(mean(U(:,vecs,:),3),[],"all"),max(mean(U(:,vecs,:),3),[],"all")];

for v=vecs
    disp("Vector: "+v)
    f=figure;
    f.Position(3:4)=[700 600];
    sgtitle("Mean eigVec"+v,'horizontalAlignment','right');
    
    % Huang
    % meanPD1=mean(squeeze(abs(U(:,v,indsPDHC{1}))),2);
    % meanPD2=mean(squeeze(abs(U(:,v,indsPDHC{2}))),2);
    % meanHC1=mean(squeeze(abs(U(:,v,indsPDHC{3}))),2);
    % meanHC2=mean(squeeze(abs(U(:,v,indsPDHC{4}))),2);
    
    % VDV
    meanPD1=mean(squeeze(U(:,v,indsPDHC{1})),2);
    meanPD2=mean(squeeze(U(:,v,indsPDHC{2})),2);
    meanHC1=mean(squeeze(U(:,v,indsPDHC{3})),2);
    meanHC2=mean(squeeze(U(:,v,indsPDHC{4})),2);

    plot_on_brain([meanPD1,meanPD2,meanHC1,meanHC2],0,CodeBook,CA, ...
        ["PD1","PD2","HC1","HC2"],0,thrVal,expSphereSize,bmeshBool);
    % always 0 threshold because we want to 
    % compare eigenvectors within subjects and between. Setting a thr
    % will make this difficut because for low eig vecs everything will be
    % close to 0.

    if(boolSavePlots)
        saveas(f,brainSubFolders{1}+"eigVec"+v+".png")
    end
    clf(f); % Clear figure to avoid a memory leak
end
disp("Single Group Mean Eigvecs Done")


% strPlots=["Low","Medium","High"];
% % for plotting with and without threshold
% 
% for boolThr=[0,1]
%     disp("BoolThr: "+boolThr)
%     disp("Plotting mean eigvecs")
%     for i=1:length(strPlots)
%         plotData=squeeze(meanU(:,:,i));%squeeze(vecData{i});
%         strPlot=strPlots(i);
% 
%         f=figure;
%         f.Position(3:4)=[700 600];
%         sgtitle("Mean mag. "+strPlot+" freq. eigvecs",'horizontalAlignment','right');
%         plot_on_brain(plotData,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolThr,thrVal,expSphereSize,bmeshBool)
%         if(boolSavePlots)
%             saveas(f,brainSubFolders{i+1}+"mean"+strPlot+"EigVecs_thr"+boolThr+".png")
%         end
%         clf(f);
% 
%         % Only plot diff in eigvecs if we have FC and not SC
%         % because SC will be same for diff PD and HC sessions
%         if(boolFC) 
% 
%             f=figure;
%             sgtitle("Abs Diff PD "+strPlot+" eigvecs",'horizontalAlignment','right')
%             plot_on_brain([plotData(:,1),plotData(:,2),abs(plotData(:,1)-plotData(:,2))], ...
%                 0,CodeBook,0,["PD1","PD2","abs diff"],boolThr,thrVal,expSphereSize,bmeshBool)
%             if(boolSavePlots)
%                 saveas(f,brainSubFolders{i+1}+"mean"+strPlot+"EigVecDiffPD_thr"+boolThr+".png")
%             end
%             clf(f);
% 
%             f=figure;
%             sgtitle("Abs Diff HC "+strPlot+" eigvecs",'horizontalAlignment','right')
%             plot_on_brain([plotData(:,3),plotData(:,4),abs(plotData(:,3)-plotData(:,4))], ...
%                 0,CodeBook,0,["HC1","HC2","abs diff"],boolThr,thrVal,expSphereSize,bmeshBool)
%             if(boolSavePlots)
%                 saveas(f,brainSubFolders{i+1}+"mean"+strPlot+"EigVecDiffHC_thr"+boolThr+".png")
%             end
%             clf(f);
%         end
%         disp("Done " + i)
%     end
% 
%     disp("Plotting mean signals")
%     for i=1:length(strPlots)
%         plotData=squeeze(meanX(:,:,i));
%         strPlot=strPlots(i);
% 
%         f=figure;
%         f.Position(3:4)=[700 600];
%         sgtitle("Mean abs. "+strPlot+" freq. signal",'horizontalAlignment','right');
%         plot_on_brain(plotData,0,CodeBook,0,["PD1","PD2","HC1","HC2"],boolThr,thrVal,expSphereSize,bmeshBool)
%         if(boolSavePlots)
%             saveas(f,brainSubFolders{i+4}+"mean"+strPlot+"X_thr"+boolThr+".png")
%         end
%         clf(f);
% 
%         f=figure;
%         sgtitle("Abs Diff PD mean "+strPlot+" signal",'horizontalAlignment','right');
%         plot_on_brain([plotData(:,1),plotData(:,2),abs(plotData(:,1)-plotData(:,2))], ...
%             0,CodeBook,0,["PD1","PD2","abs diff"],boolThr,thrVal,expSphereSize,bmeshBool)
%         if(boolSavePlots)
%             saveas(f,brainSubFolders{i+4}+"mean"+strPlot+"XDiffPD_thr"+boolThr+".png");
%         end
%         clf(f);
% 
%         f=figure;
%         sgtitle("Abs Diff HC mean "+strPlot+" signal",'horizontalAlignment','right');
%         plot_on_brain([plotData(:,3),plotData(:,4),abs(plotData(:,3)-plotData(:,4))], ...
%             0,CodeBook,0,["HC1","HC2","abs diff"],boolThr,thrVal,expSphereSize,bmeshBool)
%         if(boolSavePlots)
%             saveas(f,brainSubFolders{i+4}+"mean"+strPlot+"XDiffHC_thr"+boolThr+".png");
%         end
%         clf(f);
%         disp("Done " + i)
%     end
% end
% disp("Done with brain plots")

%%
% Test to see if the colors work as intended
% boolThr=0;
% f=figure;
% valter=abs(meanU(:,1,1)-meanU(:,2,1));
% valter(20)=0.03;
% plot_on_brain([meanU(:,1,1),meanU(:,2,1),abs(meanU(:,1,1)-meanU(:,2,1))], ...
%                 0,CodeBook,0,["PD1","PD2","abs diff"],boolThr)