%% Statistical tests
% Plots the significant box plots with p-value and effect size
% Also plots all the p-values from the statistical test
% Also possible to save the values in a table

function [testRejectNull,testPvals]=stat_tests(vals,nROI,nFreqGroups,indsPDHC, ...
                    testType,folder,savePlot,saveStr)

%testType can be t, mwu or ks
disp(testType+"-Test");

% Colors for each region
hemisphereColors=[repmat("b",5,1);repmat("r",4,1);repmat("y",3,1);...
               repmat("g",2,1);repmat("c",4,1);repmat("m",4,1)];
cortexColors=[hemisphereColors;hemisphereColors];
% Used for single regions colors and all regions colors (over all freqs)
% when plotting the p-value distributions

switch nFreqGroups

    case 2
        cortexColors=[cortexColors;cortexColors]; 
    case 3
        cortexColors=[cortexColors;cortexColors;cortexColors]; 
    otherwise

end



allBoxes=0; % allBoxes or "" for all boxplots or only significant
if(allBoxes)
    boxStr="all";
else
    boxStr="";
end


whenToPlot="removeHCSig"; %allBoxes, ifSig, removeHCSig
disp("Plotting settings: "+whenToPlot);

% Cell with each test we perform (e.g. [1,3] is PD1 vs HC1)
testInds={[3,4],[1,2],[1,3],[1,4],[2,3],[2,4]};
nTests=length(testInds);

% Cells for results from t-test
%tTestRejectNull=cell(nROI,nFreqGroups,nTests);
%tTestPvals=cell(nROI,nFreqGroups,nTests);

% Cells for results from staistical test
testRejectNull=cell(nROI,nFreqGroups,nTests);
testPvals=cell(nROI,nFreqGroups,nTests);
effectSizes={nROI,nFreqGroups,nTests};

% Strings for plotting
tableStrs=["PD1","PD2","HC1","HC2"];
switch nFreqGroups        
    case 2
        freqStrings=["Low","High"];
    case 3
        freqStrings=["Low","Medium","High"];
    otherwise

end

% Get the names of the regions for the rows in the tables
regionTable=readtable('/Users/Mortimtosh/Desktop/kth-MSc/code/kth/gravity_centers_HCP_atlas_sorted.csv');

% The regions' names
switch nFreqGroups
    case 1 % If only 1 frequency (e.g. when testing eigvals) then we have
           % spectral indices (1:44) and not regions.
        regions=1:nROI; 
        % Names for x-axis 
        regionsFreqs=(1:nROI).'; 
    case 2
        regions=table2array(regionTable(:,1));
        % Names for x-axis 
        regionsFreqs=["Low_"+regions;"High_"+regions]; 
    case 3 % When low, medium and high freqs each regions has a name
        regions=table2array(regionTable(:,1));
        % Names for x-axis 
        regionsFreqs=["Low_"+regions;"Medium_"+regions;"High_"+regions]; 
end

testFolder=folder+"/tests/";
if ~exist(testFolder, 'dir')
mkdir(testFolder)
end

for t=1:nTests
    t1=testInds{t}(1);
    t2=testInds{t}(2);

    figTest=figure('Visible','off');
    hold on;
    figTest.Position=[100 100 1300 700];
    title(testType+"-test Null rejected "+saveStr+" "+tableStrs(t1)+"vs"+tableStrs(t2),'FontSize', 12);

    boxPos=1; %Current plotting position for each pair of boxplots
    boxNames={}; % x-labels for the boxplot pairs
    boxPvals={}; % Stars/p-values for each boxplot pair
    boxEffects={}; % stars/effect size for each boxplot pair
    significantColors={}; % Color of the region for a boxplot pair
    boxDist=4; % Distance between two pairs of boxplots

    nSignificant=0;
    for f=1:nFreqGroups
        for i=1:nROI
            
            % Values for each group
            valsG1=vals(i,indsPDHC{t1},f);
            valsG2=vals(i,indsPDHC{t2},f);

            % Group mean and std and number of subject
            mean1=mean(valsG1,2);
            mean2=mean(valsG2,2);
            std1=std(valsG1,0,2);
            std2=std(valsG2,0,2);
            n1=length(indsPDHC{t1}); %n PD
            n2=length(indsPDHC{t2}); %n HC
            
            % Calculate the effect
            stdPooled=sqrt(((n1-1)*std1^2+(n2-1)*std2^2)/(n1+n2-2));
            effectSizes{i,f,t}=abs(mean1-mean2)/stdPooled;
            
            % tests with the two groups
            switch testType
                case "t"
                    [testRejectNull{i,f,t},testPvals{i,f,t}]=ttest2( ...
                    valsG1,valsG2);

                case "mwu"
                    [testPvals{i,f,t},testRejectNull{i,f,t}]=ranksum( ...
                    valsG1,valsG2);
                        
                case "ks"
                    [testRejectNull{i,f,t},testPvals{i,f,t}]=kstest2( ...
                        valsG1,valsG2);
                otherwise
                    error("Incorrect <testType>, must be t, mwu or ks")
            end

            % Switch case 
            switch whenToPlot
                % Plot regardless of the significant value 
                case "allBoxes" 
                    plotBool=1;
                % Plot only if rejected
                case "ifSig" 
                    plotBool=testRejectNull{i,f,t};
                % Plot if significant and the region is not significant 
                % between HC1 and HC2 (i.e. the region remains stable
                % between the two HC session) 
                case "removeHCSig"
                    % If we are not comparing HC1vsHC2 or PD1vsPD2
                    if(t>2)
                        % if significant and doesn't differ significantly
                        % between HC1/2
                        plotBool=testRejectNull{i,f,t} && ~testRejectNull{i,f,1};
                    else
                        plotBool=testRejectNull{i,f,t};
                    end
                otherwise
                    error("Variable <whenToPlot> is incorrectly specified")
            end
                
            if(plotBool)
            % Can set a tougher plotting condition for the p value
            %if(testPvals{i,f,t} <= 0.01)
                nSignificant=nSignificant+1;

                % Name of xlabel for significant region
                switch nFreqGroups
                    case 1
                        boxNames{end+1}="eigval "+regions(i);
                    otherwise
                        boxNames{end+1}=freqStrings(f)+"_"+regions(i);
                end
                significantColors{end+1}=validatecolor(cortexColors(i));

                % Stars for the p-value to show significance in the plot
                %if(testPvals{i,f,t}<0.0001)
                %    stars="****";
                if(testPvals{i,f,t}<0.001)
                    stars="***";
                elseif(testPvals{i,f,t}<0.01)
                    stars="**";
                else
                    stars="*";
                end
                boxPvals{end+1}=stars;

                % Stars for the effect's value
                if(effectSizes{i,f,t}>=0.8)
                    stars="***";
                elseif(effectSizes{i,f,t}>=0.5)
                    stars="**";
                elseif(effectSizes{i,f,t}>=0.2)
                    stars="*";
                else
                    stars="";
                end
                boxEffects{end+1}=stars;
                % Values for the region and freq.
                x = [valsG1.'; valsG2.'];
                % Group 1 vs Group 2
                g = [zeros(n1, 1); ones(n2, 1)];
                
                % Create boxplot
                h=boxplot(x,g, 'Color',[1,0,0;0,0,1],'positions',[boxPos,boxPos+1]);
                set(h,{'linew'},{2})
                boxPos=boxPos+boxDist;
            end

            % % TODO: Maybe remove the else statement
            % else %Plot all boxes
            %     nDecimals=4; % Changeable rounding value for plotting
            %     switch nFreqGroups
            %         case 1
            %             boxNames{end+1}="eigvalue " + regions(i);
            %         case 3
            %             boxNames{end+1}=freqStrings(f)+"_"+regions(i);
            %     end
            % 
            %     boxPvals{end+1}=round(testPvals{i,f,t},nDecimals);
            %     boxEffects{end+1}=round(effectSizes{i,f,t},nDecimals);
            %     x = [valsG1.'; valsG2.'];
            %     g = [zeros(n1, 1); ones(n2, 1)];
            % 
            %     % Create boxplot
            %     h=boxplot(x,g, 'Color',[1,0,0;0,0,1],'positions',[boxPos,boxPos+1]);
            %     set(h,{'linew'},{2})
            %     boxPos=boxPos+boxDist;
            % end
        end
    end
    
    % Display x,y-axes and labels
    fSize=10; %font size
    yLim = get(gca,'YLim');
    
    % The last x-label position value +0.5
    % ----|-|-------|-|-------*-|-|----------->
    %     1 2       5 6       8 9 10
    % boxDist*(length(boxNames)-1) lands us at *
    % +2 to get to the x-val for last box
    maxBox=boxDist*(length(boxNames)-1)+2;
    
    for i=1:length(boxNames)
        xVal=1.52+(4*(i-1));
        % Plot p value
        text(xVal,-yLim(2)*0.055,boxPvals(i),'FontSize',fSize)
        % Plot the effect size
        text(xVal,-yLim(2)*0.08,boxEffects(i),'FontSize',fSize)
        % Plot region Color
        if(nFreqGroups > 1)
            text(xVal,-yLim(2)*0.1,'#','FontSize',fSize,'Color',significantColors{i},'fontweight', 'bold')
        end
    end
    
    % Plot the region name
    xlabel("nSignificant= " + string(nSignificant))

    %set(gca,'defaultTextInterpreter','none'); % No Latex interpreter

    set(gca,'xtick',1.5:boxDist:maxBox,'xticklabel',boxNames,'FontSize',fSize)
    
    % Create legend with colored squares
    L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
    set(L, {'MarkerEdgeColor'}, {[1 0 0];[0 0 1]},...
        {'MarkerFaceColor'},{[1 0 0];[0 0 1]},... % setting the markers to filled squares
        'Marker','s'); 
    lgd=legend(L, {tableStrs(t1),tableStrs(t2)},'Location', 'northeast');
    fontsize(lgd,14,'points')
    set(lgd,'color','none')
    set(figTest, 'Visible', 'on'); 
    set(gca,'YLim', [-yLim(2)*0.1 yLim(2)]);
    if(savePlot)
        saveas(figTest,testFolder+"/"+boxStr+testType+"Test_"+saveStr+"_"+tableStrs(t1)+"vs"+tableStrs(t2)+".png")
    end
end

% Areas for one hemisphere
regionsArea=[repmat("Visual",5,1);repmat("Sensorimotor",4,1);...
               repmat("Auditory",3,1);repmat("Temporal",2,1);...
               repmat("Posterior",4,1);repmat("Anterior",4,1)];
regionsArea=[regionsArea;regionsArea]; % Both hemispheres

% If we have Low, Medium and High frequencies then we need all nROI*3
switch nFreqGroups
    case 2 % 50/50 split
        regionsArea=[regionsArea;regionsArea];
    case 3 
        regionsArea=[regionsArea;regionsArea;regionsArea];
    otherwise
end

%%% Plot all p-values from the statistical test %%%
pValsMat=cell2mat(testPvals);

for t=1:nTests
    t1=testInds{t}(1);
    t2=testInds{t}(2);

    pValues=reshape(pValsMat(:,:,t),1,[]).';
    figPvals=figure;
    figPvals.Position=[100 100 2000 700];
    scatterHist = scatterhist(1:length(pValues),pValues,'Location','SouthEast','NBins',[100,100],'Parent',uipanel('Parent',figPvals));
    delete(scatterHist(2))
    sgtitle("All p-values "+saveStr+" "+tableStrs(t1)+"vs"+tableStrs(t2))
    yline(0.05, 'b')
    set(gca,'xtick',1:length(regionsArea),'xticklabel',regionsFreqs,'TickLabelInterpreter','none')

    % Plot the color of the region as a colored'#'
    if(nFreqGroups > 1)
        yLim = get(gca,'YLim');
        for i=1:nROI*nFreqGroups
            text(i-0.5,yLim(1),'#','FontSize',10,'Color',cortexColors(i),'fontweight', 'bold')
        end
    end

    if(savePlot)
        warning('off','all')
        % Generates a warning due to scatterhist overwritting the figure
        saveas(figPvals,testFolder+"/"+testType+"AllPVals_"+saveStr+"_"+tableStrs(t1)+"vs"+tableStrs(t2)+".png")
        warning('on','all')
    end
end

% %%% Might be useful if we want to look at all tests %%%
% % Create a table of the pvalues for each region and saves the significant ones.
% 
% %tTable=table; % Create new empty table
% %tTable.Region=regionsFreqs;
% 
% testTable=table; % Create new empty table
% testTable.areas=regionsArea;
% testTable.region=regionsFreqs;
% 
% for t=1:nTests
%     t1=testInds{t}(1);
%     t2=testInds{t}(2);
%     %%disp("Test: "+tableStrs(t1)+"vs"+tableStrs(t2))
% 
%     % Reshape (44,3) to (132,1) column by column
%     % t-test table
%     %tTable.tTest=reshape(tTestPvals(:,:,t),1,[]).';%tTestPvals(:,f,t);
%     %tTable.RejectNull=reshape(tTestRejectNull(:,:,t),1,[]).';%tTestRejectNull(:,f,t);
% 
%     % test table
%     testTable.test=reshape(testPvals(:,:,t),1,[]).';%testPvals(:,f,t);
%     testTable.rejectNull=reshape(testRejectNull(:,:,t),1,[]).';%testRejectNull(:,f,t);
% 
%     %disp("tTest");
%     %tTable{cell2mat(tTable{:,3})==1,:}
% 
%     rejectTable=testTable{cell2mat(testTable{:,4})==1,:};
% 
%     writematrix(rejectTable(:,1:2),testFolder+"/"+testType+"Table_"+saveStr+"_"+tableStrs(t1)+"vs"+tableStrs(t2)+".txt");
% end

end