function [FC,meanFC,pvalsMat]=create_connectivity_matrix(data,connType,indsPDHC,boolGroupFC, ...
                                       FCtoAdj,boolFDR,boolParCorr,boolPrec)
    % data: time x nROI x nSubjs 
    % connType: String with text for type of FC to be calculated
    % indsPDHC: cell with 4 entires, each has the data rows for each session
    % boolGroupFC: (0/1) If each session/group should have the same FC
    %              if 1 then the FC is the mean over the entire group.
    % FCtoAdj: char-string of which conversion FC -> Adj to use
    % boolFDR: Whether or not to apply false discovery rate
    % boolPrec: Whether or not to calculate the precision matrix

    % FC: Functional Connectivity matrices
    % The adjacency matrix
    % pvalsMat: The correlation values from the correlation calculation 
    % corrMeans, corrStds are mean and std of the dot product timeseries
    % returned from create_corr_timeseries

    nROI=size(data,2);
    nSubjs=size(data,3);
    
    FC=zeros(nROI,nROI,nSubjs);
    meanFC=zeros(nROI,nROI,nSubjs);
    pvalsMat=zeros(nROI,nROI,nSubjs);

    switch connType

        % Correlation
        case 'corr'
            for s=1:nSubjs
                % Correlation between regions over time
                subData=data(:,:,s);

                if boolParCorr % Calculate partial correlations
                    [rho,pvals]=partialcorr(subData);
                else
                    [rho,pvals]=corrcoef(subData);

                    if(boolPrec) % Calculate Precision matrix
                        %FC(:,:,s)=inv(rho);
                        rho=inv(rho);
                        FC(:,:,s)=rho-diag(diag(rho));

                        continue;
                    end
                end
                
                pvalsMat(:,:,s)=pvals;

                % pValues=reshape(pvals,1,[]).';
                % figPvals=figure;
                % figPvals.Position=[100 100 2000 700];
                % scatterHist = scatterhist(1:length(pValues),pValues,'Location','SouthEast','NBins',[100,100],'Parent',uipanel('Parent',figPvals));
                % delete(scatterHist(2))

                FC(:,:,s)=rho-diag(diag(rho)); %set diag to 0

                % % non-significant entries are set to 0
                % TODO make more beautiful
                if(boolFDR)
                    tempFC=FC(:,:,s);
                    tempFC(pvals > 0.05)=0; 
                    FC(:,:,s)=tempFC;
                end

                % Check that matrix is symmetric
                if(~issymmetric(FC(:,:,s)))
                    error('Invalid Non-symmetric FC')
                end
                
                % TODO comment out?
                % nNonSig=sum(pvals > 0.05,"all");
                % if(nNonSig > 0)
                %     disp("Sub: "+string(s)+",FC N Non significant: "+string(nNonSig))
                %     disp("% : " + string(nNonSig/nROI^2))
                % end

            end

            %All subject in a group receive the same group mean connectivity matrix
            %if(boolGroupFC)
            for g=1:length(indsPDHC)
                % Calculate mean FC for every group
                groupFC=mean(FC(:,:,indsPDHC{g}),3);
                % The FC for all subjects in the group is set to the mean FC.
                for sub=indsPDHC{g}
                    meanFC(:,:,sub)=groupFC;
                end
            end
            %end
        
        % Correlation over all of the subjects
        case 'Avgcorr'
            for i=1:length(indsPDHC)
                % Correlation between regions over time
                groupData=data(:,:,indsPDHC{i});

                % Stack all of the timeseries in each ROI (column) 
                % for each group
                groupROIData=reshape(permute(groupData,[1,3,2]),size(groupData,1)*size(groupData,3),size(groupData,2));
                
                if boolParCorr % partial correlations
                    [rho,pvals]=partialcorr(groupROIData);
                else % normal corr and or precision matrix
                    [rho,pvals]=corrcoef(groupROIData);

%                     if(boolPrec) % Calculate Precision matrix
%                         rho=inv(rho);
%                         FC(:,:,s)=rho-diag(diag(rho));
%                         
%                         continue;
%                     end
                end
                
                for s=indsPDHC{i}
                    pvalsMat(:,:,s)=pvals;
                    
                    % pValues=reshape(pvals,1,[]).';
                    % figPvals=figure;
                    % figPvals.Position=[100 100 2000 700];
                    % scatterHist = scatterhist(1:length(pValues),pValues,'Location','SouthEast','NBins',[100,100],'Parent',uipanel('Parent',figPvals));
                    % delete(scatterHist(2))

                    FC(:,:,s)=rho-diag(diag(rho)); %set diag to 0
                    
                    % % non-significant entries are set to 0
                    % TODO make more beautiful
                    if(boolFDR)
                        tempFC=FC(:,:,s);
                        tempFC(pvals > 0.05)=0; 
                        FC(:,:,s)=tempFC;
                    end

                    % Check that matrix is symmetric
                    if(~issymmetric(FC(:,:,s)))
                        error('Invalid Non-symmetric FC')
                    end
                end
                
                % TODO comment out?
                % nNonSig=sum(pvals > 0.05,"all");
                % if(nNonSig > 0)
                %     disp("Group: "+string(i))
                %     disp(nNonSig)
                %     disp("% : " + string(nNonSig/nROI^2))
                % end
        
            end

        otherwise
            error("Wrong FC computation type")
    end
% 
%     A=FC;
%     %If FC !>= 0 then the adjacency matrix is the FC but w_ij >= 0
%     switch FCtoAdj
%         case 'abs'
%             A=abs(FC);
%         case 'geq0' % All values below 0 set to 0
%             A(FC<0)=0;
%         otherwise
%             error('No FC to Adj conversion')
%     end
% 
%     for s=1:nSubjs
%         if(~issymmetric(A(:,:,s)))
%             error("Adjacency not symmetric")
%         end
% 
%         g = digraph(A(:,:,s));
%         bins = conncomp(g, 'Type', 'weak');
%         isConnected = all(bins == 1);
%         if(~isConnected)
%             error('Adjacency not connected')
%         end
%     end
end