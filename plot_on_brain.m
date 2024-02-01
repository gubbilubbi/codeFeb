%% Author: Valter Lundegårdh

function plot_on_brain(data,std,CodeBook,CA,annot,caseTHR,THR,Exp_Sphere,bmeshAlpha)
    % data: nROIxnGroups
    % std of the data:  nROIxnGroups
    % Codebook: Struct with information about the data needed for plotting.
    % CA: 0 min and max value of the data in an array
    % annot: array of strings for the annotations of each row.
    % boolTHR: Whether or not to apply threshold to colors, i.e. set
    %          all colors less than half of the max to 0.
    % THR: Threshold value if boolTHR is 1

    % TODO: Should I do log2? Not for eigenvectors according to Huang2016
    % This is good if I want SDI though...
    %lowCC2=log2(avgLowU); 
    %highCC2=log2(avgHighU); 

    % CC2 are the values we want to plot
    % CodeBook is a cell with:ids,nROI,name,sname and center coords
    N=CodeBook.num;    % number of regions in the atlas
    CM=zeros(N,N); % For plotting connections

    % TODO Should I threshold?
    % %% adjust Cvalues for saturation (to eliminate outliers peaks)
    % saturate=1;
    % if saturate
    % thr=1;
    % CC2new=CC2;
    % CC2new(find(CC2>thr))=0;
    % CC2new(find(CC2>thr))=max(CC2new);
    % CC2new(find(CC2<-thr))=0;
    % CC2new(find(CC2<-thr))=min(CC2new);
    % CC2=CC2new;
    % end

    % Min max for the current data
    %minData=min(data,[],"all");
    maxData=max(data,[],"all");
    %disp("min data")
    %disp(minData)
    %disp("max data")
    %disp(maxData)
    
    % If no colormap specified, calculate min and max
    % CA=[-1,1];
    if(CA==0)
        CA=[0,maxData];
    end

    if(~(std==0)) %If variance is not 0 calc max variance
        maxStd=max(std,[],"all");
        disp("Max standard deviation");
        disp(maxStd);
    end
    
    Views=[1,3,2]; % Brain view (axial, coronal, sagittal)
    nData=size(data,2);

    switch nData
        % For one session
        case 1
            poss={[-0.04 0.35 0.43 0.43],[0.28 0.35 0.36 0.36],[0.62 0.35 0.36 0.36]};
        % For three sessions
        case 3
            poss={
            [-0.05 0.66 0.33 0.33],[0.19 0.66 0.32 0.32],[0.49 0.66 0.32 0.32],...
            [-0.05 0.33 0.33 0.33],[0.19 0.33 0.32 0.32],[0.49 0.33 0.32 0.32],...
            [-0.05 0 0.33 0.33],[0.19 0 0.32 0.32],[0.49 0 0.32 0.32]
            };
        case 4
            poss={
            [0.01 0.75 0.25 0.25],[0.225 0.75 0.23 0.23],[0.46 0.72 0.27 0.26],...
            [0.01 0.5 0.25 0.25],[0.225 0.51 0.23 0.23],[0.46 0.48 0.27 0.26],...
            [0.01 0.25 0.25 0.25],[0.225 0.26 0.23 0.23],[0.46 0.23 0.27 0.26],...
            [0.01 0.0 0.25 0.25],[0.225 0.01 0.23 0.23],[0.46 -0.02 0.27 0.26]
            };
        otherwise
            error("Wrong data size")
    end
    
    Colormap_nodes='jet';
    
    %For plotting connections
    T_conn=0; 
    Colormap_edges='jet';
    Gamma=0.5;
    LinearWeight=1;

    Factor_SphereSize=maxData;
    Factor_Col=maxData;
    
    % plot with normal color scheme 
    for ses=1:size(data,2)
        %CC=abs(CC2(:,ses)); %Only makes sense if we have negative values
        CC=data(:,ses);
        
        if(std==0)
            %Exp_Sphere=2; %Size of spheres exponential
            % Set to 1 means all regions have the same size
            sphereSize=1;%CC./Factor_SphereSize; %1;
        else
            %Exp_Sphere=2; %Size of spheres exponential
            sphereSize=std(:,ses)./maxStd;
        end
        
        sphereColor=data(:,ses);%./Factor_Col; %1;

        switch caseTHR

            case 0 % no threshold 
               
            case 1 %less than threshold set to 0
                sphereColor(sphereColor < CA(2)*THR)=0;

            case 2 %greater than treshold set to 0
                sphereColor(sphereColor >= THR)=1;
            otherwise % No threshold
                error("caseTHR has unknown value")
            % Filter data less than 80% of max value
            
        end
    
        for i=1:length(Views)
            %disp(poss{(ses-1)*3+i})
            axes('Position',poss{(ses-1)*3+i})
            myPlotBrainGraph(CM,CC,data(:,ses),CodeBook,T_conn,sphereSize,...
                sphereColor,Exp_Sphere,Views(i),Colormap_nodes,Colormap_edges,Gamma,...
                LinearWeight,CA,bmeshAlpha);
        end
    end
    ax = axes;
    switch nData
        case 1
            colorbar('SouthOutside');
        case 3
            colorbar(ax);
            annotation('textbox', [0.95, 0.88, 0, 0], 'string', annot(1))
            annotation('textbox', [0.95, 0.55, 0, 0], 'string', annot(2))
            annotation('textbox', [0.95, 0.23, 0, 0], 'string', annot(3))
        case 4
            colorbar(ax)
            annotation('textbox', [0.95, 0.88, 0, 0], 'string', annot(1))
            annotation('textbox', [0.95, 0.65, 0, 0], 'string', annot(2))
            annotation('textbox', [0.95, 0.4, 0, 0], 'string', annot(3))
            annotation('textbox', [0.95, 0.15, 0, 0], 'string', annot(4))
        otherwise
    end    
    clim(CA);
    ax.Visible = 'off';
end