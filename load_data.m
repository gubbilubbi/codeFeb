function data=load_data(mypath,nROI,times,idsHC)
    
    % Parameters
    % mypath: Path to where the data is stored e.g. 'scilife/meg/';
    % times: range of start to end indx for the data
    % nROI: Number of regions
    % Which of the patients in folder order that are HC e.g. [2,4] 
    
    % Returns
    % data: Matrix time x regions x subject 
    % The Subjects in the returned data are ordered as:
    % PD first sessions, PD second sessions
    % HC first sessions, HC second sessions

    
    dataZeroFolders=dir([mypath '0*']); % Only consider folders that start with 0
    % Exclude 0354 and 0332 because no .npy        
    idx = cellfun(@(x)(contains(x,'0345') || contains(x,'0332')),{dataZeroFolders.name}, 'Unif',0);  
    dnew = dataZeroFolders(~[idx{:}]);                                         
    dataFolderNames = {dnew.name}; 
    nSubjs=length(dataFolderNames);

    % Binary matrix of which subjects in order that are HC (0) and PD (1)
    boolPD=create_bool_PD(nSubjs,idsHC);
    
    % Total number of PDs and HCs
    nPD=sum(boolPD);
    nHC=sum(~boolPD);
    
    % Counters for current PD and HC subject
    indPD=1;
    indHC=1;
    
    % data is times x nROI x nSubjs*2 because 2 sessions for each subject
    data=zeros(length(times),nROI,nSubjs*2);
    % Temporary matrices for storing PD and HC data separately
    % that are later concatenated.
    pdData=zeros(length(times),nROI,nPD*2);
    hcData=zeros(length(times),nROI,nHC*2);

    strs=["-ses1-","-ses2-"];
    for i=1:length(dataFolderNames)
        for j=1:length(strs)
            dataName=dataFolderNames{i}+strs(j)+"labeltc.npy";
            subj=readNPY(dataName)'; %time x regions    

            % if PD
            if(boolPD(i))
                % if ses 1
                if(mod(j,2))
                    pdData(:,:,indPD)=subj(times,:);
                % if ses 2
                else
                    pdData(:,:,indPD+nPD)=subj(times,:);
                    indPD=indPD+1;
                end

            % if HC
            else
                % if ses 1
                if(mod(j,2))
                    hcData(:,:,indHC)=subj(times,:);
                % if ses 2
                else
                    hcData(:,:,indHC+nHC)=subj(times,:);
                    indHC=indHC+1;
                end
            end
        end
    end
    % Stack the data PDses1, PDses2, HCses1, HCses2
    data(:,:,1:2*nPD)=pdData;
    data(:,:,2*nPD+1:end)=hcData;
    disp("Data extracted")
end


function boolPD=create_bool_PD(nSubjs,idsHC)
    boolPD=ones(nSubjs,1); % PD indices set to 1
    % HC indicies set to 0
    for row=idsHC
        boolPD(row)=0;
    end
    boolPD=logical(boolPD);
end

