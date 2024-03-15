function single_psd_energy_split(U,X_hat,nScans,psdNN,psdFolder)

nROI=size(X_hat,1);
dataLength=size(X_hat,2);
nGroups=4;

nFreqs=size(psdNN,2)+1;

%%% Create 3 versions of the U matrix, containing only low-, medium- or high- frequency Laplacian eigenvectors
Ufilter=zeros([size(U),nFreqs]);
% If 33 split
if(nFreqs==3)
    for s=1:nScans
        Ufilter(:,1:psdNN(s,1),s,1)=U(:,1:psdNN(s,1),s); %low spatial frequencies  
        Ufilter(:,psdNN(s,1)+1:psdNN(s,2)-1,s,2)=U(:,psdNN(s,1)+1:psdNN(s,2)-1,s); %medium spatial frequencies
        Ufilter(:,psdNN(s,2):end,s,3)=U(:,psdNN(s,2):end,s); %high spatial frequencies
    end
% If 50 split
else
    for s=1:nScans
        Ufilter(:,1:psdNN(s,1),s,1)=U(:,1:psdNN(s,1),s); %low spatial frequencies  
        Ufilter(:,psdNN(s,1)+1:end,s,2)=U(:,psdNN(s,1)+1:end,s); %high spatial frequencies
    end
end

%%% reconstruct functional signals containing only low / high frequencies
% This both filters and inverts the data back to time domain 
% because V^-1 * H^ (H binary diag matrix) takes only the specified values
% from V^-1, i.e. low or high components.

Xfilter=zeros(dataLength,nROI,nScans,nFreqs);
energyX=zeros(nROI,nScans,nFreqs);
meanX=zeros(nROI,nGroups,nFreqs);
stdX=zeros(nROI,nGroups,nFreqs);
for f=1:nFreqs
    for s=1:nScans
        % reconstruction of coupled signal, low, medium and high freqs.
        Xfilter(:,:,s,f)=(Ufilter(:,:,s,f)*X_hat(:,:,s)).';
    end
    
    % Energy (l2-norm) over time for each region and scan
    % For SDI and statistical tests
    energyX(:,:,f)=squeeze(vecnorm(Xfilter(:,:,:,f),2,1));
    
    % For mean and std energy plots on the brain    
    for s=1:nScans
        % Mean and std of reconstructed signals over time 
        meanX(:,s,f)=mean(abs(Xfilter(:,:,s,f)),1);
        stdX(:,s,f)=std(abs(Xfilter(:,:,s,f)),0,1);
    end
end

save(psdFolder+"/meanX.mat",'meanX')
save(psdFolder+"/stdX.mat",'stdX')
save(psdFolder+"/energyX.mat",'energyX')
