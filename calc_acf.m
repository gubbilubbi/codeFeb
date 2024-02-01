function acfVals=calc_acf(meg)  
    % Takes in a timeseries and calculates the autocorrelation function
    % for a given lag (maxTime)

    % The data is assuemd to already have been z-scored

    nROI=size(meg,2);
    nSubjs=size(meg,3);
    maxLag=500; % The maxLag which can be changed
    acfVals=zeros(maxLag,nROI,nSubjs);

    for s=1:nSubjs
        for r=1:nROI
            acfVals(:,r,s)=autocorr(meg(:,r,s),maxLag);
        end
    end
end

function x = autocorr(A,maxLag)
    [row,col] = size(A);
    if (row ~= 1 && col ~= 1)
        error('The input should be a vector, not a matrix!');
    end
    if row == 1
        A = A';
    end
    
    N = length(A); % N timepoints
    
    % If no maxLag is defined then set it to the length of the timeseries
    if(~maxLag)
        maxLag=N;
    end
    x = zeros(maxLag,1);
    x(1) = sum(A.*A); % correlation of the signal with itself
    
    % For every time shift A backwards one step and calculate the
    % dot-product between itself and the unshifted signal to retrieve the
    % auto-correlation for maxLag timepoints
    for i = 2:maxLag
        B = circshift(A,-(i-1));
        B = B(1:(N-i+1));
        x(i) = sum(B.*(A(1:(N-i+1))));
    end
    % Divide by the standard deviation sig_x = (1/N)*sum((x-mean_x)^2) which
    % is the same as the autocorrelation of the same datapoint.
    x = x/x(1);
end