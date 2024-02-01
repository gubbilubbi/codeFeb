%% TODO
% Remove comments that are not used
%%

% Plot the individual param b in a*exp(b*x) fitted to the
% autocorrelation of every patient for each group.
% Used to get a qualitative estimate of how much noise each signal has

% Load the data
megACF=load("./gsp_data/dataFeb/misc/acfVals.mat").acfVals;

%regs=[2,7,8,10,13,16,19,23,24,25,26,27,30,36,41,44];
maxLag=size(megACF,1);
nROI=size(megACF,2);
nSubjs=size(megACF,3);

%set(0,'DefaultFigureVisible','off')
whiteNoiseACF=autocorr(randn(maxLag,1),maxLag);

% Fit exponential with 2 params to the megACF
fitParams=zeros(2,nROI,nSubjs);

for s=1:nSubjs
    for r=1:nROI

        % Don't fit to the first timepoint which is always = 1
        % x_vals=zeros(maxLag-1,1);
        % x_vals(:,1)=2:maxLag;
        % pop=fit(x_vals,megACF(2:end,r,s),'exp1');
        
        % %Fit to all acf values, including the first
        x_vals=zeros(maxLag,1);
        x_vals(:,1)=1:maxLag;
        pop=fit(x_vals,megACF(:,r,s),'exp1');

        % f=figure;
        % f.Position=[100 100 900 900];
        % sgtitle("Subject="+s+", Region="+r)
        % plot(megACF(2:end,r,s));
        % hold on
        % plot(whiteNoiseACF,'k');
        % hold on
        % plot(pop);
        % hold off
        % text(5,0.9,string(pop.a)+"exp("+string(pop.b+"*x)"),"color",[1,0,0],'Fontsize',16)
        
        % disp(pop)
        % popWhiteNoise=fit(x_vals,whiteNoiseACF(1:end,1),'exp1');
        % disp(popWhiteNoise)
        % disp("_------_")

        % Store the params
        fitParams(1,r,s)=pop.a; 
        fitParams(2,r,s)=pop.b;

        %saveas(f,savePlotsFolder+"/acf_subject"+subjs(s)+"_region"+reg+".png")
    end
end

maxAxes=zeros(44,1);
minAxes=ones(44,1)*min(fitParams,[],"all");
for i=1:nGroups
    f=figure;
    f.Position=[100 100 1000 1000];
    for s=indsPDHC{i}
        spider_plot(fitParams(2,:,s),"AxesLimits",[minAxes,maxAxes]');
        hold on
    end
    %saveas(f,"./scilife/data08-12-2023/misc/corr/acfValsGroup"+i+".png")
end

popWhiteNoise=fit(x_vals,whiteNoiseACF(1:end,1),'exp1');
disp(popWhiteNoise.b); %Values for fitting exp1 to whiteNoise

function x = autocorr(A,maxLag)

    [row,col] = size(A);
    if (row ~= 1 && col ~= 1)
        error('The input should be a vector, not a matrix!');
    end
    if row == 1
        A = A';
    end
    
    N = length(A);
    
    % If no maxLag is defined then set it to the length of the timeseries
    if(~maxLag)
        maxLag=N;
    end
    x = zeros(maxLag,1);
    x(1) = sum(A.*A);
    
    for ii = 2:maxLag
        B = circshift(A,-(ii-1));
        B = B(1:(N-ii+1));
        x(ii) = sum(B.*(A(1:(N-ii+1))));
    end
    x = x/x(1);
end