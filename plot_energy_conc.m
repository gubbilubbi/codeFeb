function plot_energy_conc(x,stdX,colors,plotLegends,subPlotsInds,sgtit,fileName,boolSave)

% x: matrix nROIxnGroups with mean value of the energy
% stdX: matrix nROIxnGroups with std value of the energy
% colors: string array with the colors for each nROI
% plotLegends: string array with the titles for subplots
% subPlotsInds: cell with inds for which corner each g should be plotted
% sgtit: Main title string
% fileName: string to be saved as
% boolSave: Boolean value which decides if the file is to be saved

nROI=size(x,1);
nGroups=size(x,2);

f=figure;
f.Position = [300 300 700 700];
sgtitle(sgtit)

for g=1:nGroups
    for i=1:nROI
        subplot(2,2,subPlotsInds{g})
        if(size(colors,2)>1) % TODO not sure why I made this one...
                             % is it even used?
            clrs=colors(i,g);
        else 
            clrs=colors(i);
        end
        barh(i,x(i,g),clrs);
        set(gca, 'YDir','reverse') % swaps lowest <-> highest on y-axis
        hold on
        errorbar(x(i,g),i,stdX(i,g),'k.','horizontal')
        title(plotLegends{g})
        xlim([0,max(x,[],"all")+max(stdX,[],"all")]);
        xlabel("mean l2 norm")
    end
end
qw{1} = plot(nan, 'b');
qw{2} = plot(nan, 'r');
qw{3} = plot(nan, 'y');
qw{4} = plot(nan, 'g');
qw{5} = plot(nan, 'c');
qw{6} = plot(nan, 'm');
legend([qw{:}], {'Visual','Sensorimotor','Auditory', 'Temporal','Posterior','Anterior'}, 'Position', [0.47 0.444 0.1 0.1]);

if(boolSave)
    saveas(f,fileName)
end


end