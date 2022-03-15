figure('Color','none','Position',[0 0 1600 500])
tl = tiledlayout(1,3,'TileSpacing','none','Padding','none');
title(tl,'GOV2','FontSize',17)
dataKStem=data(:,[1,2]);
dataLovins=data(:,[3,4]);
dataPorter=data(:,[5,6]);

nexttile
myplot(dataKStem,'',labelK)
nexttile
myplot(dataLovins,'',labelL)
nexttile
myplot(dataPorter,'',labelP)



function [] = myplot(data, t, labels)
    plot(data,'LineWidth',2)
    set(findall(gca,'-property','FontSize'),'FontSize',14)
    xlabel('alpha','FontSize',15)
    ylabel('Robustness','FontSize',15)
    title(t)
    xticks([ 1 2 3 4 5 6])
    xticklabels({'0', '1', '2', '3', '4', '5'})
    yticks([-5 -4 -3 -2 -1 0 1 2 3 4])
    yticklabels({'-5','-4', '-3', '-2', '-1' '0', '1', '2', '3', '4'})
    ylim([-4.2 5.2])
    xlim([0.8 6.2])
    legend(labels)
    grid on
    set(findall(gca,'-property','FontSize'),'FontSize',14)
end