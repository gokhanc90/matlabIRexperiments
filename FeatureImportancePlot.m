dataKStem = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureWSJKStem',[1 2 3]);
dataLovins = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureWSJLovins',[1 2 3]);
dataPorter = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureWSJSnowballEng',[1 2 3]);

dataKStem=sortrows(dataKStem,1,'descend');
dataLovins=sortrows(dataLovins,1,'descend');
dataPorter=sortrows(dataPorter,1,'descend');

fullK = table2array(KStemAll(KStemAll.Collection=="WSJ","Score"));
fullL = table2array(LovinsAll(LovinsAll.Collection=="WSJ","Score"));
fullP = table2array(PorterAll(PorterAll.Collection=="WSJ","Score"));

Row={fullK,11,'IncludeAll'};
dataKStem = [Row;dataKStem];


Row={fullL,11,'IncludeAll'};
dataLovins = [Row;dataLovins];


Row={fullP,11,'IncludeAll'};
dataPorter = [Row;dataPorter];

%RANK%
TKRank=table([11:-1:1]',dataKStem.Ablate,'VariableNames',{'Rank','Feature'});
TLRank=table([11:-1:1]',dataLovins.Ablate,'VariableNames',{'Rank','Feature'});
TPRank=table([11:-1:1]',dataPorter.Ablate,'VariableNames',{'Rank','Feature'});
TRankWSJ = [TKRank;TLRank;TPRank];

figure('Position',[0 0 1600 400])
tl = tiledlayout(1,3,'TileSpacing','none','Padding','none');
title(tl,'WSJ','FontSize',16)


nexttile
myplot(dataKStem,'Sel-KStem')
nexttile
myplot(dataLovins,'Sel-Lovins')
nexttile
myplot(dataPorter,'Sel-Porter')

saveas(gcf,'WSJFeatureImportance.fig');
saveas(gcf,'WSJFeatureImportance.eps','epsc');

function [] = myplot(data, t)
    barh(data.Mean)
    xlim([0.32 0.38])
    yticklabels(data.Ablate)
    title(t)
    ax = gca;
    ax.FontSize = 13;
    xlabel('Score','FontSize',15)
    ylabel('Removed Feature','FontSize',15)
    grid on
   % set(findall(gca,'-property','FontSize'),'FontSize',14)
end

