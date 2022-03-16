dataKStem = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureMQ09KStem',[1 2 3]);
dataLovins = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureMQ09Lovins',[1 2 3]);
dataPorter = featureswithLovinsS2(featureswithLovinsS2.Stemmer=='FeatureMQ09SnowballEng',[1 2 3]);

dataKStem=sortrows(dataKStem,1,'descend');
dataLovins=sortrows(dataLovins,1,'descend');
dataPorter=sortrows(dataPorter,1,'descend');


figure('Position',[0 0 1600 400])
tl = tiledlayout(1,3,'TileSpacing','none','Padding','none');
title(tl,'MQ09','FontSize',16)


nexttile
myplot(dataKStem,'Sel-KStem')
nexttile
myplot(dataLovins,'Sel-Lovins')
nexttile
myplot(dataPorter,'Sel-Porter')

function [] = myplot(data, t)
    barh(data.Mean)
    xlim([0.1 0.3])
    yticklabels(data.Ablate)
    title(t)
    
    xlabel('Score','FontSize',13)
    ylabel('Removed Feature','FontSize',13)
    grid on
   % set(findall(gca,'-property','FontSize'),'FontSize',14)
end