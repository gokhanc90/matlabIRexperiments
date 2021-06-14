load('Features.mat');
load('CW09B.mat');
load('CW12B.mat');
load('GOV2.mat');
load('MQ07.mat');
load('MQ08.mat');
load('MQ09.mat');
load('NTCIR.mat');
load('WSJ.mat');
load('FeaturesTerm.mat');
load('CollectionStats.mat');

%TWS={'BM25' 'LGD' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'PL2'};
%MEASURES={'MAP' 'NDCG100' 'NDCG20'};
TWS={'BM25'};
MEASURES={'NDCG20'};
COLLECTIONS={ 'CW09B' 'CW12B' 'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'  };

gca=figure();
t = tiledlayout(4,2,'TileSpacing','none','Padding','compact');
for tw = 1:size(TWS,2)
        for measure = 1:size(MEASURES,2)
            for coll = 1:size(COLLECTIONS,2)


                dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
                runtopic = eval(dataNameR);
                
                runNoStem = runtopic.NoStem;
                runKStem = runtopic.KStem;
                
                diff = runKStem-runNoStem;
                sortedDiff = sort(diff);
                
                numberOfStemGreater = sum(sortedDiff>0);
                perStem = numberOfStemGreater*100/size(runNoStem,1);
                perStem=round(perStem);
                numberOfNoStemGreater = sum(sortedDiff<0);
                perNoStem = numberOfNoStemGreater*100/size(runNoStem,1);
                perNoStem = round(perNoStem);
                
                perTie= 100-perNoStem-perStem;
                
                nexttile
                ax=bar(sortedDiff);
                ylim([-0.6 0.6])
                yticks(-0.6:0.2:0.6)
                
                text(0.45,0.85,COLLECTIONS{coll},'Units','normalized')
                
                text(0.45,0.60,[num2str(perTie),'%'],'Units','normalized')
                
                text(0.05,0.25,'NoStem > Stem','Units','normalized')
                text(0.05,0.60,[num2str(perNoStem),'%'],'Units','normalized')
                
                text(0.8,0.75,'Stem > NoStem','Units','normalized')
                text(0.9,0.40,[num2str(perStem),'%'],'Units','normalized')
            end
        end
end
xlabel(t,'Number of Queries','FontSize',15)
ylabel(t,'Diff. in nDCG@20','FontSize',15)
set(findall(gca,'-property','FontSize'),'FontSize',14)
