load('Features.mat');
load('CW09B.mat');
load('CW12B.mat');
load('GOV2.mat');
load('MQ07.mat');
load('MQ08.mat');
load('MQ09.mat');
load('NTCIR.mat');
load('WSJ.mat');
load('FeaturesTerm.mat')
load('CollectionStats.mat')


STEMMERS={'KStem' 'SnowbalEng' 'NoStem' 'HPS' 'Gupta' 'KStemQBS' 'SnowbalEngQBS'};
%STEMMERS={'KStem'};
%TWS={'BM25' 'LGD' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'PL2'};
TWS={ 'BM25'};
MEASURES={'MAP' 'NDCG100' 'NDCG20'};
%MEASURES={'NDCG20' };
COLLECTIONS={ 'CW09B' 'CW12B' 'GOV2' 'MQ07' 'MQ08' 'MQ09' 'NTCIR' 'WSJ'};
%COLLECTIONS={ 'CW09B'};



for tw = 1:size(TWS,2)
    for measure = 1:size(MEASURES,2)
        for coll = 1:size(COLLECTIONS,2)
            dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
            runtopic = eval(dataNameR);

            table = runtopic(:,{'KStem' 'SnowballEng' 'NoStem' 'HPS' 'Gupta' 'KStemQBS' 'SnowballEngQBS'});
            avg = mean(table2array(table));
            NoStemIndex = find(table.Properties.VariableNames == "NoStem");
            largerInd = find(avg>avg(NoStemIndex));
            table2 = table(:,largerInd);
            [p,tbl] = anova1(table2array(table2),[],'off');
            fprintf("%s %.4f %d\n",dataNameR,p, tbl{2,3})
        end
    end
end
   
fprintf("\n\n");

for tw = 1:size(TWS,2)
    for measure = 1:size(MEASURES,2)
        for coll = 1:size(COLLECTIONS,2)
            dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
            runtopic = eval(dataNameR);

            table = runtopic(:,{'KStem' 'SnowballEng' 'NoStem' 'HPS' 'Gupta' 'KStemQBS' 'SnowballEngQBS'});
            avg = mean(table2array(table));
            NoStemIndex = find(table.Properties.VariableNames == "NoStem");
            largerInd = find(avg>avg(NoStemIndex));
            table2 = table(:,largerInd);
           
                [p,tbl] = friedman(table2array(table2),1,'off');
                fprintf("%s %.4f %d\n",dataNameR,p, tbl{2,3})

            
        end
    end
end

   