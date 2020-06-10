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

TWS={'BM25' 'LGD' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'PL2'};
MEASURES={'MAP' 'NDCG100' 'NDCG20'};
COLLECTIONS={ 'CW09B' 'CW12B' 'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'  };

fileID = fopen('oracleSignificant.txt','w');

for tw = 1:size(TWS,2)
        for measure = 1:size(MEASURES,2)
            for coll = 1:size(COLLECTIONS,2)


                dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
                runtopic = eval(dataNameR);
                
                stemmers= {'KStemQBS'    'SnowballEngQBS'    'HPS'    'Gupta'    'KStem'    'SnowballEng'};
                for i=1:6
                    [p,isSig,oracle,label] = getOracle(runtopic.NoStem, runtopic.(stemmers{i}) );
                    fprintf(fileID,'%s\t%.4f\t%.2e\t%d\n',strcat(dataNameR,'.NoStem-', dataNameR, '.', stemmers{i}),mean(oracle),p,isSig );    
                end
                fprintf(fileID,'--------------\n');
            end
        end
end
