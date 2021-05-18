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


STEMMERS={'KStem' };
%STEMMERS={'KStem'};
TWS={'BM25' 'LGD' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'PL2'};
%TWS={ 'BM25'};
MEASURES={'MAP' 'NDCG100' 'NDCG20'};
%MEASURES={'NDCG20' };
COLLECTIONS={ 'CW09B' 'CW12B' 'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'};
%COLLECTIONS={ 'CW09B'};


%for s = 1:size(STEMMERS,2)
   for tw = 1:size(TWS,2)
        for measure = 1:size(MEASURES,2)
            for coll = 1:size(COLLECTIONS,2)
                
                    
                   % docCount = CollectionStats.(strcat(COLLECTIONS{coll},'DocCount'));
                   % TF = CollectionStats.(strcat(COLLECTIONS{coll},'TermCount'));
                    
                    dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
                    runtopic = eval(dataNameR);
                    

                    % Import the data
                    T = readtable(strcat("/home/ubuntu/Desktop/Data/Sonuçlar/nlp_",COLLECTIONS{coll},".txt"),'Delimiter', '\t');
                    T=T(:,1:2);
                    
                   assert(size(runtopic,1)==size(T,1))
                    
                    Joined = join(runtopic,T,'LeftKeys',1,'RightKeys',1);
                    assert(size(Joined,1)==size(T,1))
                    
                    writetable(Joined,strcat(dataNameR,'.csv'),'Delimiter',',');  
                    
                         
            end
        end
   end
%end
