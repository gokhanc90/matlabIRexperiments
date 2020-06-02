load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/Features.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CW09B.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CW12B.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/GOV2.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ07.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ08.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ09.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/NTCIR.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/WSJ.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/FeaturesTerm.mat')
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CollectionStats.mat')


STEMMERS={'SnowballEng'};
%STEMMERS={'KStem'};
%TWS={'BM25' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'LGD' 'PL2'};
TWS={ 'BM25'};
%MEASURES={'MAP' 'NDCG20' 'NDCG100'};
MEASURES={'NDCG20'};
COLLECTIONS={ 'CW09B' 'CW12B' 'NTCIR' 'GOV2'  'WSJ' 'MQ07' 'MQ08' 'MQ09'};
%COLLECTIONS={ 'CW09B'};


for s = 1:size(STEMMERS,2)
   for tw = 1:size(TWS,2)
        for measure = 1:size(MEASURES,2)
            for coll = 1:size(COLLECTIONS,2)
                    
                    docCount = CollectionStats.(strcat(COLLECTIONS{coll},'DocCount'));
                    TF = CollectionStats.(strcat(COLLECTIONS{coll},'TermCount'));
                    
                    dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
                    runtopic = eval(dataNameR);

                    
                    dataNameF = strcat('Feature',COLLECTIONS{coll},STEMMERS{s});
                    features = eval(dataNameF);

                    dataNameFTerms = strcat(dataNameF,'Term');    
                    terms=eval(dataNameFTerms);
                                
                    Joined = JoinTables(terms,features,runtopic);
                            if strcmp(COLLECTIONS{coll},'WSJ')
                                JoinedT=Joined;
                                JoinedT(:,:)=[];
                            else
                                JoinedT = extTrainData(COLLECTIONS{coll},MEASURES{measure},TWS{tw},STEMMERS{s});
                            end

                    Joined.Diff = Joined.NoStem - Joined.(STEMMERS{s});
                    JoinedT.Diff = JoinedT.NoStem - JoinedT.(STEMMERS{s});
                            
                    [p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.(STEMMERS{s}));
                    
                          [pT,isSigT,oracleT,labelT]=getOracle(JoinedT.NoStem,JoinedT.(STEMMERS{s}));
                          ScoresT=JoinedT(:,{'NoStem',STEMMERS{s}});

                    fileID = fopen('runtopic2.txt','a');

                    [m, n]=size(Joined);
                            predictionScores=zeros(m,1);
                            predictedlabel=categorical(zeros(m,1));
                            for i=1:m
                                Xtest=Joined(i,:);
                                Xtrain=Joined([1:i-1,i+1:end],:);
    %                            Xtrain=vertcat(Xtrain,JoinedT);
                                
                                Xtrain=Xtrain(Xtrain.Diff~=0,:);

                                [trainedModel, validationRMSE, yfit,predictorNames] = trainLinearInteractionReg(Xtrain,Xtest);
                                if yfit > 0
                                    predictionScores(i)=Xtest.NoStem; 
                                    predictedlabel(i)='0';
                                else
                                    predictionScores(i)=Xtest.(STEMMERS{s}); 
                                    predictedlabel(i)='1';
                                end
                            end
                            Scores=Joined(:,{'NoStem',STEMMERS{s}});
                            [ms, significant, m1, m2, oracle] = AverageNDCG(table2array(Scores),predictedlabel);
                            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s %s\n','InteractionReg',...
                                ms,significant,m1,m2,oracle,dataNameR,dataNameF,strjoin(predictorNames));
 
                
            end
        end
   end
end

function Joined = JoinTables(terms,features,runtopic)
    terms.IdfAdvDF=terms.idfs .* terms.advanceDF;
    terms.CtiAdvDF=terms.ctis .* terms.advanceDF;
    terms.idfRatio = terms.idfStem ./ terms.idfs;
    termAgg1 = groupsummary(terms,'QueryID',{'mean','max','sum','var','min'},'IdfAdvDF');
    termAgg2 = groupsummary(terms,'QueryID',{'mean','max','sum','var'},'CtiAdvDF');
    termAgg3 = groupsummary(terms,'QueryID',{'mean','max','sum','var','min'},'idfRatio');
    termAgg4 = groupsummary(terms,'QueryID',{'mean','max','var'},'idfStem');

    gamma1Table = Gamma1(terms(:,{'QueryID','word','idfStem'}));

    Joined = join(features,runtopic,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg1,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg2,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg3,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg4,'LeftKeys',1,'RightKeys',1);
    Joined.idfRatioMinMax = Joined.min_idfRatio ./ Joined.max_idfRatio;
    
    dftfTable = chi2(terms(:,{'QueryID','DFNoStem','TFNoStem','DF','TF'}));

    corrTable = IdfOrderDist(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
    lstmstTable = LSTMST(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
    Joined = join(Joined,corrTable,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,lstmstTable,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,gamma1Table,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,dftfTable,'LeftKeys',1,'RightKeys',1);

    
end

%min/max
function dftfTable = chi2(DFTFs)
    QIDs = unique(DFTFs.QueryID);
    
    dftfTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','Chi2DFTF'});
    
    for i=1:size(QIDs,1)
        terms = DFTFs(DFTFs.QueryID==QIDs(i),:);
        obs1=[terms.DFNoStem' terms.TFNoStem'];
        obs2=[terms.DF' terms.TF'];
        c= chi2Val([obs1;obs2]);

        dftfTable.Chi2DFTF(dftfTable.QueryID == QIDs(i))=c;
    end
        
end


function corrTable = IdfOrderDist(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    corrTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','correlationTerms'});
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'idfs'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'idfStem'); % sort the table by 'DOB'
        
        %if(lev(termsNoStem.Position, termsStem.Position)/wordCount > 0.1) 
        c=corr(termsNoStem.Position,termsStem.Position,'Type','Spearman');
  
        corrTable.correlationTerms(corrTable.QueryID == QIDs(i))=c;
      
    end
        
end

function lstmstTable = LSTMST(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    lstmstTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','lstmstTerms'});
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount = size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'idfs','ascend'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'idfStem','ascend'); % sort the table by 'DOB'
        
        c=0;
        if(termsNoStem.Position(1) == termsStem.Position(1)) 
            c=c+0.5;
        end
        if(termsNoStem.Position(wordCount) == termsStem.Position(wordCount)) 
            c=c+0.5;
        end
        
        lstmstTable.lstmstTerms(lstmstTable.QueryID == QIDs(i))=c;
    end
        
end


%min/max
function gamma1Table = Gamma1(stemTerms)
    QIDs = unique(stemTerms.QueryID);
    
    gamma1Table = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','GammaStem'});
    
    for i=1:size(QIDs,1)
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount = size(termsStem,1);
        order=1:wordCount;
        termsStem.Position = order';
        
        termsStem = sortrows(termsStem, 'idfStem','ascend'); % sort the table by 'DOB'
        
        c=termsStem.idfStem(1)/termsStem.idfStem(wordCount);

        gamma1Table.GammaStem(gamma1Table.QueryID == QIDs(i))=c;
    end
        
end


function JoinedT = extTrainData(coll,measure,tw,stemmer)

    if strcmp(coll, 'CW09B')
        dataNameRT = strcat('MQ09','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ09',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);
    elseif strcmp(coll, 'MQ09')
        dataNameRT = strcat('CW09B','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','CW09B',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);
    elseif strcmp(coll, 'CW12B')
        dataNameRT = strcat('NTCIR','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','NTCIR',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);
    elseif strcmp(coll, 'NTCIR')
        dataNameRT = strcat('CW12B','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','CW12B',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);
    elseif strcmp(coll, 'GOV2')
        dataNameRT = strcat('MQ07','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ07',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);

        dataNameRT = strcat('MQ08','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ08',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT);
        
        JoinedT=vertcat(JoinedT,JoinedT2);
    elseif strcmp(coll, 'MQ07')
        dataNameRT = strcat('GOV2','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','GOV2',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);

        dataNameRT = strcat('MQ08','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ08',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT);
        JoinedT=vertcat(JoinedT,JoinedT2);
    elseif strcmp(coll,'MQ08')
        dataNameRT = strcat('GOV2','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','GOV2',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT);

        dataNameRT = strcat('MQ07','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ07',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT);
        JoinedT=vertcat(JoinedT,JoinedT2);
    end
                    
end


