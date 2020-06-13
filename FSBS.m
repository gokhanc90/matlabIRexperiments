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
MEASURES={'NDCG20'  'NDCG100' 'MAP'};
COLLECTIONS={ 'CW09B' 'CW12B' 'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'  };

fileID = fopen('oracleSignificant.txt','w');

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
                    
                    Joined = JoinTables(terms,features,runtopic,TF);
                    
                    SelectedFeatures=Joined(:,{... 
                    'Gamma','Omega','AvgPMI','MaxPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF','MaxIDF','MeanCTI',... %11
                    'VarCTI','MaxCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',... %19
                    'MaxSCQ','SumSCQ','MeanCommonality','VarCommonality','SCCS','MeanSCCQ','VarSCCQ',... %26
                    'MeanAdvanceTF','MaxAdvanceTF','VarAdvanceTF',... %29
                    'MeanAdvanceDF','MaxAdvanceDF','VarAdvanceDF',... %32
                    'AvgQL',...
                    TWS{tw},...
                     'WordCount',... %35
                     'sum_CtiAdvDF','mean_CtiAdvDF','max_CtiAdvDF'... %38
                     'sum_idfRatio','mean_idfRatio','max_idfRatio','min_idfRatio',... %42
                     'mean_idfStem','max_idfStem','var_idfStem',... %45
                     'idfRatioMinMax','correlationTermsIdf','lstmstTermsDF','GammaStem'... %49
                     'sum_IdfAdvDF','mean_IdfAdvDF','max_IdfAdvDF','min_IdfAdvDF', ... %53
                     'Chi2DFTF','correlationTermsIctf','lstmstTermsTF',... %56
                     'sum_TFIdfNoStem','sum_TFIdfStem'
                     });
                 
                    SelectedFeatures.WordCount=double(SelectedFeatures.WordCount)+1;
                    
                    SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);
                    
                    
                    [p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.(STEMMERS{s}));
                    Scores=Joined(:,{'NoStem',STEMMERS{s}});

                    
                    
                    Y=[table2array(Scores) label];
                    X=table2array(SelectedFeatures);
                    
                    trainSetInx = Y(:,1) ~= Y(:,2) ; %Discard All same all zero 
                    X=X(trainSetInx,:);
                    Y=Y(trainSetInx,:);
                    
                    [numSample, featureSize]=size(X);


                    ff=zeros(1,featureSize);

                    opts = statset('display','iter','UseParallel',1);
                   
                       [inmodel,history] = sequentialfs(@criteriaFunCubicKNN,X,Y,'cv',numSample,'direction','forward','options',opts);
                      ff=ff+inmodel;
                    

                    fb=zeros(1,featureSize);
                    
                       [inmodel,history] = sequentialfs(@criteriaFunCubicKNN,X,Y,'cv',numSample,'direction','backward','options',opts);
                       fb=fb+inmodel;
                    

                    fileID = fopen(strcat('ff','criteriaFunCubicKNN',STEMMERS{s},dataNameR,'.txt'),'w');
                    fprintf(fileID,'%i\t',ff);
                    fclose(fileID);

                    fileID = fopen(strcat('fb','criteriaFunCubicKNN',STEMMERS{s},dataNameR,'.txt'),'w');
                    fprintf(fileID,'%i\t',fb);
                    fclose(fileID);
                end

                
            end
    end
end


function Joined = JoinTables(terms,features,runtopic,TF)
                   
    terms.IdfAdvDF=terms.idfs .* terms.advanceDF;
    terms.CtiAdvDF=terms.ctis .* terms.advanceDF;
    terms.idfRatio = terms.idfStem ./ terms.idfs;
    terms.ictfStem = log(TF./terms.TF); 
    terms.TFIdfNoStem = terms.TFNoStem .* terms.idfs;
    terms.TFIdfStem = terms.TF .* terms.idfStem;
    termAgg1 = groupsummary(terms,'QueryID',{'mean','max','sum','var','min'},'IdfAdvDF');
    termAgg2 = groupsummary(terms,'QueryID',{'mean','max','sum','var'},'CtiAdvDF');
    termAgg3 = groupsummary(terms,'QueryID',{'mean','max','sum','var','min'},'idfRatio');
    termAgg4 = groupsummary(terms,'QueryID',{'mean','max','var'},'idfStem');
    termAgg5 = groupsummary(terms,'QueryID',{'sum'},'TFIdfNoStem');
    termAgg6 = groupsummary(terms,'QueryID',{'sum'},'TFIdfStem');

    gamma1Table = Gamma1(terms(:,{'QueryID','word','idfStem'}));

    Joined = join(features,runtopic,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg1,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg2,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg3,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg4,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg5,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,termAgg6,'LeftKeys',1,'RightKeys',1);
    Joined.idfRatioMinMax = Joined.min_idfRatio ./ Joined.max_idfRatio;
    
    dftfTable = chi2(terms(:,{'QueryID','DFNoStem','TFNoStem','DF','TF'}));

    corrTable = IdfOrderDist(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
    corrTableICTF = ictfOrderDist(terms(:,{'QueryID','word','ictfs'}), terms(:,{'QueryID','word','ictfStem'}));
    lstmstTable = LSTMST(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
    lstmstTableTF = LSTMSTTF(terms(:,{'QueryID','word','ictfs'}), terms(:,{'QueryID','word','ictfStem'}));
    Joined = join(Joined,corrTable,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,lstmstTable,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,gamma1Table,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,dftfTable,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,corrTableICTF,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,lstmstTableTF,'LeftKeys',1,'RightKeys',1);

    
end



%min/max
function dftfTable = chi2(DFTFs)
    QIDs = unique(DFTFs.QueryID);
    
    dftfTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','Chi2DFTF'});
    
    [Nn,edgesn,binDN] = histcounts(DFTFs.DFNoStem,'BinMethod','fd');
    [Ns,edgess,binDS] = histcounts(DFTFs.DF,'BinMethod','fd');
    
    [Nn,edgesn,binTN] = histcounts(DFTFs.TFNoStem,'BinMethod','fd');
    [Ns,edgess,binTS] = histcounts(DFTFs.TF,'BinMethod','fd');
    
    DFTFs.BinDN = binDN;
    DFTFs.binDS = binDS;
    
    DFTFs.binTN = binTN;
    DFTFs.binTS = binTS;
    
    for i=1:size(QIDs,1)
        terms = DFTFs(DFTFs.QueryID==QIDs(i),:);
%         obs1=[terms.DFNoStem' terms.TFNoStem'];
%         obs2=[terms.DF' terms.TF'];
        obs1=[terms.BinDN' terms.binTN'];
        obs2=[terms.binDS' terms.binTS'];
        c= chi2Val([obs1;obs2]);
        p=1-chi2cdf(c,length(obs1)-1);
        dftfTable.Chi2DFTF(dftfTable.QueryID == QIDs(i))=p;
    end
        
end

% 
% function corrTable = IdfOrderDist(noStemTerms, stemTerms) 
%     QIDNoStem = unique(noStemTerms.QueryID);
%     QIDStem = unique(noStemTerms.QueryID);
%     
%     if ~isequal(QIDNoStem,QIDStem) 
%         error('QIDS not equal')
%     end
%     
%     QIDs = QIDNoStem;
%     clear QIDNoStem
%     clear QIDStem
%     corrTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','correlationTerms'});
%     
%     for i=1:size(QIDs,1)
%         termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
%         termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
%         
%         wordCount=size(termsNoStem,1);
%         order=1:wordCount;
%         termsNoStem.Position = order';
%         termsStem.Position = order';
%         
%        
%         termsNoStem = sortrows(termsNoStem, 'idfs'); % sort the table by 'DOB'
%         
%         termsStem = sortrows(termsStem, 'idfStem'); % sort the table by 'DOB'
%         
%         %if(lev(termsNoStem.Position, termsStem.Position)/wordCount > 0.1) 
%         c=corr(termsNoStem.Position,termsStem.Position,'Type','Spearman');
%   
%         corrTable.correlationTerms(corrTable.QueryID == QIDs(i))=c;
%       
%     end
%         
% end
% 
% function lstmstTable = LSTMST(noStemTerms, stemTerms)
%     QIDNoStem = unique(noStemTerms.QueryID);
%     QIDStem = unique(noStemTerms.QueryID);
%     
%     if ~isequal(QIDNoStem,QIDStem) 
%         error('QIDS not equal')
%     end
%     
%     QIDs = QIDNoStem;
%     clear QIDNoStem
%     clear QIDStem
%     lstmstTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','lstmstTerms'});
%     
%     for i=1:size(QIDs,1)
%         termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
%         termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
%         
%         wordCount = size(termsNoStem,1);
%         order=1:wordCount;
%         termsNoStem.Position = order';
%         termsStem.Position = order';
%         
%        
%         termsNoStem = sortrows(termsNoStem, 'idfs','ascend'); % sort the table by 'DOB'
%         
%         termsStem = sortrows(termsStem, 'idfStem','ascend'); % sort the table by 'DOB'
%         
%         c=0;
%         if(termsNoStem.Position(1) == termsStem.Position(1)) 
%             c=c+0.5;
%         end
%         if(termsNoStem.Position(wordCount) == termsStem.Position(wordCount)) 
%             c=c+0.5;
%         end
%         
%         lstmstTable.lstmstTerms(lstmstTable.QueryID == QIDs(i))=c;
%     end
%         
% end
% 

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


function JoinedT = extTrainData(coll,measure,tw,stemmer,TF)

    if strcmp(coll, 'CW09B')
        
        dataNameRT = strcat('MQ09','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ09',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);
    elseif strcmp(coll, 'MQ09')
        
        dataNameRT = strcat('CW09B','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','CW09B',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);
    elseif strcmp(coll, 'CW12B')
        
        dataNameRT = strcat('NTCIR','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','NTCIR',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);
    elseif strcmp(coll, 'NTCIR')
        
        dataNameRT = strcat('CW12B','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','CW12B',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);
    elseif strcmp(coll, 'GOV2')
        
        dataNameRT = strcat('MQ07','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ07',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);

        dataNameRT = strcat('MQ08','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ08',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT,TF);
        
        JoinedT=vertcat(JoinedT,JoinedT2);
    elseif strcmp(coll, 'MQ07')
        
        dataNameRT = strcat('GOV2','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','GOV2',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);

        dataNameRT = strcat('MQ08','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ08',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT,TF);
        JoinedT=vertcat(JoinedT,JoinedT2);
    elseif strcmp(coll,'MQ08')
        
        dataNameRT = strcat('GOV2','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','GOV2',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT = JoinTables(termsT,featuresT,runtopicT,TF);

        dataNameRT = strcat('MQ07','_',measure,'_',tw);
        runtopicT = evalin('base',dataNameRT);

        dataNameFT = strcat('Feature','MQ07',stemmer);
        featuresT = evalin('base',dataNameFT);

        dataNameFTermsT = strcat(dataNameFT,'Term');    
        termsT=evalin('base',dataNameFTermsT);

        JoinedT2 = JoinTables(termsT,featuresT,runtopicT,TF);
        JoinedT=vertcat(JoinedT,JoinedT2);
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
    
    lstmstTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','lstmstTermsDF'});
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount = size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'idfs','ascend'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'idfStem','ascend'); % sort the table by 'DOB'
        
        
        if(termsNoStem.Position(1) == termsStem.Position(1) || termsNoStem.Position(wordCount) == termsStem.Position(wordCount)) 
            predictedLabel=1;
        else
            predictedLabel=0;
        end
        
        lstmstTable.lstmstTermsDF(lstmstTable.QueryID == QIDs(i))=predictedLabel;
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
    corrTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','correlationTermsIdf'});
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
        if( c> 0.5 ) 
            predictedLabel=1;
        else
            predictedLabel=0;
        end
        
        corrTable.correlationTermsIdf(corrTable.QueryID == QIDs(i))=predictedLabel;
    end
        
end



function lstmstTable = LSTMSTTF(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    lstmstTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','lstmstTermsTF'});
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount = size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'ictfs','ascend'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'ictfStem','ascend'); % sort the table by 'DOB'
        
        if(termsNoStem.Position(1) == termsStem.Position(1) || termsNoStem.Position(wordCount) == termsStem.Position(wordCount)) 
            predictedLabel=1;
        else
            predictedLabel=0;
        end
        
        lstmstTable.lstmstTermsTF(lstmstTable.QueryID == QIDs(i))=predictedLabel;
    end
        
end



function corrTable = ictfOrderDist(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    corrTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','correlationTermsIctf'});
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'ictfs'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'ictfStem'); % sort the table by 'DOB'
        
        %if(lev(termsNoStem.Position, termsStem.Position)/wordCount > 0.1) 
        c=corr(termsNoStem.Position,termsStem.Position,'Type','Spearman');
        if( c> 0.5 ) 
            predictedLabel=1;
        else
            predictedLabel=0;
        end
        
        corrTable.correlationTermsIctf(corrTable.QueryID == QIDs(i))=predictedLabel;
    end
        
end



