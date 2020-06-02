load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/Features.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CW09B.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CW12B.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/GOV2.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ07.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ08.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/MQ09.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/NTCIR.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/WSJ.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/FeaturesTerm.mat');
load('/home/ubuntu/Documents/MATLAB/matlabIRexperiments/CollectionStats.mat');



STEMMERS={'KStem'  };
%STEMMERS={'KStem'};
%TWS={'BM25' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'LGD' 'PL2'};
TWS={'LGD'};
MEASURES={'NDCG20'};
%MEASURES={'NDCG20'};
COLLECTIONS={ 'CW09B' 'CW12B'  'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'};


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
                    terms.ictfStem = log(TF./terms.TF);
                    Joined = JoinTables(terms,features,runtopic);
                    
                    [p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.(STEMMERS{s}));
                    Scores=Joined(:,{'NoStem',STEMMERS{s}});
                    Y=[table2array(Scores) label];
                    X=label;    
                    trainSetInx = Y(:,1) ~= Y(:,2);
                    %%%%%%%%%%%%%%%%%%
%                     diff=abs(Joined.NoStem-Joined.(STEMMERS{s}));
%                     trainSetInx = diff > std(diff) / sqrt(size(diff,1));
% 
%                     a=Joined(trainSetInx,:);
%                     l=label(trainSetInx);
%                     idx=ismember(terms.QueryID,a.QueryID);
%                    terms=terms(idx,:);
                    fileID = fopen('Heuristic.txt','a');
               
                
%                     predictedLabels = IdfOrderChange(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','IdfOrderChange',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
% 
%                     predictedLabels = DFOrderChange(terms(:,{'QueryID','word','DFNoStem'}), terms(:,{'QueryID','word','DF'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','DFOrderChange',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
%                     
%                     bins=[50 100 150 250 500 750 1000 1250 1500 2000 5000 7500 10000];
%                     for bin=bins
%                         predictedLabels = DFOrderBinning(terms(:,{'QueryID','word','DFNoStem'}), terms(:,{'QueryID','word','DF'}),docCount,bin);
%                         [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                         fprintf(fileID,'Func: %s %d Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','DFOrderBinning',bin,...
%                             ms,significant,m1,m2,oracle,dataNameR,dataNameF);
%                     end
%                     
%                     predictedLabels = DFOrderBinningScott(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','DFOrderBinningScott',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
%                     X=[X, predictedLabels];
%                     
                     predictedLabels = IdfOrderDist(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
                    [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
                    fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','IdfOrderDist',...
                        ms,significant,m1,m2,oracle,dataNameR,dataNameF);
                    X=[X, predictedLabels];
                    
                     predictedLabels = LSTMST(terms(:,{'QueryID','word','idfs'}), terms(:,{'QueryID','word','idfStem'}));
                    [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
                    fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','LSTMST',...
                        ms,significant,m1,m2,oracle,dataNameR,dataNameF);
                    X=[X, predictedLabels];
                   
%                    
                   
%                     predictedLabels = IctfOrderChange(terms(:,{'QueryID','word','ictfs'}), terms(:,{'QueryID','word','ictfStem'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','IctfOrderChange',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
% 
%                     predictedLabels = TFOrderChange(terms(:,{'QueryID','word','TFNoStem'}), terms(:,{'QueryID','word','TF'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','TFOrderChange',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
                    
%                     bins=[50 100 150 250 500 750 1000 1250 1500 2000 5000 7500 10000];
%                     for bin=bins
%                         predictedLabels = TFOrderBinning(terms(:,{'QueryID','word','TFNoStem'}), terms(:,{'QueryID','word','TF'}),TF,bin);
%                         [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                         fprintf(fileID,'Func: %s %d Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','TFOrderBinning',bin,...
%                             ms,significant,m1,m2,oracle,dataNameR,dataNameF);
%                     end
%                     
%                     predictedLabels = TFOrderBinningScott(terms(:,{'QueryID','word','TFNoStem'}), terms(:,{'QueryID','word','TF'}));
%                     [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
%                     fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','TFOrderBinningScott',...
%                         ms,significant,m1,m2,oracle,dataNameR,dataNameF);
%                     X=[X, predictedLabels];
%                     
                     predictedLabels = ictfOrderDist(terms(:,{'QueryID','word','ictfs'}), terms(:,{'QueryID','word','ictfStem'}));
                    [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
                    fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','ictfOrderDist',...
                        ms,significant,m1,m2,oracle,dataNameR,dataNameF);
                    X=[X, predictedLabels];
                    
                     predictedLabels = LSTMSTTF(terms(:,{'QueryID','word','ictfs'}), terms(:,{'QueryID','word','ictfStem'}));
                    [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedLabels);
                    fprintf(fileID,'Func: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s\n','LSTMSTTF',...
                        ms,significant,m1,m2,oracle,dataNameR,dataNameF);
                    X=[X, predictedLabels];
                    X2= X(trainSetInx,:);
                    Y=X2(:,1);
                    X2=X2(:,2:end);
                    fclose(fileID);
            end
        end
   end
end
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

    Joined = join(Joined,gamma1Table,'LeftKeys',1,'RightKeys',1);
    Joined = join(Joined,dftfTable,'LeftKeys',1,'RightKeys',1);

    
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

function predictedLabels = LSTMST(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
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
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = IdfOrderChange(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        order=1:size(termsNoStem,1);
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'idfs'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'idfStem'); % sort the table by 'DOB'
        
        if(termsNoStem.Position == termsStem.Position) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end



function predictedLabels = DFOrderChange(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        order=1:size(termsNoStem,1);
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'DFNoStem'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'DF'); % sort the table by 'DOB'
        
        if(termsNoStem.Position == termsStem.Position) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = IdfOrderDist(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
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
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = DFOrderBinningScott(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    
    [Nn,edgesn,binn] = histcounts(noStemTerms.DFNoStem,'BinMethod','fd');
    [Ns,edgess,bins] = histcounts(stemTerms.DF,'BinMethod','fd');
    
    noStemTerms.Bin = binn;
    stemTerms.Bin = bins;
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
        dif = abs(termsNoStem.Bin-termsStem.Bin);
        
        
        termsNoStem = sortrows(termsNoStem, 'Bin'); 
        
        termsStem = sortrows(termsStem, 'Bin'); 
        
        c=compareTieAware(termsNoStem, termsStem);
        
        if(c) %if(sum(dif)>wordCount) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = DFOrderBinning(noStemTerms, stemTerms, docCount, numberOfBins) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
        binNoStem=zeros(wordCount,1);
        binStem=zeros(wordCount,1);
        for k=1:wordCount
            dfNoStem = termsNoStem.DFNoStem(k);
            dfStem = termsStem.DF(k);
            binNoStem(k) = floor(dfNoStem*numberOfBins/docCount);
            binStem(k) = floor(dfStem*numberOfBins/docCount);
        end
 
        termsNoStem.Bin = binNoStem;
        termsStem.Bin = binStem;
        
        termsNoStem = sortrows(termsNoStem, 'Bin'); 
        
        termsStem = sortrows(termsStem, 'Bin'); 
        
        c=compareTieAware(termsNoStem, termsStem);
        
        if(c) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end

function c=compareTieAware(termsNoStem, termsStem)
    wordCount = size(termsNoStem,1);
    c=0;
    for i=1:wordCount
            p1=i;
            p2=i;
            s1=1;
            s2=1;
            while 1
                if termsNoStem.Position(p1) == termsStem.Position(p2)
                    p1=p1+1;
                    p2=p2+1;
                    break;
                else 
                    if(p1==wordCount || p2==wordCount ) 
                        c=0;
                        return;
                    end
                    if termsNoStem.Bin(p1) == termsNoStem.Bin(p1 + s1) 
                        temp = termsNoStem(p1,:);
                        termsNoStem(p1,:) = termsNoStem(p1+s1,:);
                        termsNoStem(p1+s1,:) = temp;
                        s1 = s1+1;
                    elseif termsStem.Bin(p2) == termsStem.Bin(p2 + s2) 
                        temp = termsStem(p2,:);
                        termsStem(p2,:) = termsStem(p2 + s2,:);
                        termsStem(p2 + s2,:) = temp;
                        s2 = s2+1;
                    
                    else
                        c=0;
                        return;
                    end
                end
            end
            
            if i==wordCount 
                c=1;
                return;
            end
    end
    c=0;
    return;

end

%  TF   %%%%%%%%%%%%%

function predictedLabels = LSTMSTTF(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
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
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = IctfOrderChange(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        order=1:size(termsNoStem,1);
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'ictfs'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'ictfStem'); % sort the table by 'DOB'
        
        if(termsNoStem.Position == termsStem.Position) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end



function predictedLabels = TFOrderChange(noStemTerms, stemTerms)
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        order=1:size(termsNoStem,1);
        termsNoStem.Position = order';
        termsStem.Position = order';
        
       
        termsNoStem = sortrows(termsNoStem, 'TFNoStem'); % sort the table by 'DOB'
        
        termsStem = sortrows(termsStem, 'TF'); % sort the table by 'DOB'
        
        if(termsNoStem.Position == termsStem.Position) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = ictfOrderDist(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
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
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = TFOrderBinningScott(noStemTerms, stemTerms) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    
    [Nn,edgesn,binn] = histcounts(noStemTerms.TFNoStem,'BinMethod','fd');
    [Ns,edgess,bins] = histcounts(stemTerms.TF,'BinMethod','fd');
    
    noStemTerms.Bin = binn;
    stemTerms.Bin = bins;
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
        dif = abs(termsNoStem.Bin-termsStem.Bin);
        
        
        termsNoStem = sortrows(termsNoStem, 'Bin'); 
        
        termsStem = sortrows(termsStem, 'Bin'); 
        
        c=compareTieAware(termsNoStem, termsStem);
        
        if(c) %if(sum(dif)>wordCount) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end


function predictedLabels = TFOrderBinning(noStemTerms, stemTerms, TF, numberOfBins) 
    QIDNoStem = unique(noStemTerms.QueryID);
    QIDStem = unique(noStemTerms.QueryID);
    if ~isequal(QIDNoStem,QIDStem) 
        error('QIDS not equal')
    end
    
    QIDs = QIDNoStem;
    clear QIDNoStem
    clear QIDStem
    predictedLabels=zeros(size(QIDs,1),1);
    
    for i=1:size(QIDs,1)
        termsNoStem = noStemTerms(noStemTerms.QueryID==QIDs(i),:);
        termsStem = stemTerms(stemTerms.QueryID==QIDs(i),:);
        
        wordCount=size(termsNoStem,1);
        order=1:wordCount;
        termsNoStem.Position = order';
        termsStem.Position = order';
        
        binNoStem=zeros(wordCount,1);
        binStem=zeros(wordCount,1);
        for k=1:wordCount
            tfNoStem = termsNoStem.TFNoStem(k);
            tfStem = termsStem.TF(k);
            binNoStem(k) = floor(tfNoStem*numberOfBins/TF);
            binStem(k) = floor(tfStem*numberOfBins/TF);
        end
 
        termsNoStem.Bin = binNoStem;
        termsStem.Bin = binStem;
        
        termsNoStem = sortrows(termsNoStem, 'Bin'); 
        
        termsStem = sortrows(termsStem, 'Bin'); 
        
        c=compareTieAware(termsNoStem, termsStem);
        
        if(c) 
            predictedLabels(i)=1;
        else
            predictedLabels(i)=0;
        end
    end
        
end
