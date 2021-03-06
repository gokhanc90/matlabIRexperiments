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


STEMMERS={'KStem'};
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
                
                    %Sf=[1   2   3   5   6   7   8   9  15   16 17  18  20  22  25  26   30  33  36]; %CW12B
                    %Sf=[3   4   9  10  11  12  15  16  19  22  24  27  28  29  32  33  34  37]; %GOV2
                    %Sf=[4   6   7   9  10  11  13  14  20  21  23  24  25  28  29  34  36  37]; %MQ08
                    %Sf=[1   4   5  11  13  15  17  19  20  21  22  23  24  27  30  31  32  33]; %MQ09
                    %Sf=[2   7   9  12  16  20  22  23  25  27  28  29  30  34  36  37]; %NTCIR
                    %Sf=[1   4   6   7  10  11  12  19  20  22  29  33  34  35  36]; %CW09B
                    %Sf=[5   6   8   9  10  11  13  18  20  26  27  36  37 38]; %WSJ
                    %Sf=[ 1 5  7   9  10  11  13  19  20  24  26  27  28  30  37  38  39  42  43];
                    %Sf=[ 1 5 6 7 8  9  10  11 12  19  20  24  26  27  28 30 31 33 36 37  38  39 42  43 45 47];
                    
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
%                             if strcmp(COLLECTIONS{coll},'WSJ')
%                                 JoinedT=Joined;
%                                 JoinedT(:,:)=[];
%                             else
%                                 JoinedT = extTrainData(COLLECTIONS{coll},MEASURES{measure},TWS{tw},STEMMERS{s});
%                             end
                    Sf=[ 1 5 6 7 8  9  10  11 12  19  20  24  26  27  28 30 31 33 36 37  38  39 42  43 45 47 ];
                    
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
                     'idfRatioMinMax','correlationTerms','lstmstTerms','GammaStem'... %49
                     'sum_IdfAdvDF','mean_IdfAdvDF','max_IdfAdvDF','min_IdfAdvDF', ... %53
                     'Chi2DFTF'
                     });
%                             SelectedFeaturesT = JoinedT(:,SelectedFeatures.Properties.VariableNames);
                    
%                     SelectedFeaturesT.WordCount=double(SelectedFeaturesT.WordCount)+1;
                    SelectedFeatures.WordCount=double(SelectedFeatures.WordCount)+1;
                    
                    %SelectedFeaturesT=fillmissing(SelectedFeaturesT,'constant',0);
                    SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);
                    
%                            SelectedFeaturesT=SelectedFeaturesT(:,Sf);
                    SelectedFeatures=SelectedFeatures(:,Sf);
                    
                    
                    [p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.(STEMMERS{s}));
                    Scores=Joined(:,{'NoStem',STEMMERS{s}});
                    
 %                         [pT,isSigT,oracleT,labelT]=getOracle(JoinedT.NoStem,JoinedT.(STEMMERS{s}));
 %                         ScoresT=JoinedT(:,{'NoStem',STEMMERS{s}});

                    % 
                    %  [idx,weights] = relieff(table2array(SelectedFeatures),table2array(Label),20);
                    %   poz=weights>0
                    %   pf=idx(poz)
                    %   vn=SelectedFeatures.Properties.VariableNames(pf)
                    %   SelectedFeatures = SelectedFeatures(:,vn);
                    %  
                    %---------------------------------------
                    %runtopic=Filtered(:,5:8);

                    option=2; %0:combination 1:remove add else: all

                    fileID = fopen('runtopic2.txt','a');
                    % functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunEnsembleRUSBoost ,...
                    % @criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
                    % 	@criteriaFunMediumKNN ,@criteriaFunSVM }; % , @criteriaFunCubicKNN, @criteriaFunEnsembleRUSBoost ,@criteriaFunDiscriminateQuadratic 

                    functions={@criteriaFunCubicKNN};


                    Y=[table2array(Scores) label];
%                            YT=[table2array(ScoresT) labelT];

                    [m, n]=size(SelectedFeatures);

                    if option==0
                        for S =[2:4,34:n-1]
                            C=nchoosek(SelectedFeatures.Properties.VariableNames,uint16(S));
                            for Ci = 1:length(C)
                            %X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
                            X=table2array(SelectedFeatures(:,C(Ci,:)));
                            for K = 1 : length(functions)

                                predictionScores=zeros(m,1);
                                predictedlabel=categorical(zeros(m,1));
                                for i=1:m
                                    Xtest=X(i,:);
                                    Xtrain=X([1:i-1,i+1:end],:);

                                    Ytest=Y(i,:);
                                    Ytrain=Y([1:i-1,i+1:end],:);

                                    [criterion, ms, significant, m1, m2,m3, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                    predictionScores(i)=ms; 
                                    predictedlabel(i)=labels;
                                end

                                [ms, significant, m1, m2, oracle ] = AverageNDCG(table2array(Scores),predictedlabel);
                                fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f FeatureSize: %f Feature: %s\n',func2str(functions{K}),...
                                    ms,significant,m1,m2,oracle, S, string(strjoin(C(Ci,:))));
                            end
                           %    runtopic(:,S+4)=table(predictionScores);
                           %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
                            end
                        end
                    elseif  option==1   
                        for S =1:n
                            %X=table2array(SelectedFeatures(:,S));
                            X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
                            for K = 1 : length(functions)

                                predictionScores=zeros(m,1);
                                predictedlabel=categorical(zeros(m,1));
                                for i=1:m
                                    Xtest=X(i,:);
                                    Xtrain=X([1:i-1,i+1:end],:);

                                    Ytest=Y(i,:);
                                    Ytrain=Y([1:i-1,i+1:end],:);

                                     trainSetInx = Ytrain(:,1) ~= Ytrain(:,2) ; %Discard All same all zero 
                                  Xtrain=Xtrain(trainSetInx,:);
                                  Ytrain=Ytrain(trainSetInx,:);

                                    
                                    [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                    predictionScores(i)=ms; 
                                    predictedlabel(i)=labels;
                                end

                                [ms, significant, m1, m2, oracle ] = AverageNDCG(table2array(Scores),predictedlabel);
                                fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f Discard: %s %s %s\n',func2str(functions{K}),...
                                    ms,significant,m1,m2,oracle, SelectedFeatures.Properties.VariableNames{S},dataNameR,dataNameF);
                            end
                           %    runtopic(:,S+4)=table(predictionScores);
                           %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};

                        end
                    else
                        X=table2array(SelectedFeatures);
%                                 XT=table2array(SelectedFeaturesT);
%                         mdl = rica(X,6,'IterationLimit',1e3,'Standardize',1);
%                         X = transform(mdl,X);
%                         X=normalize(X);
                        for K = 1 : length(functions)

                            predictionScores=zeros(m,1);
                            predictedlabel=categorical(zeros(m,1));
                            for i=1:m
                                Xtest=X(i,:);
                                Xtrain=X([1:i-1,i+1:end],:);
%                                         Xtrain=[Xtrain;XT];
                                
                                Ytest=Y(i,:);
                                Ytrain=Y([1:i-1,i+1:end],:);
%                                         Ytrain=[Ytrain;YT];
                                
%                                 diff=abs(Ytrain(:,1)-Ytrain(:,2));
%                                  trainSetInx = diff > std(diff) / sqrt(size(diff,1));
% %                                 
                                trainSetInx = Ytrain(:,1) ~= Ytrain(:,2) ; %Discard All same all zero 
                                  Xtrain=Xtrain(trainSetInx,:);
                                  Ytrain=Ytrain(trainSetInx,:);

                                [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                predictionScores(i)=ms; 
                                predictedlabel(i)=labels;
                            end

                            [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedlabel);
                            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s %s\n',func2str(functions{K}),...
                                ms,significant,m1,m2,oracle,dataNameR,dataNameF,strjoin(SelectedFeatures.Properties.VariableNames));

                       %     runtopic(:,end)=table(predictionScores);
                        %    runtopic.Properties.VariableNames{end}='All';

                        end
                    end
                
                
                
                
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
% 
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
%     
%     lstmstTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','lstmstTermsDF'});
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
%         
%         if(termsNoStem.Position(1) == termsStem.Position(1) || termsNoStem.Position(wordCount) == termsStem.Position(wordCount)) 
%             predictedLabel=1;
%         else
%             predictedLabel=0;
%         end
%         
%         lstmstTable.lstmstTermsDF(lstmstTable.QueryID == QIDs(i))=predictedLabel;
%     end
%         
% end
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
%     corrTable = array2table([QIDs zeros(size(QIDs))],'VariableNames',{'QueryID','correlationTermsIdf'});
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
%         if( c> 0.5 ) 
%             predictedLabel=1;
%         else
%             predictedLabel=0;
%         end
%         
%         corrTable.correlationTermsIdf(corrTable.QueryID == QIDs(i))=predictedLabel;
%     end
%         
% end
% 


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


