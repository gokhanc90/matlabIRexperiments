%-----INIT----------------
%{'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF',...
%'MeanCTI','VarCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
%'SCCSNoStem','MeanSCCQNoStem','VarSCCQNoStem',...
%'MeanCommonalityKStem','VarCommonalityKStem','SCCSKStem','MeanSCCQKStem','VarSCCQKStem',...
%'MeanCommonalitySnowball','VarCommonalitySnowball','SCCSSnowball','MeanSCCQSnowball','VarSCCQSnowball',...
%'Chi2DFTF','Chi2IdfIctf','Chi2SCQ','Chi2Commonalities','Chi2SCCQ',...
%'MeanSCCS','HarmMeanSCCS','MeanVarCommonality','MeanVarSCCQ',...
%'BM25CollNoStem','BM25CollKStem','BM25CollSnowball'...
%}


 %Filtered=MQ07TypeQ(MQ07TypeQ.AllSameAllZero == '0',:);
 SelectedFeatures=Filtered(:,{ 
'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF',...
'MeanCTI','VarCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
'SCCSNoStem','MeanSCCQNoStem','VarSCCQNoStem',...
'Chi2DFTF','Chi2IdfIctf','Chi2SCQ',...
'BM25CollNoStem','BM25CollKStem','BM25CollSnowball',...
'AdvanceKStem','AdvanceSnowball','BM25AdvKStem','BM25AdvSnowball'...
 });
 SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);

Label=Filtered(:,4);
oracleFiltered=Filtered(:,5);
ScoresFiltered=Filtered(:,[6 7 8]);

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

fileID = fopen('runtopic.txt','a');

Y=double(table2array(Label))-1;
[m, n]=size(SelectedFeatures);

%Activation functions%
act={'sig','sin','hardlim','tribas','radbas'};
neurons={3, 10, 20, 50, 75, 100, 150};
C={2^(-10), 2^(-5), 2^(5), 2^(10), 2^(20), 2^(35), 2^(50)};


if option==0
    for S =[2:4,34:n-1]
        C=nchoosek(SelectedFeatures.Properties.VariableNames,uint16(S));
        for Ci = 1:length(C)
        %X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
        X=table2array(SelectedFeatures(:,C(Ci,:)));
        X=normalize(X);
       
            predictionScores=zeros(m,1);
            predictedlabel=categorical(zeros(m,1));
            for i=1:m
                Xtest=X(i,:);
                Xtrain=X([1:i-1,i+1:end],:);

                Ytest=Y(i,:);
                Ytrain=Y([1:i-1,i+1:end],:);

                
                %[TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Predictions] = ELM(Xtrain,Ytrain, Xtest,Ytest, 1, 50, 'tribas');
                [TrainingTime, TestingTime, test_sensitivity, test_specificity, test_gmean,Predictions] = ELM_regularized_NXN(Xtrain,Ytrain, Xtest,Ytest, 1, 3, 'hardlim',2^(30));
                %predictionScores(i)=ms; 
                predictedlabel(i)=Predictions;
            end

            [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f FeatureSize: %f Feature: %s\n',ELM_regularized_NXN,...
                ms,significant,m1,m2,m3,oracle, S, string(strjoin(C(Ci,:))));
       
       %    runtopic(:,S+4)=table(predictionScores);
       %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
        end
    end
elseif  option==1   
    for S =1:n
        %X=table2array(SelectedFeatures(:,S));
        X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
        X=normalize(X);
       
            predictionScores=zeros(m,1);
            predictedlabel=categorical(zeros(m,1));
            for i=1:m
                Xtest=X(i,:);
                Xtrain=X([1:i-1,i+1:end],:);

                Ytest=Y(i,:);
                Ytrain=Y([1:i-1,i+1:end],:);

                %[TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Predictions] = ELM(Xtrain,Ytrain, Xtest,Ytest, 1, 50, 'tribas');
                [TrainingTime, TestingTime, test_sensitivity, test_specificity, test_gmean,Predictions] = ELM_regularized_LXL(Xtrain,Ytrain, Xtest,Ytest, 1, 100, 'hardlim',2^(30));
                %predictionScores(i)=ms; 
                predictedlabel(i)=Predictions;
            end

            [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f Discard: %s\n','ELM_regularized_LXL',...
                ms,significant,m1,m2,m3,oracle, SelectedFeatures.Properties.VariableNames{S});
       
       %    runtopic(:,S+4)=table(predictionScores);
       %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
        
    end
else
    X=table2array(SelectedFeatures);
    X=normalize(X);
        for nr=1:length(neurons)
            for ac=1:length(act)
                for c1=1:length(C)
                    predictionScores=zeros(m,1);
                    predictedlabel=categorical(zeros(m,1));
                    for i=1:m
                        Xtest=X(i,:);
                        Xtrain=X([1:i-1,i+1:end],:);

                        Ytest=Y(i,:);
                        Ytrain=Y([1:i-1,i+1:end],:);

                        %[TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Predictions] = ELM(Xtrain,Ytrain, Xtest,Ytest, 1, 50, 'tribas');
                        [TrainingTime, TestingTime, test_sensitivity, test_specificity, test_gmean,Predictions] = ELM_regularized_LXL(Xtrain,Ytrain, Xtest,Ytest, 1, neurons{nr}, act{ac},C{c1});
                        %predictionScores(i)=ms; 
                        predictedlabel(i)=Predictions;

                    end
                    [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
                    fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f Neuron: %f ActF: %s C: %f %s\n','ELM_regularized_NXN',...
                        ms,significant,m1,m2,m3,oracle, neurons{nr}, act{ac},C{c1},'All');
                end
            end
        end
   %     runtopic(:,end)=table(predictionScores);
    %    runtopic.Properties.VariableNames{end}='All';

end


