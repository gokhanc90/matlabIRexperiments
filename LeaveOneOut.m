%-----INIT----------------
% {'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF',...
% 'MeanCTI','VarCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
% 'SCCSNoStem','MeanSCCQNoStem','VarSCCQNoStem','MeanCommonalityKStem','VarCommonalityKStem',...
% 'SCCSKStem','MeanSCCQKStem','VarSCCQKStem','MeanCommonalitySnowball','VarCommonalitySnowball',...
% 'SCCSSnowball','MeanSCCQSnowball','VarSCCQSnowball',...
% 'Chi2DFTF','Chi2IdfIctf','Chi2SCQ','Chi2Commonalities','Chi2SCCQ',...
% 'MeanSCCS','HarmMeanSCCS','MeanVarCommonality','MeanVarSCCQ'}

 Filtered=MQ09TypeQ(MQ09TypeQ.AllSameAllZero == '0',:);
 SelectedFeatures=Filtered(:,{'MeanCTI','VarSCCQNoStem','VarSCCQKStem','VarSCCQNoStem', ...
 });
 SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);
Label=Filtered(:,4);
oracleFiltered=Filtered(:,5);
ScoresFiltered=Filtered(:,[6 7 8]);
%---------------------------------------
%runtopic=Filtered(:,5:8);

option=2; %0:combination 1:remove add else: all

fileID = fopen('runtopic.txt','a');
% functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunDiscriminateQuadratic ,@criteriaFunEnsembleRUSBoost ,...
% @criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
% 	@criteriaFunMediumKNN ,@criteriaFunSVM };

functions={@criteriaFunGaussianNaiveBayes };


Y=[table2array(ScoresFiltered) double(table2array(Label))-1];

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

            [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f FeatureSize: %f Feature: %s\n',func2str(functions{K}),...
                ms,significant,m1,m2,m3,oracle, S, string(strjoin(C(Ci,:))));
        end
       %    runtopic(:,S+4)=table(predictionScores);
       %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
        end
    end
elseif  option==1   
    for S =1:n
        X=table2array(SelectedFeatures(:,S));
        %X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
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

            [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f Allow: %s\n',func2str(functions{K}),...
                ms,significant,m1,m2,m3,oracle, SelectedFeatures.Properties.VariableNames{S});
        end
       %    runtopic(:,S+4)=table(predictionScores);
       %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
        
    end
else
    X=table2array(SelectedFeatures);
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

        [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
        fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f  %s\n',func2str(functions{K}),...
            ms,significant,m1,m2,m3,oracle, 'All');
        
   %     runtopic(:,end)=table(predictionScores);
    %    runtopic.Properties.VariableNames{end}='All';

    end
end


