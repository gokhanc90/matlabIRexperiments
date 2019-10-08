%-----INIT----------------
 %Filtered=MQ09TypeQ(MQ09TypeQ.AllSameAllZero == '0',:);
 SelectedFeatures=Filtered(:,[10:end]);
 SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);
Label=Filtered(:,4);
oracleFiltered=Filtered(:,5);
ScoresFiltered=Filtered(:,[6 7 8]);
%---------------------------------------
%runtopic=Filtered(:,5:8);
AllFeatures=0;
fileID = fopen('runtopic.txt','a');
 %functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunDiscriminateQuadratic ,@criteriaFunEnsembleRUSBoost ,...
 %	@criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
 %	@criteriaFunMediumKNN ,@criteriaFunSVM };

functions={@criteriaFunGaussianNaiveBayes };

% Features=[1 4 10 45 18 35 16 44 39 38 37 8];
% Y=[table2array(ScoresFiltered) double(table2array(Label))-1];
% X=table2array(SelectedFeatures(:,Features));

Y=[table2array(ScoresFiltered) double(table2array(Label))-1];

[m, n]=size(SelectedFeatures);

if AllFeatures==1
    for S =1:n
        X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
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
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f Discard: %s\n',func2str(functions{K}),...
                ms,significant,m1,m2,m3,oracle, SelectedFeatures.Properties.VariableNames{S});

  %          runtopic(:,S+4)=table(predictionScores);
   %         runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
        end
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
        fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f Discard: %s\n',func2str(functions{K}),...
            ms,significant,m1,m2,m3,oracle, 'All');
        
  %      runtopic(:,end)=table(predictionScores);
  %      runtopic.Properties.VariableNames{end}='All';

    end
end


