
% functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunDiscriminateQuadratic ,@criteriaFunEnsembleRUSBoost ,...
% 	@criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
% 	@criteriaFunMediumKNN ,@criteriaFunSVM };

functions={@criteriaFunGaussianNaiveBayes };

% Features=[1 4 10 45 18 35 16 44 39 38 37 8];
% Y=[table2array(ScoresFiltered) double(table2array(Label))-1];
% X=table2array(SelectedFeatures(:,Features));

Y=[table2array(ScoresFiltered) double(table2array(Label))-1];

[m, n]=size(SelectedFeatures);

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

            [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
            predictionScores(i)=ms; 
            predictedlabel(i)=labels;
        end
        %mean(predictionScores)

        %[h,p] = ttest(predictionScores,table2array(ScoresFiltered(:,2)),'Alpha',0.05);
        [ms, significant, m1, m2, oracle ] = AverageNDCG(table2array(ScoresFiltered),predictedlabel);
        fprintf('MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f\n',func2str(functions{K}),...
            ms,significant,m1,m2,oracle);
        
        runtopic(:,S+3)=table(predictionScores);
        runtopic.Properties.VariableNames{S+3}=strcat('Minus',num2str(S));
    end
end


