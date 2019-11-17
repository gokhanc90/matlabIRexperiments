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

Filtered=MQ07TypeQ(MQ07TypeQ.AllSameAllZero == '0',:);
 SelectedFeatures=Filtered(:,{ 
'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF',...
'MeanCTI','VarCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
'SCCSNoStem','MeanSCCQNoStem','VarSCCQNoStem',...
'MeanCommonalityKStem','VarCommonalityKStem','SCCSKStem','MeanSCCQKStem','VarSCCQKStem',...
'MeanCommonalitySnowball','VarCommonalitySnowball','SCCSSnowball','MeanSCCQSnowball','VarSCCQSnowball',...
'Chi2DFTF','Chi2IdfIctf','Chi2SCQ','Chi2Commonalities','Chi2SCCQ',...
'BM25CollNoStem','BM25CollKStem','BM25CollSnowball',...
'AdvanceKStem','AdvanceSnowball','BM25AdvKStem','BM25AdvSnowball'...
 });
 SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);

Label=Filtered(:,4);
oracleFiltered=Filtered(:,5);
ScoresFiltered=Filtered(:,[6 7 8]);

% act={'euclidean','chebychev','cosine','seuclidean','cityblock','minkowski','correlation','spearman','hamming','jaccard'};
% 
% l=table2array(Label);
% goodrows = not(any(isnan(table2array(SelectedFeatures)),2));
% 
% v = double(categorical(l));
% c = full(sparse(1:numel(v),v,ones(size(v)),numel(v),3));
% X=normalize(table2array(SelectedFeatures));
% for a=1:length(act)
%     Y = tsne(X,'Algorithm','exact','Distance',act{a},'Standardize',1,'NumDimensions',3,'LearnRate',500);
%     subplot(4,3,a)
%     scatter3(Y(:,1),Y(:,2),Y(:,3),5,c,'filled');
%     title(act{a})
% end

%----------------------RICA----------------------
functions={@criteriaFunGaussianNaiveBayes };


Y=[table2array(ScoresFiltered) double(table2array(Label))-1];

[m, n]=size(SelectedFeatures);

X=table2array(SelectedFeatures);
mdl = rica(X,10,'IterationLimit',1e6,'Standardize',1);
X = transform(mdl,X);
%X=normalize(X);
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

    [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(Y(:,[1 2 3]),predictedlabel);
    fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f KStemMean: %f SnowballMean: %f Oracle: %f  %s\n',func2str(functions{K}),...
        ms,significant,m1,m2,m3,oracle, 'All');
end
    

