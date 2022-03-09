function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunFineTree(trainX,trainY,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY(:,3));
[Z,mu,sigma] = zscore(trainX);
testX = (testX-mu) ./ sigma;
trainX = normalize(trainX);
% %FineTree%
Mdl = fitctree(...
    trainX, ...
    LabelsY, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', size(LabelsY,1), ...
    'Surrogate', 'off', ...
    'Prune','off', ...
    'PruneCriterion','impurity', ...
    'ScoreTransform','logit',...
    'ClassNames', categorical({'0'; '1'}));
% 

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

