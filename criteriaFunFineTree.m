function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunFineTree(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

% %FineTree%
Mdl = fitctree(...
    trainX2, ...
    LabelsY, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', size(LabelsY,1), ...
    'Surrogate', 'all', ...
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

