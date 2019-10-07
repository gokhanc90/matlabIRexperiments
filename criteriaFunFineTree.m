function [criterion, ms, significant, m1, m2, m3,oracle,labels ] = criteriaFunFineTree(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,4));

% %FineTree%
Mdl = fitctree(...
    trainX2, ...
    LabelsY, ...
    'SplitCriterion', 'gdi', ...
    'MaxNumSplits', 30, ...
    'Surrogate', 'off', ...
    'ClassNames', categorical({'0'; '1'; '2'}));
% 

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

