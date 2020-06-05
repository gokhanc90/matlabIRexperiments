function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunEnsembleRUSBoost(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

% %Ensemble RUSBoost%
template = templateTree(...
    'MaxNumSplits', size(LabelsY,1));
Mdl = fitcensemble(...
    trainX2, ...
    LabelsY, ...
    'Method', 'RUSBoost', ...
    'NumLearningCycles', 150, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', categorical({'0'; '1'}));
%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

