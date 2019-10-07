function [criterion, ms, significant, m1, m2, m3,oracle,labels ] = criteriaFunEnsembleRUSBoost(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,4));

% %Ensemble RUSBoost%
template = templateTree(...
    'MaxNumSplits', 20);
Mdl = fitcensemble(...
    trainX2, ...
    LabelsY, ...
    'Method', 'RUSBoost', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'LearnRate', 0.1, ...
    'ClassNames', categorical({'0'; '1'; '2'}));
%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

