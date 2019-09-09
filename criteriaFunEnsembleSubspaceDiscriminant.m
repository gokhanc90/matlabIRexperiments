function criterion = criteriaFunEnsembleSubspaceDiscriminant(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

%  %Ensemble Subspace Discriminant%
subspaceDimension = max(1, min(9, size(trainX2,2) - 1));
Mdl = fitcensemble(...
    trainX2, ...
    LabelsY, ...
    'Method', 'Subspace', ...
    'NumLearningCycles', 30, ...
    'Learners', 'discriminant', ...
    'NPredToSample', subspaceDimension, ...
    'ClassNames', categorical({'0'; '1'}));
 

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

