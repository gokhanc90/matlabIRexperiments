function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunDiscriminateQuadratic(trainX,trainY,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY(:,3));

% %Discriminate Quadratic%
Mdl = fitcdiscr(...
    trainX, ...
    LabelsY, ...
    'DiscrimType', 'pseudolinear', ...
    'FillCoeffs', 'off', ...
    'ScoreTransform','logit', ...
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

