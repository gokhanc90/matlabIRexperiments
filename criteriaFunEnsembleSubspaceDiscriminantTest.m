function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunEnsembleSubspaceDiscriminantTest(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

%t = templateDiscriminant(...
%    'Delta',0,...
%    'DiscrimType','pseudoquadratic'... %'linear' (default) | 'quadratic' | 'diaglinear' | 'diagquadratic' | 'pseudolinear' | 'pseudoquadratic'
%    ); 

t = templateKNN('NumNeighbors',11,'Standardize',1,'Distance','Minkowski','Exponent',3);

%  %Ensemble Subspace Discriminant%
subspaceDimension = max(1, min(9, size(trainX2,2) - 1));
Mdl = fitcensemble(...
    trainX2, ...
    LabelsY, ...
    'Method', 'Subspace', ...
    'NumLearningCycles', 5, ...
    'Learners', t, ...
    'NPredToSample', subspaceDimension, ...
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

