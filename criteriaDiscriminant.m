function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaDiscriminant(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

%trainX2 = normalize(trainX2);
%    'Delta',0,...
%    'DiscrimType','pseudoquadratic'... %'linear' (default) | 'quadratic' | 'diaglinear' | 'diagquadratic' | 'pseudolinear' | 'pseudoquadratic'

Mdl = fitcdiscr(...
    trainX2, ...
    LabelsY, ...
    'FillCoeffs', 'off', ...
    'DiscrimType','pseudolinear',...
    'ClassNames', categorical({'0'; '1'})...
); 

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

