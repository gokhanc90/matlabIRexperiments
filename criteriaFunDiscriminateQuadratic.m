function criterion = criteriaFunDiscriminateQuadratic(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

% %Discriminate Quadratic%
Mdl = fitcdiscr(...
    trainX2, ...
    LabelsY, ...
    'DiscrimType', 'quadratic', ...
    'FillCoeffs', 'off', ...
    'ClassNames', categorical({'0'; '1'}));

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

