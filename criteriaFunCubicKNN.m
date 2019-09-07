function criterion = criteriaFunCubicKNN(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,3));

% %Cubic KNN%
Mdl = fitcknn(...
    trainX2, ...
    LabelsY, ...
    'Distance', 'Minkowski', ...
    'Exponent', 3, ...
    'NumNeighbors', 10, ...
    'DistanceWeight', 'Equal', ...
    'Standardize', true, ...
    'ClassNames', categorical({'0'; '1'}));

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
 [ms, significant, m1, m2, oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

