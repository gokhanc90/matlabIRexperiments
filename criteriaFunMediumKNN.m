function [criterion, ms, significant, m1, m2,m3, oracle,labels ] = criteriaFunMediumKNN(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY2(:,4));

% %Medium KNN%    
Mdl = fitcknn(...
    trainX2, ...
    LabelsY, ...
    'Distance', 'Euclidean', ...
    'Exponent', [], ...
    'NumNeighbors', 10, ...
    'DistanceWeight', 'Equal', ...
    'Standardize', true, ...
    'ClassNames', categorical({'0'; '1'; '2'}));    
    

%-------------------------------------------------%

[pred,ci] = predict(Mdl,testX);
if ~iscategorical(pred) 
    error('Predictions must be categorical')
end
labels=pred;
 [ms, significant, m1, m2, m3,oracle ] = AverageNDCG(testY,pred);


criterion=1-ms;
end

