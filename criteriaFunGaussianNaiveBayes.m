function [criterion, ms, significant, m1, m2, oracle,labels ] = criteriaFunGaussianNaiveBayes(trainX,trainY,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here
LabelsY=categorical(trainY(:,3));


% 
%Gaussian Naive Bayes%
%distributionNames =  repmat({'normal'}, 1, size(trainX2,2));

Mdl = fitcnb(...
        trainX, ...
        LabelsY, ...
        'DistributionNames', {'normal','normal','normal','normal','normal','kernel','normal','normal','kernel'}, ...
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

