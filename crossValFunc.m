function [ms, significant, m1, m2, m3,oracle ]  = crossValFunc(trainX2,trainY2,testX,testY)
%CRITERIAFUN Summary of this function goes here
%   Detailed explanation goes here

LabelsY=categorical(trainY2(:,4));
%Mdl = fitcknn(trainX2,categorical(trainY2(:,3)),'NumNeighbors',10,'Standardize',1,'Distance','minkowski');
%Mdl = fitcdiscr(trainX2,categorical(trainY2(:,1)));
Mdl = fitcsvm(...
    trainX2, ...
    LabelsY, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 2, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 10, ...
    'Standardize', true, ...
    'ClassNames', categorical({'0'; '1'; '2'}));

[label,score,cost] = predict(Mdl,testX);

[ms, significant, m1, m2, m3,oracle ] = AverageNDCG(testY,label);

end



