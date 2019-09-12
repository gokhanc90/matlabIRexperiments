function [expVals] = expectedVal(obsVectors)
%Calculates expected values for chi2 
%Each row points to one observation

rowSums=sum(obsVectors,2);
colSums=sum(obsVectors,1);
totalSum=sum(obsVectors,'all');

expVals=rowSums.*colSums;
expVals=expVals./totalSum;
end

