function [criterion, ms, significant, m1, m2,oracle, labels  ] = knnIR(trainX,trainY,testX,testY)
%KNNIR Summary of this function goes here
%   Detailed explanation goes here
k=5;
[Z,mu,sigma] = zscore(trainX);
testX = (testX-mu) ./ sigma;
trainX = normalize(trainX);
[m,n]=size(testX);
preds=categorical(zeros(m,1));
for i=1:m
 %   sim = similarity(trainX,testX(i,:));
   
 %   simLabel=[sim trainY];
    
    [CIDX,dist] = knnsearch(trainX, testX(i,:),'k',k,...
                       'Distance','minkowski','p',3);
    
 %   sorted = sortrows(simLabel,1,'ascend');
 %   neighbors = sorted(1:k,:);
    neighbors = trainY(CIDX,:);
    stems=neighbors(neighbors(:,3)==1,:);
    nostems=neighbors(neighbors(:,3)==0,:);

    numberOfstems = size(stems,1);
    sumScorestem = sum(stems(:,2));

    numberOfnostems = size(nostems,1);
    sumScorenostem = sum(nostems(:,1));

    
    if  sumScorestem > sumScorenostem
        preds(i) = '1';
    end
end


labels=preds;    
[ms, significant, m1, m2, oracle ] = AverageNDCG(testY,preds);
criterion=1-ms;
end


function sim = similarity(trainX,test)
    sim = pdist2(trainX,test,'minkowski',3);
end
