function [criterion, ms, significant, m1, m2,oracle, labels  ] = knnIR(trainX,trainY,testX,testY,k,Exp)
%KNNIR Summary of this function goes here
%   Detailed explanation goes here
%k=9;
[Z,mu,sigma] = zscore(trainX);
testX = (testX-mu) ./ sigma;
trainX = normalize(trainX);
[m,n]=size(testX);
[numberOfTrainSample,n]=size(trainX);
preds=categorical(zeros(m,1));
for i=1:m
 %   sim = similarity(trainX,testX(i,:));
   
 %   simLabel=[sim trainY];
    
    [CIDX,dist] = knnsearch(trainX, testX(i,:),'k',k,...
                       'Distance','minkowski','p',Exp);
    
 %   sorted = sortrows(simLabel,1,'ascend');
 %   neighbors = sorted(1:k,:);
    neighbors = trainY(CIDX,:);
    stems=neighbors(neighbors(:,3)==1,:);
    nostems=neighbors(neighbors(:,3)==0,:);

    numberOfstemsN = size(stems,1);
    sumScorestemN = sum(stems(:,2));

    numberOfnostemsN = size(nostems,1);
    sumScorenostemN = sum(nostems(:,1));

    if abs(numberOfstemsN - numberOfnostemsN) > k*0.15
        if  numberOfstemsN > numberOfnostemsN
            preds(i) = '1';
        end
    else
        % Look farest neighbor
        [CIDX,dist] = knnsearch(trainX, testX(i,:),'k',numberOfTrainSample,...
                       'Distance','minkowski','p',Exp);
    
         %   sorted = sortrows(simLabel,1,'ascend');
         %   neighbors = sorted(end-k:end,:);
         CIDX = CIDX(length(CIDX)-k+1:length(CIDX));
        neighbors = trainY(CIDX,:);
        stems=neighbors(neighbors(:,3)==1,:);
        nostems=neighbors(neighbors(:,3)==0,:);

        numberOfstemsF = size(stems,1);
        sumScorestemF = sum(stems(:,2));

        numberOfnostemsF = size(nostems,1);
        sumScorenostemF = sum(nostems(:,1));
        
        % reverse assignment
         if  abs(numberOfstemsF - numberOfnostemsF) > k*0.15
            if  numberOfnostemsF > numberOfstemsF
                preds(i) = '1';
            end
         else
            %Uncertain situation
            if sumScorestemN/numberOfstemsN  > sumScorenostemN/numberOfnostemsN 
                preds(i) = '1';
            end
         end
        
    end
end


labels=preds;    
[ms, significant, m1, m2, oracle ] = AverageNDCG(testY,preds);
criterion=1-ms;
end


function sim = similarity(trainX,test)
    sim = pdist2(trainX,test,'minkowski',3);
end
