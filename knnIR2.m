function [criterion, ms, significant, m1, m2,oracle, labels  ] = knnIR2(trainX,trainY,testX,testY,k,Exp)
%KNNIR Summary of this function goes here
%   Detailed explanation goes here
%k=9;
[Z,mu,sigma] = zscore(trainX);
testX = (testX-mu) ./ sigma;
trainX = normalize(trainX);
[m,n]=size(testX);
[numberOfTrainSample,n]=size(trainX);
preds=categorical(zeros(m,1));
priors=[sum(trainY(:,3)==0)/size(trainY(:,3),1) sum(trainY(:,3)==1)/size(trainY(:,3),1)];
[~,cls] = max(priors);
classNames={'0','1'};
priorLabel = classNames(cls);
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

    numberOfstems = size(stems,1);
    sumScorestem = sum(stems(:,2));

    numberOfnostems = size(nostems,1);
    sumScorenostem = sum(nostems(:,1));
    
    
%     alphaSE=0.2;
%     while abs(numberOfstems - numberOfnostems) == 1
%          if alphaSE > 2 
%              break;
%          end
%          diff=abs(trainY(:,1)-trainY(:,2));
%          trainSetInx = diff < alphaSE*(std(diff)/sqrt(size(diff,1)));
%          trainX=trainX(~trainSetInx,:);
%          trainY=trainY(~trainSetInx,:);
%   
%          
%           [CIDX,dist] = knnsearch(trainX, testX(i,:),'k',k,...
%                        'Distance','minkowski','p',Exp);
%     
%         neighbors = trainY(CIDX,:);
%         stems=neighbors(neighbors(:,3)==1,:);
%         nostems=neighbors(neighbors(:,3)==0,:);
% 
%         numberOfstems = size(stems,1);
%         sumScorestem = sum(stems(:,2));
% 
%         numberOfnostems = size(nostems,1);
%         sumScorenostem = sum(nostems(:,1));
%         alphaSE=alphaSE+0.20;
%         
%     end
    
    
    if numberOfstems > numberOfnostems
        preds(i)='1';
    elseif numberOfstems == numberOfnostems     
        preds(i) = priorLabel;
    end
end


labels=preds;    
[ms, significant, m1, m2, oracle ] = AverageNDCG(testY,preds);
criterion=1-ms;
end


function sim = similarity(trainX,test)
    sim = pdist2(trainX,test,'minkowski',3);
end
