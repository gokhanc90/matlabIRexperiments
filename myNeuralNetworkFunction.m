function [criterion, ms, significant, m1, m2,oracle, labels  ]  = myNeuralNetworkFunction(trainX,trainY,testX,testY)
[Z,mu,sigma] = zscore(trainX);
testX = (testX-mu) ./ sigma;
trainX = normalize(trainX);

[m,n]=size(testX);
preds=categorical(zeros(m,1));

for i=1:m
    net = patternnet(10);
    net.trainParam.showWindow = 0;
    net = train(net,trainX',[ ~trainY(:,3) trainY(:,3)]');
    y = net(testX');
    class = vec2ind(y);
    if class==2
        preds(i)='1';
    end
end

labels=preds;    
[ms, significant, m1, m2, oracle ] = AverageNDCG(testY,preds);
criterion=1-ms;


end