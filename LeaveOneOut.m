Y=[table2array(ScoresFiltered) double(table2array(Label))-1];
X=table2array(ExtractedUnion);
vals = crossval(@crossValFunc,X,Y,'leaveout',1);
mean(vals)
mean(table2array(ScoresFiltered))



%---MY Method
[m n]=size(X);
predictionScores=zeros(m:1);
for i=1:m
   Xtest=X(i,:);
   Xtrain=X([1:i-1,i+1:end],:);
   
   Ytest=Y(i,:);
   Ytrain=Y([1:i-1,i+1:end],:);
   
   [ms, significant, m1, m2, oracle ]  = crossValFunc(Xtrain,Ytrain,Xtest,Ytest);
   predictionScores(i)=ms; 
    
end

mean(predictionScores)

[h,p] = ttest(predictionScores',table2array(ScoresFiltered(:,2)),'Alpha',0.05)