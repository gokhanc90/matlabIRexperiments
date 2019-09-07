trainY=[table2array(ScoresFiltered) double(table2array(Label))-1];
trainX=table2array(SelectedFeatures);

[numSample, featureSize]=size(trainX);
ff=zeros(1,46);

keepin=[45];

opts = statset('display','iter','TolFun',0.0);
for i=1:10
    [inmodel,history] = sequentialfs(@criteriaFun,trainX,trainY,'cv',5,'keepin',keepin,'direction','forward','options',opts)
    ff=ff+inmodel;
end

 fb=zeros(1,46);
 for i=1:10
     [inmodel,history] = sequentialfs(@criteriaFun,trainX,trainY,'cv',5,'keepin',keepin,'direction','backward','options',opts)
     fb=fb+inmodel;
 end