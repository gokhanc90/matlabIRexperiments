trainY=[table2array(ScoresFiltered) double(table2array(Label))-1];
trainX=table2array(SelectedFeatures);
[numSample, featureSize]=size(trainX);

%
functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunDiscriminateQuadratic ,@criteriaFunEnsembleRUSBoost ,...
	@criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
	@criteriaFunMediumKNN ,@criteriaFunSVM };

keepin=[43 44 45 46];
for K = 1 : length(functions)
	for i = 1 : length(keepin)
        keep=keepin(i);   	

        ff=zeros(1,46);

        opts = statset('display','iter','TolFun',0.0);
        for i=1:10
           [inmodel,history] = sequentialfs(functions{K},trainX,trainY,'cv',5,'keepin',keep,'direction','forward','options',opts)
          ff=ff+inmodel;
        end

        fb=zeros(1,46);
        for i=1:10
           [inmodel,history] = sequentialfs(functions{K},trainX,trainY,'cv',5,'keepin',keep,'direction','backward','options',opts)
           fb=fb+inmodel;
        end

        fileID = fopen(strcat('ff',func2str(functions{K}),num2str(keep),'.txt'),'w');
        fprintf(fileID,'%i\t',ff);
        fclose(fileID);

        fileID = fopen(strcat('fb',func2str(functions{K}),num2str(keep),'.txt'),'w');
        fprintf(fileID,'%i\t',fb);
        fclose(fileID);
	end
end




