STEMMERS={'SnowballEng' 'KStem'};
TWS={'BM25' 'DFIC' 'DFRee' 'DLH13' 'DLM' 'DPH' 'LGD' 'PL2'};
MEASURES={'MAP' 'NDCG20' 'NDCG100'};
COLLECTIONS={'CW12B' 'CW09B' 'NTCIR' 'GOV2' 'WSJ' 'MQ07' 'MQ08' 'MQ09'};

for coll = 1:size(COLLECTIONS,2)
   for tw = 1:size(TWS,2)
        for measure = 1:size(MEASURES,2)
            for s = 1:size(STEMMERS,2)
                
                
                
                    dataNameR = strcat(COLLECTIONS{coll},'_',MEASURES{measure},'_',TWS{tw});
                    runtopic = eval(dataNameR);

                    dataNameF = strcat('Feature',COLLECTIONS{coll},STEMMERS{s});
                    features = eval(dataNameF);


                    %-----INIT----------------
                    %{'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF',...
                    %'MeanCTI','VarCTI','MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
                    %'SCCSNoStem','MeanSCCQNoStem','VarSCCQNoStem',...
                    %'MeanCommonalityKStem','VarCommonalityKStem','SCCSKStem','MeanSCCQKStem','VarSCCQKStem',...
                    %'MeanCommonalitySnowball','VarCommonalitySnowball','SCCSSnowball','MeanSCCQSnowball','VarSCCQSnowball',...
                    %'Chi2DFTF','Chi2IdfIctf','Chi2SCQ','Chi2Commonalities','Chi2SCCQ',...
                    %'MeanSCCS','HarmMeanSCCS','MeanVarCommonality','MeanVarSCCQ',...
                    %'BM25CollNoStem','BM25CollKStem','BM25CollSnowball',...
                    %'AdvanceKStem','AdvanceSnowball','BM25AdvKStem','BM25AdvSnowball'...
                    %}
                    Joined = join(features,runtopic,'LeftKeys',1,'RightKeys',1);
                    % Filtered=MQ07TypeQ(MQ07TypeQ.AllSameAllZero == '0',:);
                     SelectedFeatures=Joined(:,{... 
                    'Gamma','Omega','AvgPMI','SCS','MeanICTF','VarICTF','MeanIDF','VarIDF','MeanCTI','VarCTI',...
                    'MeanSkew','VarSkew','MeanKurt','VarKurt','MeanSCQ','VarSCQ',...
                    'MeanCommonality','VarCommonality','SCCS','MeanSCCQ','VarSCCQ','Advance',...
                    TWS{tw}
                    %'DLM',...
                    %'LGD',...
                    %'PL2',...
                    %'DFIC',...
                    %'DPH',...
                    %'DLH13',...
                    %'DFRee',...
                     });
                     SelectedFeatures=fillmissing(SelectedFeatures,'constant',0);
                    SelectedFeatures=SelectedFeatures(:,Sf);
                    [p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.(STEMMERS{s}));
                    Scores=Joined(:,{'NoStem',STEMMERS{s}});

                    % 
                    %  [idx,weights] = relieff(table2array(SelectedFeatures),table2array(Label),20);
                    %   poz=weights>0
                    %   pf=idx(poz)
                    %   vn=SelectedFeatures.Properties.VariableNames(pf)
                    %   SelectedFeatures = SelectedFeatures(:,vn);
                    %  
                    %---------------------------------------
                    %runtopic=Filtered(:,5:8);

                    option=2; %0:combination 1:remove add else: all

                    fileID = fopen('runtopic.txt','a');
                    % functions={@criteriaFunCoarseKNN,@criteriaFunCubicKNN ,@criteriaFunDiscriminateQuadratic ,@criteriaFunEnsembleRUSBoost ,...
                    % @criteriaFunEnsembleSubspaceDiscriminant ,@criteriaFunEnsembleSubspaceKNN ,@criteriaFunFineTree ,@criteriaFunGaussianNaiveBayes ,...
                    % 	@criteriaFunMediumKNN ,@criteriaFunSVM };

                    functions={@criteriaFunGaussianNaiveBayes };


                    Y=[table2array(Scores) label];

                    [m, n]=size(SelectedFeatures);

                    if option==0
                        for S =[2:4,34:n-1]
                            C=nchoosek(SelectedFeatures.Properties.VariableNames,uint16(S));
                            for Ci = 1:length(C)
                            %X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
                            X=table2array(SelectedFeatures(:,C(Ci,:)));
                            for K = 1 : length(functions)

                                predictionScores=zeros(m,1);
                                predictedlabel=categorical(zeros(m,1));
                                for i=1:m
                                    Xtest=X(i,:);
                                    Xtrain=X([1:i-1,i+1:end],:);

                                    Ytest=Y(i,:);
                                    Ytrain=Y([1:i-1,i+1:end],:);

                                    [criterion, ms, significant, m1, m2,m3, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                    predictionScores(i)=ms; 
                                    predictedlabel(i)=labels;
                                end

                                [ms, significant, m1, m2, oracle ] = AverageNDCG(table2array(Scores),predictedlabel);
                                fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f FeatureSize: %f Feature: %s\n',func2str(functions{K}),...
                                    ms,significant,m1,m2,oracle, S, string(strjoin(C(Ci,:))));
                            end
                           %    runtopic(:,S+4)=table(predictionScores);
                           %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};
                            end
                        end
                    elseif  option==1   
                        for S =1:n
                            %X=table2array(SelectedFeatures(:,S));
                            X=table2array(SelectedFeatures(:,[1:S-1,S+1:end]));
                            for K = 1 : length(functions)

                                predictionScores=zeros(m,1);
                                predictedlabel=categorical(zeros(m,1));
                                for i=1:m
                                    Xtest=X(i,:);
                                    Xtrain=X([1:i-1,i+1:end],:);

                                    Ytest=Y(i,:);
                                    Ytrain=Y([1:i-1,i+1:end],:);

                                    [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                    predictionScores(i)=ms; 
                                    predictedlabel(i)=labels;
                                end

                                [ms, significant, m1, m2, oracle ] = AverageNDCG(table2array(Scores),predictedlabel);
                                fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f Discard: %s\n',func2str(functions{K}),...
                                    ms,significant,m1,m2,oracle, SelectedFeatures.Properties.VariableNames{S});
                            end
                           %    runtopic(:,S+4)=table(predictionScores);
                           %    runtopic.Properties.VariableNames{S+4}=SelectedFeatures.Properties.VariableNames{S};

                        end
                    else
                        X=table2array(SelectedFeatures);
                        %mdl = rica(X,10,'IterationLimit',1e3,'Standardize',1);
                        %X = transform(mdl,X);
                        %X=normalize(X);
                        for K = 1 : length(functions)

                            predictionScores=zeros(m,1);
                            predictedlabel=categorical(zeros(m,1));
                            for i=1:m
                                Xtest=X(i,:);
                                Xtrain=X([1:i-1,i+1:end],:);

                                Ytest=Y(i,:);
                                Ytrain=Y([1:i-1,i+1:end],:);

                                diff=abs(Ytrain(:,1)-Ytrain(:,2));
                                trainSetInx = diff > std(diff) / sqrt(size(diff,1));
                                
                                Xtrain=Xtrain(trainSetInx,:);
                                Ytrain=Ytrain(trainSetInx,:);

                                [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                                predictionScores(i)=ms; 
                                predictedlabel(i)=labels;
                            end

                            [ms, significant, m1, m2, oracle ] = AverageNDCG(Y(:,[1 2]),predictedlabel);
                            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f  %s %s %s\n',func2str(functions{K}),...
                                ms,significant,m1,m2,oracle,dataNameR,dataNameF,strjoin(SelectedFeatures.Properties.VariableNames));

                       %     runtopic(:,end)=table(predictionScores);
                        %    runtopic.Properties.VariableNames{end}='All';

                        end
                    end
                
                
                
                
            end
        end
   end
end



