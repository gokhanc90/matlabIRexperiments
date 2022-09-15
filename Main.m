
% load('CW09BExpData.mat');
% load('CW12BExpData.mat');
% load('NTCIRExpData.mat');
% load('GOV2ExpData.mat');
% load('WSJExpData.mat');
% load('MQ07ExpData.mat');
% load('MQ08ExpData.mat');
 load('MQ09ExpData.mat');

STEMMERS={'KStem'  'Lovins' 'SnowballEng'};
COLLECTIONS={'CW09B' 'CW12B' 'NTCIR' 'GOV2'  'WSJ' 'MQ07' 'MQ08' 'MQ09'};
%RISK GRAPH
gca=figure();
t = tiledlayout(2,3,'TileSpacing','none','Padding','compact');
%RISK GRAPH END

% Change here
coll = 'MQ09';
docCount = MQ09CollectionStats.(strcat(coll,'DocCount'));
TF = MQ09CollectionStats.(strcat(coll,'TermCount'));
StemmerScores = MQ09Scores;
% Change here END

for  s = 1:size(STEMMERS,2)
        
        stemmer = STEMMERS{s};


        [p,isSig,oracle,label]=getOracle(StemmerScores.NoStem,StemmerScores.(STEMMERS{s}));
        Scores=StemmerScores(:,{'NoStem',STEMMERS{s}});


        functions={@knnIR2};

        fileID = fopen('Result.txt','a');


        Y=[table2array(Scores) label];

        
        Features = strcat(coll,STEMMERS{s},'Features');
        Features = eval(Features);
                
        [m, n]=size(Features);



        X=table2array(Features);

        for K = 1 : length(functions)
            predictionScores=zeros(m,1);
            predictedlabel=categorical(zeros(m,1));
            for i=1:m
                Xtest=X(i,:);
                Xtrain=X([1:i-1,i+1:end],:);

                Ytest=Y(i,:);
                Ytrain=Y([1:i-1,i+1:end],:);

               trainSetInx = Ytrain(:,1) ~= Ytrain(:,2) ; %Discard All same all zero 
               Xtrain=Xtrain(trainSetInx,:);
               Ytrain=Ytrain(trainSetInx,:);

                [criterion, ms, significant, m1, m2, oracle,labels ]  = functions{K}(Xtrain,Ytrain,Xtest,Ytest);
                predictionScores(i)=ms; 
                predictedlabel(i)=labels;
            end

            % Accuracy
            o = categorical(Y(:,3));
            diffInx = Y(:,1) ~= Y(:,2);
            predictedD=predictedlabel(diffInx,:);
            oD=o(diffInx,:);
            Tie=sum(~diffInx)
            TP = sum(predictedD==oD);
           
            TPNo = sum(predictedD==oD & oD=='0')
            TPS = sum(predictedD==oD & oD=='1')


            NoAct=sum('0'==oD)

            SAct=sum('1'==oD)
            accuarcy = (TP/(NoAct+SAct))*100
            [ms, significant, m1, m2, oracle,p, bestSingle, sellArr ] = AverageNDCG(Y(:,[1 2]),predictedlabel);
            fprintf(fileID,'MLFunc: %s Mean: %f Sig: %d NoStemMean: %f StemMean: %f Oracle: %f %0.2f %s %s %s\n',func2str(functions{K}),...
                ms,significant,m1,m2,oracle,p,coll,STEMMERS{s},strjoin(Features.Properties.VariableNames));

            TrisklistSellStem=[sellArr';Y(:,2)'];
            TrisklistSellNoStem=[sellArr';Y(:,1)'];
            TrisklistNoStemStem=[Y(:,1)';Y(:,2)'];
            [h,p]=ttest(sellArr,Y(:,1),'Alpha',0.05);
            [h,p]=ttest(sellArr,Y(:,2),'Alpha',0.05);


            TriskNoSvsSell_0=TRisk(Y(:,1),sellArr,0)
            TriskSvsSell_0=TRisk(Y(:,2),sellArr,0)
            TriskNovsStem_0=TRisk(Y(:,1),Y(:,2),0)

            [scores] = riskscore([Y(:,1)'; sellArr'],[0,1,5],{'NoStem','Sel'},'measure','trisk','baseline',[1])
            [scores] = riskscore([Y(:,2)'; sellArr'],[0,1,5],{'Stem','Sel'},'measure','trisk','baseline',[1])

            [scores] = riskscore([Y(:,1)'; sellArr'],[0,1,5],{'NoStem','Sel'},'measure','grisk','baseline',[1])
            [scores] = riskscore([Y(:,2)'; sellArr'],[0,1,5],{'Stem','Sel'},'measure','grisk','baseline',[1])

            fprintf(fileID,'alpha0_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_0,TriskSvsSell_0,TriskNovsStem_0);

            TriskNoSvsSell_1=TRisk(Y(:,1),sellArr,1)
            TriskSvsSell_1=TRisk(Y(:,2),sellArr,1)
            TriskNovsStem_1=TRisk(Y(:,1),Y(:,2),1)
            fprintf(fileID,'alpha1_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_1,TriskSvsSell_1,TriskNovsStem_1);

            TriskNoSvsSell_2=TRisk(Y(:,1),sellArr,2)
            TriskSvsSell_2=TRisk(Y(:,2),sellArr,2)
            TriskNovsStem_2=TRisk(Y(:,1),Y(:,2),2)
            fprintf(fileID,'alpha2_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_2,TriskSvsSell_2,TriskNovsStem_2);

            TriskNoSvsSell_3=TRisk(Y(:,1),sellArr,3)
            TriskSvsSell_3=TRisk(Y(:,2),sellArr,3)
            TriskNovsStem_3=TRisk(Y(:,1),Y(:,2),3)
            fprintf(fileID,'alpha3_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_3,TriskSvsSell_3,TriskNovsStem_3);

            TriskNoSvsSell_4=TRisk(Y(:,1),sellArr,4)
            TriskSvsSell_4=TRisk(Y(:,2),sellArr,4)
            TriskNovsStem_4=TRisk(Y(:,1),Y(:,2),4)
            fprintf(fileID,'alpha4_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_4,TriskSvsSell_4,TriskNovsStem_4);

            TriskNoSvsSell_5=TRisk(Y(:,1),sellArr,5)
            TriskSvsSell_5=TRisk(Y(:,2),sellArr,5)
            TriskNovsStem_5=TRisk(Y(:,1),Y(:,2),5)
            fprintf(fileID,'alpha5_N_S_StNo:\t%0.4f\t%0.4f\t%0.4f \n',TriskNoSvsSell_5,TriskSvsSell_5,TriskNovsStem_5);


            % RISK GRAPH


            %Stem vs NoStem

            diff = Y(:,2)-Y(:,1);

            sortedDiff = sort(diff);

            numberOfStemGreater = sum(sortedDiff>0);
            perStem = numberOfStemGreater*100/size(sortedDiff,1);
            perStem=round(perStem);
            numberOfNoStemGreater = sum(sortedDiff<0);
            perNoStem = numberOfNoStemGreater*100/size(sortedDiff,1);
            perNoStem = round(perNoStem);

            perTie= 100-perNoStem-perStem;

            nexttile(s+3)
            ax=bar(sortedDiff);
            ylim([-0.8 0.8])
            yticks(-0.8:0.2:0.8)

            text(0.38,0.85,strcat(coll,'-',STEMMERS{s}),'Units','normalized','FontSize',13)

            text(0.45,0.60,[num2str(perTie),'%'],'Units','normalized','FontSize',13)

            text(0.05,0.25,'NoStem > Stem','Units','normalized','FontSize',13)
            text(0.05,0.60,[num2str(perNoStem),'%'],'Units','normalized','FontSize',13)

            text(0.75,0.78,'Stem > NoStem','Units','normalized','FontSize',13)
            text(0.9,0.40,[num2str(perStem),'%'],'Units','normalized','FontSize',13)

            % Sel vs NoStem
            diff = sellArr-Y(:,1);
            sortedDiff = sort(diff);

            numberOfStemGreater = sum(sortedDiff>0);
            perStem = numberOfStemGreater*100/size(sortedDiff,1);
            perStem=round(perStem);
            numberOfNoStemGreater = sum(sortedDiff<0);
            perNoStem = numberOfNoStemGreater*100/size(sortedDiff,1);
            perNoStem = round(perNoStem);

            perTie= 100-perNoStem-perStem;

            nexttile
            ax=bar(sortedDiff);
            ylim([-0.8 0.8])
            yticks(-0.8:0.2:0.8)

            text(0.38,0.85,strcat(coll,'-Sel',STEMMERS{s}),'Units','normalized','FontSize',13)

            text(0.45,0.60,[num2str(perTie),'%'],'Units','normalized','FontSize',13)

            text(0.05,0.25,'NoStem > Sel','Units','normalized','FontSize',13)
            text(0.05,0.60,[num2str(perNoStem),'%'],'Units','normalized','FontSize',13)

            text(0.75,0.78,'Sel > NoStem','Units','normalized','FontSize',13)
            text(0.9,0.40,[num2str(perStem),'%'],'Units','normalized','FontSize',13)

        end
                
          
   xlabel(t,'Number of Queries','FontSize',15)
   ylabel(t,'Diff. in nDCG@20','FontSize',15)
   set(findall(gca,'-property','FontSize'),'FontSize',14)
end



