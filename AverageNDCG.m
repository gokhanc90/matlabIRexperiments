% First column NoStem, second is Stem
% 1 means stem, 0 means noStem
function [ms, significant, m1, m2, m3, oracle ] = AverageNDCG(Scores,selection)
    [m,n]=size(Scores);
    sellArr=zeros(m,1);
    for i=1:m
        if selection(i,1)=='0'
            sellArr(i,1)=Scores(i,1);
        elseif selection(i,1)=='1'
            sellArr(i,1)=Scores(i,2);
        else
            sellArr(i,1)=Scores(i,3);
        end
    end
    
    oracleArr=zeros(m,1);
    for i=1:m
        if Scores(i,1)>=Scores(i,2) && Scores(i,1)>=Scores(i,3)
            oracleArr(i,1)=Scores(i,1);
        elseif Scores(i,2)>=Scores(i,3)
            oracleArr(i,1)=Scores(i,2);
        else
            oracleArr(i,1)=Scores(i,3);
        end
    end
    
    ms=mean(sellArr);
    m1=mean(Scores(:,1));
    m2=mean(Scores(:,2));
    m3=mean(Scores(:,3));
    oracle=mean(oracleArr);
    
    significant=0;
    [h1,p] = ttest(sellArr,Scores(:,1),'Alpha',0.05);
    
    
    [h2,p] = ttest(sellArr,Scores(:,2),'Alpha',0.05);
    
    
    
    [h3,p] = ttest(sellArr,Scores(:,3),'Alpha',0.05);
    
    if h1==1 && h2==1 && h3==1
        significant=1;
    end
    
end