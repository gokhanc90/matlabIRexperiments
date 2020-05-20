% First column NoStem, second is Stem
% 1 means stem, 0 means noStem
function [ms, significant, m1, m2, oracle ] = AverageNDCG(Scores,selection)
    [m,n]=size(Scores);
    sellArr=zeros(m,1);
    for i=1:m
        if selection(i,1)=='0' ||  selection(i,1) == 0
            sellArr(i,1)=Scores(i,1);
        else
            sellArr(i,1)=Scores(i,2);
        end
    end
    
    oracleArr=zeros(m,1);
    for i=1:m
        if Scores(i,1)>=Scores(i,2) 
            oracleArr(i,1)=Scores(i,1);
        else
            oracleArr(i,1)=Scores(i,2);
        end
    end
    
    ms=mean(sellArr);
    m1=mean(Scores(:,1));
    m2=mean(Scores(:,2));
    oracle=mean(oracleArr);
    
    significant=0;
    [h1,p] = ttest(sellArr,Scores(:,1),'Alpha',0.05);
    
    
    [h2,p] = ttest(sellArr,Scores(:,2),'Alpha',0.05);
    
    
    if  ms>m1 && ms>m2
        significant=2;
    end
    
    if h1==1 && h2==1 && ms>m1 && ms>m2
        significant=1;
    end
    
end