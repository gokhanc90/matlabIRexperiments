% First column NoStem, second is Stem
% 1 means stem, 0 means noStem
function [ms, significant, m1, m2, oracle, p, bestSingle, sellArr ] = AverageNDCG(Scores,selection)
    [m,n]=size(Scores);
    sellArr=zeros(m,1);
    for i=1:m
        if selection(i,1)=='0'
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
    if m1 > m2
        bestSingle = Scores(:,1);
    else
        bestSingle = Scores(:,2);
    end
    
    significant=0;
    [h1,p1] = ttest(sellArr,Scores(:,1),'Alpha',0.05);
    
    
    [h2,p2] = ttest(sellArr,Scores(:,2),'Alpha',0.05);
    
    if p1 > p2
        p = p1;
    else
        p=p2;
    end
    
    if  (ms<m1 && ms<m2) && (h1==1 && h2==1)
        significant=-1;
    end
    
    if  ms<m1 || ms<m2
        significant=0;
    end
    
    
    if  ms>m1 && ms>m2
        significant=1;
    end
    
    if h1==1 && h2==1 && ms>m1 && ms>m2
        significant=2;
    end
    
end