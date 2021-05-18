% Generate Oracle and test whether it (ORACLE) is statistical significant to the
% both given systems.
%First one is nostem
function [p,isSig,oracle,label] = getOracle(l1,l2)
    [m,n]=size(l1);
    oracle=zeros(m,1);
    label=zeros(m,1);
    for i=1:m
        if l1(i,1)>l2(i,1)
            oracle(i,1)=l1(i,1);
        elseif l1(i,1) < l2(i,1)
            oracle(i,1)=l2(i,1);
            label(i,1)=1;
        else
            oracle(i,1)=l1(i,1);
            label(i,1)=2;
        end
    end

    isSig=0;
    [h1,p1] = ttest(oracle,l1,'Alpha',0.05);
    
    
    [h2,p2] = ttest(oracle,l2,'Alpha',0.05);
    
    p=p1;
    if p2>p1 
        p=p2;
    end
    
    if h1==1 && h2==1
        isSig=1;
    end
end

