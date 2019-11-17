ids=unique(KStem(:,1));
chiVals=zeros(height(ids),3);
for i=1:height(ids)
   %obs1=NoStem(NoStem.QueryID==table2array(ids(i,1)),{'Advance',});
   obs2=KStem(KStem.QueryID==table2array(ids(i,1)),{'BM25Adv',});
   obs3=Snowball(Snowball.QueryID==table2array(ids(i,1)),{'BM25Adv',});
   
   %obs1=[obs1.BM25Coll ];
   obs2=[obs2.BM25Adv ];
   obs3=[obs3.BM25Adv ];
   
   %obss=[obs1';obs2';obs3'];
   %chi2=chi2Val(obss);
   
   chiVals(i,:)=[table2array(ids(i,1))  mean(obs2) mean(obs3)];
end