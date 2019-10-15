ids=unique(NoStemTerm(:,1));
chiVals=zeros(height(ids),2);
for i=1:height(ids)
   obs1=NoStemTerm(NoStemTerm.QueryID==table2array(ids(i,1)),{'sccq'});
   obs2=KStemTerm(KStemTerm.QueryID==table2array(ids(i,1)),{'sccq'});
   obs3=SnowballTerm(SnowballTerm.QueryID==table2array(ids(i,1)),{'sccq'});
   
   obs1=[obs1.sccq];
   obs2=[obs2.sccq];
   obs3=[obs3.sccq];
   
   obss=[obs1';obs2';obs3'];
   chi2=chi2Val(obss);
   
   chiVals(i,:)=[table2array(ids(i,1)) chi2];
end