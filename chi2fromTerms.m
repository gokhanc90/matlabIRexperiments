ids=unique(MQ09NoStemTerm(:,1));
chiVals=zeros(height(ids),2);
for i=1:height(ids)
   obs1=MQ09NoStemTerm(MQ09NoStemTerm.QueryID==table2array(ids(i,1)),{'DF','TF'});
   obs2=MQ09KStemTerm(MQ09KStemTerm.QueryID==table2array(ids(i,1)),{'DF','TF'});
   obs3=MQ09SnowballTerm(MQ09SnowballTerm.QueryID==table2array(ids(i,1)),{'DF','TF'});
   
   obs1=[obs1.DF;obs1.TF];
   obs2=[obs2.DF;obs2.TF];
   obs3=[obs3.DF;obs3.TF];
   
   obss=[obs1';obs2';obs3'];
   chi2=chi2Val(obss);
   
   chiVals(i,:)=[table2array(ids(i,1)) chi2];
end