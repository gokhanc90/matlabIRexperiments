Joined = join(FeatureCW09BKStem,CW09B_NDCG20_BM25,'LeftKeys',1,'RightKeys',1);
[p,isSig,oracle,label]=getOracle(Joined.NoStem,Joined.KStem);
[idx,C] = kmeans(clstr(:,2:end),2);
clustered=[idx clstr];
qidLabel=[Joined.QueryID label];
clustered=[zeros(size(clustered,1),1) clustered];
[m,n] = size(clustered);
for i=1:m
    labelind=find(qidLabel(:,1)==clustered(i,3),1);
    label=qidLabel(labelind,2);
    clustered(i,1)=label;
end
s=clustered(clustered(:,1)==1,:);
no=clustered(clustered(:,1)==0,:);