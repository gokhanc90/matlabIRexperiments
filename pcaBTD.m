function [v,vtilde,d,z,ztilde,ssexp,map,mapsorted,mapranks,rad,radsorted,radranks,corrxz]=pca(x,istilde)
[n,p]=size(x);
r=rank(x);

x_=x;
if (nargin > 1) && (istilde == 1)
    x_=(x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);
end;

[v,d]=eig(x_'*x_);
v=fliplr(v);
%d=rot90(d,2);
d=fliplr(d);
d=flipud(d);

%standardized eigenvector for bi-plots
vtilde=(v-repmat(mean(v),size(v,1),1))./repmat(std(v),size(v,1),1);

%adjust eigenvalues
temp=diag(d);
d=diag([temp(1:r);zeros(p-r,1)]);

%principal component scores
z=x*v;
if (r < n)
    z(:,r+1:end)=0;
end;

%standardized eigenvector for bi-plots
ztilde=(z-repmat(mean(z),size(z,1),1))./repmat(std(z),size(z,1),1);

%sum of squares explained
ssexp=diag(d)./sum(diag(d));

%MAP distances
map=mean(x,2)';

[mapsorted,mapranks]=sort(map);
mapsorted=flipud(mapsorted);
mapranks=flipud(mapranks);

%relative analytic distance RAD
rad=sqrt(diag(x*(diag(var(x))^(-1))*x'));

[radsorted,radranks]=sort(rad);
radsorted=flipud(radsorted);
radranks=flipud(radranks);

%correlation between prin. comps and original variables
x_tilde=(x-repmat(mean(x),size(x,1),1))./repmat(std(x),size(x,1),1);
z_tilde=(z-repmat(mean(z),size(z,1),1))./repmat(std(z),size(z,1),1);
corrxz=x_tilde'*z_tilde/(size(x,1)-1);

%show component scores plot
plot(ztilde(:,1),ztilde(:,2),'ko');
plot(z(:,1),z(:,2),'+');

%show component scores bi-plot
plot(ztilde(:,1),ztilde(:,2),'ko');
hold on;
plot(vtilde(:,1),vtilde(:,2),'k+');
hold off;

title('PCA Plot');
xlabel(cellstr(strcat('First Principle Axis (', num2str(ssexp(1) * 100, 4), '%)')));
ylabel(cellstr(strcat('Second Principle Axis (', num2str(ssexp(2) * 100, 4), '%)')));
