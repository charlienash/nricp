function [error,Reallignedsource]=ICPmanu_allign2(target,source,Indices_edgesS,Indices_edgesT)

[IDX1(:,1),IDX1(:,2)]=knnsearch(target,source);
[IDX2(:,1),IDX2(:,2)]=knnsearch(source,target);
IDX1(:,3)=1:length(source(:,1));
IDX2(:,3)=1:length(target(:,1));

[C,ia]=setdiff(IDX2(:,1),Indices_edgesS);
IDX2=IDX2(ia,:);

[C,ia]=setdiff(IDX1(:,1),Indices_edgesT);
IDX1=IDX1(ia,:);

m1=mean(IDX1(:,2));
s1=std(IDX1(:,2));
IDX2=IDX2(IDX2(:,2)<(m1+1.96*s1),:);

Datasetsource=vertcat(source(IDX1(:,3),:),source(IDX2(:,1),:));
Datasettarget=vertcat(target(IDX1(:,1),:),target(IDX2(:,3),:));

[error,Reallignedsource,transform] = procrustes(Datasettarget,Datasetsource);
Reallignedsource=transform.b*source*transform.T+repmat(transform.c(1,1:3),size(source,1),1);


