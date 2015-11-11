function [error,Reallignedsource,transform]=rigidICP(target,source,flag,Indices_edgesS,Indices_edgesT)

% This function rotates, translates and scales a 3D pointcloud "source" of N*3 size (N points in N rows, 3 collumns for XYZ)
% to fit a similar shaped point cloud "target" again of N by 3 size
% 
% The output shows the minimized value of dissimilarity measure in "error", the transformed source data set and the 
% transformation, rotation, scaling and translation in transform.T, transform.b and transform.c such that
% Reallignedsource = b*source*T + c;


if flag==0
[Prealligned_source,Prealligned_target,transformtarget ]=Preall(target,source);
else
    Prealligned_source=source;
    Prealligned_target=target;
end

display ('error')
errortemp(1,:)=0;
index=2;
[errortemp(index,:),Reallignedsourcetemp]=ICPmanu_allign2(Prealligned_target,Prealligned_source,Indices_edgesS,Indices_edgesT);

while (abs(errortemp(index-1,:)-errortemp(index,:)))>0.000001
[errortemp(index+1,:),Reallignedsourcetemp]=ICPmanu_allign2(Prealligned_target,Reallignedsourcetemp,Indices_edgesS,Indices_edgesT);
index=index+1;
d=errortemp(index,:)

end

error=errortemp(index,:);

if flag==0
Reallignedsource=Reallignedsourcetemp*transformtarget.T+repmat(transformtarget.c(1,1:3),length(Reallignedsourcetemp(:,1)),1);
[d,Reallignedsource,transform] = procrustes(Reallignedsource,source);
else
   [d,Reallignedsource,transform] = procrustes(Reallignedsourcetemp,source);
end
