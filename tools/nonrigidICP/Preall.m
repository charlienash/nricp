function [Prealligned_source,Prealligned_target,transformtarget ]=Preall(target,source)

% This function performs a first and rough pre-alligment of the data as starting position for the iterative allignment and scaling procedure

% Initial positioning of the data is based on alligning the coordinates of the objects -which are assumed to be close/similar in shape- following principal component analysis

[COEFF,Prealligned_source] = princomp(source);

[COEFF,Prealligned_target] = princomp(target);

% the direction of the axes is than evaluated and corrected if necesarry.
Maxtarget=max(Prealligned_source)-min(Prealligned_source);
Maxsource=max(Prealligned_target)-min(Prealligned_target);
D=Maxtarget./Maxsource;
D=[D(1,1) 0 0;0 D(1,2) 0; 0 0 D(1,3)];
RTY=Prealligned_source*D;

load R
for i=1:8
    T=R{1,i};
    T=RTY*T;
    [bb DD]=knnsearch(T,Prealligned_target);
    MM(i,1)=sum(DD);
end

[M I]=min(MM);
 T=R{1,I};
 Prealligned_source=Prealligned_source*T;
 
 [d,Z,transformtarget] = procrustes(target,Prealligned_target);