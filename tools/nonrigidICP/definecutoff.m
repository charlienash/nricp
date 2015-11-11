function [aver, stdevui] = definecutoff( vold, fold )

fk1 = fold(:,1);
fk2 = fold(:,2);
fk3 = fold(:,3);

numverts = size(vold,1);
numfaces = size(fold,1);

D1=sqrt(sum((vold(fk1,:)-vold(fk2,:)).^2,2));
D2=sqrt(sum((vold(fk1,:)-vold(fk3,:)).^2,2));
D3=sqrt(sum((vold(fk2,:)-vold(fk3,:)).^2,2));

aver=mean([D1; D2; D3]);
stdevui=std([D1; D2]);

