function[Indices_edges]=detectedges(V,F)

fk1 = F(:,1);
fk2 = F(:,2);
fk3 = F(:,3);

ed1=sort([fk1 fk2 ]')';
ed2=sort([fk1 fk3 ]')';
ed3=sort([fk2 fk3 ]')';

%single edges
ed=[ed1 ;ed2 ;ed3];
[etemp1,ia,ic]=unique(ed,'rows','stable');
esingle=ed(ia,:);

%dubbles
edouble=removerows(ed,ia);

C = setdiff(esingle,edouble,'rows');

Indices_edges=reshape(C,size(C,1)*2,1);

