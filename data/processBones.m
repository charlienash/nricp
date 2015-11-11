%% Process bones
cd /home/charlie/nricp

%% Source
load data/EXAMPLE1.mat
Source.faces = sourceF;
Source.vertices = sourceV;

save('data/boneSource.mat', 'Source');

%% Target
Target.faces = targetF;
Target.vertices = targetV;

save('data/boneTarget.mat', 'Target');