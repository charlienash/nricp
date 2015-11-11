%% Process Faces
cd /home/charlie/nricp

%% Source
load data/face02.mat
Source.faces = surface.TRIV;
Source.vertices = [surface.X surface.Y surface.Z];
Source.normals = compute_normal(Source.vertices, Source.faces)';

save('data/faceSource.mat', 'Source');

%% Target
load data/face03.mat
Target.faces = surface.TRIV;
Target.vertices = [surface.X surface.Y surface.Z];
Target.normals = compute_normal(Target.vertices, Target.faces)';

save('data/faceTarget.mat', 'Target');